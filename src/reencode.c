/*
 * ices0 - Icecast Source Client
 * 5-band broadcast leveller / equalizer reencoding stage.
 *
 * Signal path:
 *
 *    PCM in ──► [5-band split] ──► [per-band static EQ + gentle compressor]
 *                                          │
 *                                          ▼
 *                                   [sum bands back]
 *                                          │
 *                                          ▼
 *                                 [slow overall AGC]
 *                                          │
 *                                          ▼
 *                                  [hard-knee limiter]
 *                                          │
 *                                          ▼
 *                                       PCM out
 *
 * The per-band stage applies a FIXED equalizer gain (user-configurable
 * below in `band_eq_db[]`) plus gentle 2:1 compression to even out band
 * dynamics. The overall AGC operates much more slowly than the band
 * compressors (~10 s vs ~500 ms) so it doesn't fight the bands or pump.
 * The final limiter catches transients without colouring sustained tones.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <lame/lame.h>

#include <shout/shout.h>

#include "icestypes.h"
#include "ices_config.h"
#include "reencode.h"

#define REENCODE_VERBOSE 1

#if REENCODE_VERBOSE
#  define LOG(...) do { fprintf(stderr, "[reencode] " __VA_ARGS__); \
                        fputc('\n', stderr); fflush(stderr); } while (0)
#else
#  define LOG(...) do {} while (0)
#endif

/* =========================================================================
 * USER-TUNABLE EQ CURVE
 *
 * Static per-band equalizer in dB. Positive = boost, negative = cut.
 * Bands (Hz): [0-100] [100-400] [400-1500] [1500-5000] [5000+]
 * Defaults: small smile curve, gentle bass and treble lift.
 * ========================================================================= */

static double band_eq_db[5] = {
    +2.0,    /* sub/bass     ( <100 Hz)   slight lift */
     0.0,    /* low-mid      (100-400 Hz) flat        */
    -1.0,    /* mid          (400-1.5k)   slight cut for clarity */
     0.0,    /* presence     (1.5k-5k)    flat        */
    +1.5     /* air/treble   ( >5k)       slight lift */
};

/* Overall target output loudness (RMS, as a fraction of full-scale).
 * 0.20 ≈ -14 dBFS RMS, a sensible broadcast loudness with headroom for
 * transients. The slow AGC pulls every track toward this target. */
static const double overall_target_rms = 0.20;

/* Per-band compressor: gentle 2:1 above this threshold, no expansion below. */
static const double band_comp_threshold = 0.10;   /* per-band RMS */
static const double band_comp_ratio     = 2.0;    /* 2:1 */

/* =========================================================================
 * 5-BAND SPLIT (cascaded one-pole; bands sum to input by construction)
 * ========================================================================= */

#define NBANDS 5

typedef struct {
    double p1_0, p1_1, p2_0, p2_1, p3_0, p3_1, p4_0, p4_1;
} FilterState;

static FilterState state_l, state_r;

static double xover[4] = {
    100.0  / 44100.0, 400.0  / 44100.0,
    1500.0 / 44100.0, 5000.0 / 44100.0
};

/* Pre-computed linear EQ gains from band_eq_db[]. */
static double band_eq_lin[5] = { 1.0, 1.0, 1.0, 1.0, 1.0 };

/* Per-band RMS estimators (slowly tracking, used for compression). */
static double band_rms[5] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
static double band_comp_gain[5] = { 1.0, 1.0, 1.0, 1.0, 1.0 };

/* Overall AGC: tracks signal RMS over ~10 s, computes one makeup gain. */
static double overall_rms     = 0.0;
static double overall_gain    = 1.0;

/* Smoothing coefficients (recomputed when sample rate changes). */
static double band_rms_coef   = 1.0 / (0.500 * 44100.0);   /* ~500 ms */
static double overall_rms_coef= 1.0 / (10.0  * 44100.0);   /* ~10 s   */
static double overall_gain_coef= 1.0 / (5.0  * 44100.0);   /* ~5 s    */

/* Limiter state: fast attack, slow release, hard ceiling. */
static double limiter_gain = 1.0;
static const double limiter_ceiling = 0.97;        /* -0.26 dBFS */
static double limiter_atk_coef = 1.0 / (0.001 * 44100.0);  /* ~1 ms attack */
static double limiter_rel_coef = 1.0 / (0.100 * 44100.0);  /* ~100 ms release */

static void recompute_eq(void) {
    for (int b = 0; b < NBANDS; b++)
        band_eq_lin[b] = pow(10.0, band_eq_db[b] / 20.0);
}

static void recompute_coeffs(unsigned int samplerate) {
    if (samplerate < 8000) samplerate = 44100;
    double sr = (double)samplerate;

    xover[0] = 100.0  / sr;
    xover[1] = 400.0  / sr;
    xover[2] = 1500.0 / sr;
    xover[3] = 5000.0 / sr;
    for (int i = 0; i < 4; i++) if (xover[i] > 0.45) xover[i] = 0.45;

    band_rms_coef     = 1.0 / (0.500 * sr);
    overall_rms_coef  = 1.0 / (10.0  * sr);
    overall_gain_coef = 1.0 / (5.0   * sr);
    limiter_atk_coef  = 1.0 / (0.001 * sr);
    limiter_rel_coef  = 1.0 / (0.100 * sr);
}

static void split_5_bands(double input, FilterState *s, double *bands) {
    s->p1_0 += xover[0] * (input    - s->p1_0);
    s->p1_1 += xover[0] * (s->p1_0  - s->p1_1);  double lp1 = s->p1_1;

    s->p2_0 += xover[1] * (input    - s->p2_0);
    s->p2_1 += xover[1] * (s->p2_0  - s->p2_1);  double lp2 = s->p2_1;

    s->p3_0 += xover[2] * (input    - s->p3_0);
    s->p3_1 += xover[2] * (s->p3_0  - s->p3_1);  double lp3 = s->p3_1;

    s->p4_0 += xover[3] * (input    - s->p4_0);
    s->p4_1 += xover[3] * (s->p4_0  - s->p4_1);  double lp4 = s->p4_1;

    bands[0] = lp1;
    bands[1] = lp2 - lp1;
    bands[2] = lp3 - lp2;
    bands[3] = lp4 - lp3;
    bands[4] = input - lp4;
}

/* Apply gentle compression to a band, returning the gain multiplier. */
static double band_compressor(int b, double sample_l, double sample_r) {
    /* Track per-band RMS-ish envelope. */
    double inst = 0.5 * (sample_l * sample_l + sample_r * sample_r);
    band_rms[b] += (inst - band_rms[b]) * band_rms_coef;

    double env = sqrt(band_rms[b]);
    double target_gain = 1.0;

    if (env > band_comp_threshold) {
        /* 2:1 above threshold. Above-threshold-by ratio:
         * desired_env = thresh + (env - thresh) / ratio
         * gain = desired_env / env */
        double over = env - band_comp_threshold;
        double desired = band_comp_threshold + over / band_comp_ratio;
        target_gain = desired / env;
    }

    /* Smooth gain change (~50 ms attack, gentle release). */
    double smoothing = (target_gain < band_comp_gain[b])
                       ? band_rms_coef * 4.0   /* faster attack */
                       : band_rms_coef * 1.0;  /* slower release */
    band_comp_gain[b] += (target_gain - band_comp_gain[b]) * smoothing;
    return band_comp_gain[b];
}

static void process_audio_channels(int16_t *left, int16_t *right, int nsamples) {
    double bands_l[NBANDS], bands_r[NBANDS];

    for (int i = 0; i < nsamples; i++) {
        double sample_l = (double)left[i]  / 32768.0;
        double sample_r = (double)right[i] / 32768.0;

        split_5_bands(sample_l, &state_l, bands_l);
        split_5_bands(sample_r, &state_r, bands_r);

        double mix_l = 0.0, mix_r = 0.0;

        for (int b = 0; b < NBANDS; b++) {
            double comp_gain = band_compressor(b, bands_l[b], bands_r[b]);
            double total_band_gain = band_eq_lin[b] * comp_gain;
            mix_l += bands_l[b] * total_band_gain;
            mix_r += bands_r[b] * total_band_gain;
        }

        /* Slow overall AGC: track summed RMS, nudge gain toward target. */
        double overall_inst = 0.5 * (mix_l * mix_l + mix_r * mix_r);
        overall_rms += (overall_inst - overall_rms) * overall_rms_coef;

        double current_loudness = sqrt(overall_rms);
        if (current_loudness > 0.005) {
            double desired_gain = overall_target_rms / current_loudness;
            /* Soft-limit AGC range to ±12 dB so a near-silent or already-loud
             * track doesn't get mangled. */
            if (desired_gain > 4.0)  desired_gain = 4.0;
            if (desired_gain < 0.25) desired_gain = 0.25;
            overall_gain += (desired_gain - overall_gain) * overall_gain_coef;
        }

        mix_l *= overall_gain;
        mix_r *= overall_gain;

        /* Look-back limiter: when signal exceeds ceiling, drop gain fast;
         * release slowly. This catches transients transparently. */
        double peak = fabs(mix_l) > fabs(mix_r) ? fabs(mix_l) : fabs(mix_r);
        double lim_target = 1.0;
        if (peak * limiter_gain > limiter_ceiling) {
            lim_target = limiter_ceiling / peak;
        }
        double lim_coef = (lim_target < limiter_gain)
                          ? limiter_atk_coef : limiter_rel_coef;
        limiter_gain += (lim_target - limiter_gain) * lim_coef;
        if (limiter_gain > 1.0) limiter_gain = 1.0;
        if (limiter_gain < 0.1) limiter_gain = 0.1;

        mix_l *= limiter_gain;
        mix_r *= limiter_gain;

        /* Hard safety clip (should be vanishingly rare). */
        if (mix_l >  0.999) mix_l =  0.999;
        if (mix_l < -0.999) mix_l = -0.999;
        if (mix_r >  0.999) mix_r =  0.999;
        if (mix_r < -0.999) mix_r = -0.999;

        left[i]  = (int16_t)(mix_l * 32767.0);
        right[i] = (int16_t)(mix_r * 32767.0);
    }
}

static void reset_dsp_state(void) {
    memset(&state_l, 0, sizeof(state_l));
    memset(&state_r, 0, sizeof(state_r));
    for (int b = 0; b < NBANDS; b++) {
        band_rms[b]       = 0.0;
        band_comp_gain[b] = 1.0;
    }
    /* Don't reset overall_gain across tracks — that's the whole point of
     * inter-track levelling. It should slowly carry over and adapt. */
    overall_rms = overall_target_rms * overall_target_rms;  /* sane starting point */
    limiter_gain = 1.0;
}

/* =========================================================================
 * LAME / decoder state  (unchanged from previous working version)
 * ========================================================================= */

static lame_global_flags *gfp         = NULL;
static hip_t              hip_decoder = NULL;
static int cur_in_rate  = 0;
static int cur_in_chans = 0;
static int cur_out_rate  = 0;
static int cur_out_chans = 0;
static int cur_out_brate = 0;
static input_stream_t *pending_source = NULL;

static int lame_configure(int in_rate, int in_chans,
                          int out_rate, int out_chans, int out_brate)
{
    if (gfp) { lame_close(gfp); gfp = NULL; }
    gfp = lame_init();
    if (!gfp) { LOG("lame_init() failed"); return -1; }

    lame_set_in_samplerate (gfp, in_rate);
    lame_set_num_channels  (gfp, in_chans);
    lame_set_out_samplerate(gfp, out_rate > 0 ? out_rate : in_rate);
    lame_set_mode          (gfp, out_chans == 1 ? MONO : JOINT_STEREO);
    lame_set_brate         (gfp, out_brate > 0 ? out_brate : 128);
    lame_set_quality       (gfp, 2);
    lame_set_bWriteVbrTag  (gfp, 0);

    if (lame_init_params(gfp) < 0) {
        LOG("lame_init_params() failed for %dHz/%dch -> %dHz/%dch @ %dkbps",
            in_rate, in_chans, out_rate, out_chans, out_brate);
        lame_close(gfp); gfp = NULL;
        return -1;
    }

    cur_in_rate   = in_rate;
    cur_in_chans  = in_chans;
    cur_out_rate  = out_rate > 0 ? out_rate : in_rate;
    cur_out_chans = out_chans;
    cur_out_brate = out_brate > 0 ? out_brate : 128;

    LOG("LAME configured: in %dHz/%dch -> out %dHz/%dch @ %dkbps",
        cur_in_rate, cur_in_chans, cur_out_rate, cur_out_chans, cur_out_brate);
    return 0;
}

static int absorb_stream_config(ices_stream_t *stream) {
    if (!stream) return 0;
    int new_brate = stream->bitrate         > 0 ? stream->bitrate         : cur_out_brate;
    int new_rate  = stream->out_samplerate  > 0 ? stream->out_samplerate  : cur_out_rate;
    int new_chans = stream->out_numchannels > 0 ? stream->out_numchannels : cur_out_chans;
    if (new_brate == cur_out_brate && new_rate == cur_out_rate
        && new_chans == cur_out_chans)
        return 0;
    cur_out_brate = new_brate;
    cur_out_rate  = new_rate;
    cur_out_chans = new_chans;
    return 1;
}

void ices_reencode_initialize(void) {
    extern ices_config_t ices_config;
    ices_stream_t *stream;

    for (stream = ices_config.streams; stream; stream = stream->next) {
        if (stream->reencode) { ices_config.reencode = 1; break; }
    }
    if (!ices_config.reencode) {
        LOG("no <Stream> has <Reencode>1</Reencode>; reencoder disabled");
        return;
    }

    if (hip_decoder) hip_decode_exit(hip_decoder);
    hip_decoder = hip_decode_init();
    if (!hip_decoder) LOG("hip_decode_init() failed");

    reset_dsp_state();
    recompute_eq();
    recompute_coeffs(44100);

    /* Pretty-print the active EQ curve. */
    LOG("5-band leveller initialised");
    LOG("EQ curve: %.1f / %.1f / %.1f / %.1f / %.1f dB  (sub/lowmid/mid/pres/air)",
        band_eq_db[0], band_eq_db[1], band_eq_db[2], band_eq_db[3], band_eq_db[4]);
    LOG("Target loudness: %.3f RMS (%.1f dBFS), band comp %.1f:1 above %.2f",
        overall_target_rms, 20.0 * log10(overall_target_rms),
        band_comp_ratio, band_comp_threshold);
}

void ices_reencode_reset(input_stream_t* source) {
    reset_dsp_state();

    if (hip_decoder) hip_decode_exit(hip_decoder);
    hip_decoder = hip_decode_init();
    if (!hip_decoder) LOG("hip_decode_init() failed during reset");

    pending_source = source;

    if (source) {
        LOG("reset: track=%s  source=%uHz/%uch @ %ukbps",
            source->path ? source->path : "(null)",
            source->samplerate, source->channels, source->bitrate);
        if (source->samplerate > 0)
            recompute_coeffs(source->samplerate);
    }
}

int ices_reencode_decode(unsigned char* buf, size_t blen, size_t olen,
                         int16_t* left, int16_t* right)
{
    (void)olen;
    if (!hip_decoder) return -1;
    int n = hip_decode(hip_decoder, buf, blen, left, right);
    if (n < 0) return -1;
    return n;
}

int ices_reencode(ices_stream_t* stream, int nsamples,
                  int16_t* left, int16_t* right,
                  unsigned char* outbuf, int outbuf_sz)
{
    if (nsamples <= 0) return 0;

    int out_changed = absorb_stream_config(stream);
    int need_reconfig = 0;
    int in_rate  = cur_in_rate;
    int in_chans = cur_in_chans;

    if (pending_source) {
        int src_rate  = (int)pending_source->samplerate;
        int src_chans = (int)pending_source->channels;
        if (src_rate  <= 0) src_rate  = 44100;
        if (src_chans <= 0) src_chans = 2;
        if (src_rate != cur_in_rate || src_chans != cur_in_chans) {
            in_rate  = src_rate;
            in_chans = src_chans;
            need_reconfig = 1;
        }
        pending_source = NULL;
    }

    if (out_changed || !gfp) {
        if (in_rate  <= 0) in_rate  = 44100;
        if (in_chans <= 0) in_chans = 2;
        need_reconfig = 1;
    }

    if (need_reconfig) {
        int out_rate  = cur_out_rate  > 0 ? cur_out_rate  : in_rate;
        int out_chans = cur_out_chans > 0 ? cur_out_chans : in_chans;
        int out_brate = cur_out_brate > 0 ? cur_out_brate : 128;
        if (lame_configure(in_rate, in_chans, out_rate, out_chans, out_brate) < 0) {
            LOG("LAME reconfigure failed; aborting track");
            return -2;
        }
    }

    if (!gfp) return -2;

    process_audio_channels(left, right, nsamples);

    int n = lame_encode_buffer(gfp, left, right, nsamples, outbuf, outbuf_sz);
    if (n == -1) return -1;
    if (n < -1) {
        LOG("lame_encode_buffer returned %d (samples=%d outbuf_sz=%d)",
            n, nsamples, outbuf_sz);
        return n;
    }
    return n;
}

int ices_reencode_flush(ices_stream_t* stream,
                        unsigned char *outbuf, int outbuf_sz)
{
    (void)stream;
    if (!gfp) return 0;
    int n = lame_encode_flush(gfp, outbuf, outbuf_sz);
    if (n < 0) return 0;
    return n;
}

void ices_reencode_shutdown(void) {
    if (hip_decoder) { hip_decode_exit(hip_decoder); hip_decoder = NULL; }
    if (gfp)         { lame_close(gfp);              gfp         = NULL; }
}
