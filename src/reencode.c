/*
 * ices0 - Icecast Source Client
 * 5-band broadcast leveller / equalizer reencoding stage.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
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

/* Set to 1 if you want a few one-line status messages on stderr per track.
 * Leave it on for the first run after this patch so you can confirm the
 * config it picks up matches your ices.conf. Set to 0 once happy. */
#define REENCODE_VERBOSE 1

#if REENCODE_VERBOSE
#  define LOG(...) do { fprintf(stderr, "[reencode] " __VA_ARGS__); \
                        fputc('\n', stderr); fflush(stderr); } while (0)
#else
#  define LOG(...) do {} while (0)
#endif

/* =========================================================================
 * 5-BAND MULTIBAND LEVELLER
 *
 * Crossover frequencies (Hz): 100, 400, 1500, 5000
 *   band 0:    0 -  100 Hz   (sub / bass)
 *   band 1:  100 -  400 Hz   (low-mid)
 *   band 2:  400 - 1500 Hz   (mid)
 *   band 3: 1500 - 5000 Hz   (high-mid / presence)
 *   band 4: 5000+      Hz   (air / treble)
 *
 * Each band has its own slow AGC-style leveller (~2 s adaptation) so quiet
 * tracks come up and loud tracks come down without pumping.
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

static const double target_energy[NBANDS] = { 0.18, 0.20, 0.22, 0.18, 0.12 };
static double dynamic_gains[NBANDS] = { 1.0, 1.0, 1.0, 1.0, 1.0 };
static double band_envelope[NBANDS] = { 0.0, 0.0, 0.0, 0.0, 0.0 };

static double env_coef  = 1.0 / (0.030 * 44100.0);   /* ~30 ms envelope */
static double gain_coef = 1.0 / (2.0   * 44100.0);   /* ~2  s gain adapt */

static double softclip(double x) {
    if (x >  1.5) x =  1.5;
    if (x < -1.5) x = -1.5;
    double y = x - (x * x * x) / 3.0;
    if (y >  0.99) y =  0.99;
    if (y < -0.99) y = -0.99;
    return y;
}

static void recompute_coeffs(unsigned int samplerate) {
    if (samplerate < 8000) samplerate = 44100;
    double sr = (double)samplerate;

    xover[0] = 100.0  / sr;
    xover[1] = 400.0  / sr;
    xover[2] = 1500.0 / sr;
    xover[3] = 5000.0 / sr;
    for (int i = 0; i < 4; i++) if (xover[i] > 0.45) xover[i] = 0.45;

    env_coef  = 1.0 / (0.030 * sr);
    gain_coef = 1.0 / (2.0   * sr);
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

static void process_audio_channels(int16_t *left, int16_t *right, int nsamples) {
    double bands_l[NBANDS], bands_r[NBANDS];

    for (int i = 0; i < nsamples; i++) {
        double sample_l = (double)left[i]  / 32768.0;
        double sample_r = (double)right[i] / 32768.0;

        split_5_bands(sample_l, &state_l, bands_l);
        split_5_bands(sample_r, &state_r, bands_r);

        double mix_l = 0.0, mix_r = 0.0;

        for (int b = 0; b < NBANDS; b++) {
            double inst = 0.5 * (fabs(bands_l[b]) + fabs(bands_r[b]));
            band_envelope[b] += (inst - band_envelope[b]) * env_coef;

            if (band_envelope[b] > 0.001) {
                double desired = target_energy[b] / band_envelope[b];
                dynamic_gains[b] += (desired - dynamic_gains[b]) * gain_coef;
            }

            if (dynamic_gains[b] > 4.0)  dynamic_gains[b] = 4.0;
            if (dynamic_gains[b] < 0.25) dynamic_gains[b] = 0.25;

            mix_l += bands_l[b] * dynamic_gains[b];
            mix_r += bands_r[b] * dynamic_gains[b];
        }

        mix_l = softclip(mix_l);
        mix_r = softclip(mix_r);

        left[i]  = (int16_t)(mix_l * 32767.0);
        right[i] = (int16_t)(mix_r * 32767.0);
    }
}

/* =========================================================================
 * LAME / decoder state
 * ========================================================================= */

static lame_global_flags *gfp         = NULL;
static hip_t              hip_decoder = NULL;

/* Cached input format (what the current MP3 file is). */
static int cur_in_rate  = 0;
static int cur_in_chans = 0;

/* Cached output format (what ices.conf wants). 0 = "not known yet" — we
 * pick these up from the ices_stream_t* on the first ices_reencode call. */
static int cur_out_rate  = 0;
static int cur_out_chans = 0;
static int cur_out_brate = 0;

/* Has the source for the current track changed since last reset? */
static input_stream_t *pending_source = NULL;

static int lame_configure(int in_rate, int in_chans,
                          int out_rate, int out_chans, int out_brate)
{
    if (gfp) { lame_close(gfp); gfp = NULL; }

    gfp = lame_init();
    if (!gfp) { LOG("lame_init() failed"); return -1; }

    lame_set_in_samplerate (gfp, in_rate);
    /* LAME's num_channels is what stream.c will hand us per call. stream.c
     * always passes two int16 buffers (left+right); for a mono source it
     * duplicates left into right. So we always tell LAME stereo input
     * unless ices.conf explicitly forces mono output AND the source is
     * already mono. Simpler and always correct: in_chans matches the
     * channel layout of the left/right buffers the loop passes us. */
    lame_set_num_channels  (gfp, in_chans);
    lame_set_out_samplerate(gfp, out_rate > 0 ? out_rate : in_rate);
    lame_set_mode          (gfp, out_chans == 1 ? MONO : JOINT_STEREO);
    lame_set_brate         (gfp, out_brate > 0 ? out_brate : 128);
    lame_set_quality       (gfp, 2);
    lame_set_bWriteVbrTag  (gfp, 0);   /* live stream — no Xing header */

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

    LOG("LAME configured: in %dHz/%dch  ->  out %dHz/%dch @ %dkbps",
        cur_in_rate, cur_in_chans, cur_out_rate, cur_out_chans, cur_out_brate);
    return 0;
}

/* Apply output config from the ices_stream_t to our cache. Returns 1 if
 * anything changed, 0 if the cache already matched. */
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

    /* CRITICAL: propagate per-stream reencode flag to the global config flag.
     * stream.c gates the entire reencode path on config->reencode; without
     * this loop, <Reencode>1</Reencode> in ices.conf does nothing because
     * the XML parser only sets stream->reencode, not config->reencode. */
    for (stream = ices_config.streams; stream; stream = stream->next) {
        if (stream->reencode) {
            ices_config.reencode = 1;
            break;
        }
    }
    if (!ices_config.reencode) {
        LOG("no <Stream> has <Reencode>1</Reencode>; reencoder disabled");
        return;
    }

    /* Don't init LAME here — we don't know the output format yet (it
     * arrives via the ices_stream_t on the first ices_reencode call).
     * LAME is brought up lazily then. Just bring up decoder + leveller. */
    if (hip_decoder) hip_decode_exit(hip_decoder);
    hip_decoder = hip_decode_init();
    if (!hip_decoder) LOG("hip_decode_init() failed");

    memset(&state_l, 0, sizeof(state_l));
    memset(&state_r, 0, sizeof(state_r));
    for (int b = 0; b < NBANDS; b++) {
        dynamic_gains[b] = 1.0;
        band_envelope[b] = 0.0;
    }
    recompute_coeffs(44100);

    LOG("5-band leveller initialised (LAME init deferred until first track)");
}

void ices_reencode_reset(input_stream_t* source) {
    /* Reset leveller state between tracks. */
    memset(&state_l, 0, sizeof(state_l));
    memset(&state_r, 0, sizeof(state_r));
    for (int b = 0; b < NBANDS; b++) {
        dynamic_gains[b] = 1.0;
        band_envelope[b] = 0.0;
    }

    /* Reset the MP3 decoder so partial-frame state from the previous file
     * doesn't poison the next one. */
    if (hip_decoder) hip_decode_exit(hip_decoder);
    hip_decoder = hip_decode_init();
    if (!hip_decoder) LOG("hip_decode_init() failed during reset");

    /* Remember the source. We can't reconfigure LAME yet (we don't know
     * the output format until ices_reencode is called with the stream
     * pointer), so defer the reconfigure to first encode call. */
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
    if (n < 0) {
        /* Real decode error. Returning -1 lets stream.c abort the track
         * cleanly. Returning 0 here would spin the loop forever. */
        return -1;
    }
    return n;
}

int ices_reencode(ices_stream_t* stream, int nsamples,
                  int16_t* left, int16_t* right,
                  unsigned char* outbuf, int outbuf_sz)
{
    if (nsamples <= 0) return 0;

    /* Pick up output config from ices.conf (carried on ices_stream_t).
     * This is the only place that pointer reaches us. */
    int out_changed = absorb_stream_config(stream);

    /* If a track change is pending, reconfigure LAME now that we know
     * both the input format (from pending_source) and the output format
     * (from stream). */
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
            return -2;   /* < -1: fatal, stream.c will abort the track */
        }
    }

    if (!gfp) return -2;

    /* 5-band leveller in-place. */
    process_audio_channels(left, right, nsamples);

    int n = lame_encode_buffer(gfp, left, right, nsamples, outbuf, outbuf_sz);

    if (n == -1) {
        /* Output buffer too small. Returning -1 tells stream.c to grow
         * the buffer and retry — it specifically checks for this case. */
        return -1;
    }
    if (n < -1) {
        LOG("lame_encode_buffer returned %d (samples=%d outbuf_sz=%d)",
            n, nsamples, outbuf_sz);
        return n;   /* propagate fatal */
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
