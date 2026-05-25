/*
 * ices0 - Icecast Source Client
 * Copyright (C) 2000-2004 The Icecast Team <team@icecast.org>
 *
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
 * Each band has its own slow AGC-style leveller so quiet tracks come up
 * and loud tracks come down, *without* pumping. Per-band coefficients are
 * recomputed when the sample rate changes (per track).
 * ========================================================================= */

#define NBANDS 5

typedef struct {
    /* two-pole IIR low-pass states for each crossover */
    double p1_0, p1_1;
    double p2_0, p2_1;
    double p3_0, p3_1;
    double p4_0, p4_1;
} FilterState;

static FilterState state_l, state_r;

/* One-pole coefficients (normalised) for each crossover. Recomputed in reset. */
static double xover[4] = {
    100.0  / 44100.0,
    400.0  / 44100.0,
    1500.0 / 44100.0,
    5000.0 / 44100.0
};

/* Target RMS-ish energy per band. These are tuned to a typical loudness curve
 * (a little less air than mids) and are mostly arbitrary — the leveller is
 * slow enough that the overall loudness, not the band balance, dominates. */
static const double target_energy[NBANDS] = { 0.18, 0.20, 0.22, 0.18, 0.12 };

/* Adaptive per-band gain, smoothed envelope, and per-band attack speed. */
static double dynamic_gains[NBANDS]    = { 1.0, 1.0, 1.0, 1.0, 1.0 };
static double band_envelope[NBANDS]    = { 0.0, 0.0, 0.0, 0.0, 0.0 };

/* Smoothing coefficients (recomputed per sample rate in reset).
 * envelope is fast (~30 ms), gain adapts slowly (~2 s) — broadcast-style. */
static double env_coef   = 1.0 / (0.030 * 44100.0);   /* ~30 ms */
static double gain_coef  = 1.0 / (2.0   * 44100.0);   /* ~2 s   */

/* Simple soft-knee limiter on the sum (avoids hard clipping when bands stack).
 * y = x - x^3/3, clamped to [-1, +1]. Transparent below ~0.6, gentle above. */
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

    /* clamp so a low sample rate doesn't push a crossover above Nyquist */
    for (int i = 0; i < 4; i++) {
        if (xover[i] > 0.45) xover[i] = 0.45;
    }

    env_coef  = 1.0 / (0.030 * sr);
    gain_coef = 1.0 / (2.0   * sr);
}

static void split_5_bands(double input, FilterState *s, double *bands) {
    /* Cascaded one-pole low-passes (two stages each = 12 dB/oct).
     * Bands are differences between successive low-passes. */

    s->p1_0 += xover[0] * (input    - s->p1_0);
    s->p1_1 += xover[0] * (s->p1_0  - s->p1_1);
    double lp1 = s->p1_1;

    s->p2_0 += xover[1] * (input    - s->p2_0);
    s->p2_1 += xover[1] * (s->p2_0  - s->p2_1);
    double lp2 = s->p2_1;

    s->p3_0 += xover[2] * (input    - s->p3_0);
    s->p3_1 += xover[2] * (s->p3_0  - s->p3_1);
    double lp3 = s->p3_1;

    s->p4_0 += xover[3] * (input    - s->p4_0);
    s->p4_1 += xover[3] * (s->p4_0  - s->p4_1);
    double lp4 = s->p4_1;

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
            /* Mono band energy estimate, smoothed by env_coef. */
            double inst = 0.5 * (fabs(bands_l[b]) + fabs(bands_r[b]));
            band_envelope[b] += (inst - band_envelope[b]) * env_coef;

            /* Only update gain when the band is actually active — this
             * prevents the gain running away on silence between tracks. */
            if (band_envelope[b] > 0.001) {
                double desired = target_energy[b] / band_envelope[b];
                dynamic_gains[b] += (desired - dynamic_gains[b]) * gain_coef;
            }

            /* Sane per-band gain range. */
            if (dynamic_gains[b] > 4.0) dynamic_gains[b] = 4.0;
            if (dynamic_gains[b] < 0.25) dynamic_gains[b] = 0.25;

            mix_l += bands_l[b] * dynamic_gains[b];
            mix_r += bands_r[b] * dynamic_gains[b];
        }

        /* Soft-knee limiter on the recombined signal. */
        mix_l = softclip(mix_l);
        mix_r = softclip(mix_r);

        /* Hard safety clip (should never actually trigger after softclip). */
        if (mix_l >  0.999) mix_l =  0.999;
        if (mix_l < -0.999) mix_l = -0.999;
        if (mix_r >  0.999) mix_r =  0.999;
        if (mix_r < -0.999) mix_r = -0.999;

        left[i]  = (int16_t)(mix_l * 32767.0);
        right[i] = (int16_t)(mix_r * 32767.0);
    }
}

/* =========================================================================
 * ICES0 REENCODER INTERFACE
 * ========================================================================= */

static lame_global_flags *gfp         = NULL;
static hip_t              hip_decoder = NULL;

/* Cached config so we know what LAME is currently configured for; if the
 * next track differs, we tear down and re-init LAME. */
static int cur_in_rate    = 0;
static int cur_in_chans   = 0;
static int cur_out_rate   = 0;
static int cur_out_chans  = 0;
static int cur_out_brate  = 0;

static int lame_configure(int in_rate, int in_chans,
                          int out_rate, int out_chans, int out_brate)
{
    if (gfp) {
        lame_close(gfp);
        gfp = NULL;
    }

    gfp = lame_init();
    if (!gfp) {
        fprintf(stderr, "reencode: lame_init() failed\n");
        return -1;
    }

    lame_set_in_samplerate (gfp, in_rate);
    lame_set_num_channels  (gfp, in_chans);
    lame_set_out_samplerate(gfp, out_rate > 0 ? out_rate : in_rate);
    lame_set_mode          (gfp, out_chans == 1 ? MONO : JOINT_STEREO);
    lame_set_brate         (gfp, out_brate > 0 ? out_brate : 128);
    lame_set_quality       (gfp, 2);
    /* Don't write Xing/Info header to the middle of a live stream. */
    lame_set_bWriteVbrTag  (gfp, 0);

    if (lame_init_params(gfp) < 0) {
        fprintf(stderr, "reencode: lame_init_params() failed\n");
        lame_close(gfp);
        gfp = NULL;
        return -1;
    }

    cur_in_rate   = in_rate;
    cur_in_chans  = in_chans;
    cur_out_rate  = out_rate > 0 ? out_rate : in_rate;
    cur_out_chans = out_chans;
    cur_out_brate = out_brate > 0 ? out_brate : 128;

    return 0;
}

void ices_reencode_initialize(void) {
    /* Default startup config; will be re-done per track in reset(). */
    if (lame_configure(44100, 2, 44100, 2, 128) < 0) {
        fprintf(stderr, "reencode: initial LAME configuration failed\n");
    }

    hip_decoder = hip_decode_init();
    if (!hip_decoder) {
        fprintf(stderr, "reencode: hip_decode_init() failed\n");
    }

    /* Filter / leveller state. */
    memset(&state_l, 0, sizeof(state_l));
    memset(&state_r, 0, sizeof(state_r));
    for (int b = 0; b < NBANDS; b++) {
        dynamic_gains[b] = 1.0;
        band_envelope[b] = 0.0;
    }
    recompute_coeffs(44100);

    fprintf(stderr, "reencode: 5-band leveller initialised\n");
}

void ices_reencode_reset(input_stream_t* source) {
    /* Reset filter and adaptive-gain state between tracks so a quiet outro
     * doesn't leave the next track 4x too loud. */
    memset(&state_l, 0, sizeof(state_l));
    memset(&state_r, 0, sizeof(state_r));
    for (int b = 0; b < NBANDS; b++) {
        dynamic_gains[b] = 1.0;
        band_envelope[b] = 0.0;
    }

    /* Reset the MP3 decoder so it doesn't carry partial-frame state from the
     * previous file into the next one. */
    if (hip_decoder) {
        hip_decode_exit(hip_decoder);
    }
    hip_decoder = hip_decode_init();

    /* Reconfigure LAME from the source's actual rate / channels — without
     * this, anything that isn't 44.1 kHz stereo causes lame_encode_buffer
     * to return errors on every call. */
    if (source) {
        int in_rate  = (int)source->samplerate;
        int in_chans = (int)source->channels;
        if (in_rate  <= 0) in_rate  = 44100;
        if (in_chans <= 0) in_chans = 2;

        /* If ices.conf overrides output, those are already in cur_out_*.
         * Default to keeping the input format. */
        int out_rate  = cur_out_rate  > 0 ? cur_out_rate  : in_rate;
        int out_chans = cur_out_chans > 0 ? cur_out_chans : in_chans;
        int out_brate = cur_out_brate > 0 ? cur_out_brate : 128;

        if (in_rate != cur_in_rate || in_chans != cur_in_chans) {
            if (lame_configure(in_rate, in_chans,
                               out_rate, out_chans, out_brate) < 0) {
                fprintf(stderr,
                        "reencode: failed to reconfigure LAME for %d Hz / %d ch\n",
                        in_rate, in_chans);
            }
        }

        recompute_coeffs((unsigned int)in_rate);
    }
}

/* Decode MP3 bytes into PCM int16 left/right.
 *
 * Return value convention (matches the rest of ices0):
 *   > 0  : number of samples per channel decoded into left[]/right[]
 *   == 0 : no samples yet (decoder needs more input bytes — normal)
 *   < 0  : fatal decode error (caller will abort the track)
 */
int ices_reencode_decode(unsigned char* buf, size_t blen, size_t olen,
                         int16_t* left, int16_t* right) {
    if (!hip_decoder) {
        return -1;
    }

    /* hip_decode wants the per-channel output buffer size in samples. */
    int max_samples = (int)(olen / sizeof(int16_t));
    (void)max_samples;  /* hip_decode doesn't take a size; left/right must be sized
                         * generously by the caller, which stream.c does
                         * (INPUT_BUFSIZ * 45 samples each). */

    int n = hip_decode(hip_decoder, buf, blen, left, right);

    if (n < 0) {
        /* Real decode error. Returning 0 here would make the caller spin
         * forever; returning -1 lets stream.c flag the track as failed
         * and move on cleanly. */
        return -1;
    }
    return n;
}

int ices_reencode(ices_stream_t* stream, int nsamples,
                  int16_t* left, int16_t* right,
                  unsigned char* outbuf, int outbuf_sz)
{
    (void)stream;
    if (!gfp || nsamples <= 0) {
        return 0;
    }

    /* Apply 5-band leveller in-place. */
    process_audio_channels(left, right, nsamples);

    int n = lame_encode_buffer(gfp, left, right, nsamples, outbuf, outbuf_sz);

    if (n == -1) {
        /* Output buffer too small — stream.c handles this by growing the
         * buffer and retrying. Must return -1 exactly, not 0. */
        return -1;
    }
    if (n < -1) {
        /* Any other negative is a real LAME error; propagate so stream.c
         * aborts the track instead of looping. */
        fprintf(stderr, "reencode: lame_encode_buffer error %d\n", n);
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
    if (hip_decoder) {
        hip_decode_exit(hip_decoder);
        hip_decoder = NULL;
    }
    if (gfp) {
        lame_close(gfp);
        gfp = NULL;
    }
}
