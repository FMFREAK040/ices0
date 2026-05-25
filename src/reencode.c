/*
 * ices0 - Icecast Source Client
 * Copyright (C) 2000-2004 The Icecast Team <team@icecast.org>
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

// Include libshout headers first so icestypes.h knows what shout_t is!
#include <shout/shout.h>

// Match the exact header layout of your repository branch
#include "icestypes.h"
#include "ices_config.h"
#include "reencode.h"

/* =========================================================================
 * PROFESSIONAL 5-BAND MULTIBAND DYNAMIC BROADCAST PROCESSOR
 * ========================================================================= */

#define FREQ_BAND_1 (100.0 / 44100.0)
#define FREQ_BAND_2 (400.0 / 44100.0)
#define FREQ_BAND_3 (1500.0 / 44100.0)
#define FREQ_BAND_4 (5000.0 / 44100.0)

typedef struct {
    double p1_0, p1_1;
    double p2_0, p2_1;
    double p3_0, p3_1;
    double p4_0, p4_1;
} FilterState;

static FilterState state_l = {0}, state_r = {0};

static const double target_energy[5] = {0.16, 0.18, 0.15, 0.12, 0.09};
static double dynamic_gains[5] = {1.0, 1.0, 1.0, 1.0, 1.0};

static void split_5_bands(double input, FilterState *state, double *bands) {
    state->p1_0 += FREQ_BAND_1 * (input - state->p1_0);
    state->p1_1 += FREQ_BAND_1 * (state->p1_0 - state->p1_1);
    double lp1 = state->p1_1;

    state->p2_0 += FREQ_BAND_2 * (input - state->p2_0);
    state->p2_1 += FREQ_BAND_2 * (state->p2_0 - state->p2_1);
    double lp2 = state->p2_1;

    state->p3_0 += FREQ_BAND_3 * (input - state->p3_0);
    state->p3_1 += FREQ_BAND_3 * (state->p3_0 - state->p3_1);
    double lp3 = state->p3_1;

    state->p4_0 += FREQ_BAND_4 * (input - state->p4_0);
    state->p4_1 += FREQ_BAND_4 * (state->p4_0 - state->p4_1);
    double lp4 = state->p4_1;

    bands[0] = lp1;                 
    bands[1] = lp2 - lp1;           
    bands[2] = lp3 - lp2;           
    bands[3] = lp4 - lp3;           
    bands[4] = input - lp4;         
}

static void process_audio_channels(int16_t *left, int16_t *right, int nsamples) {
    double bands_l[5], bands_r[5];

    for (int i = 0; i < nsamples; i++) {
        double sample_l = (double)left[i]  / 32768.0;
        double sample_r = (double)right[i] / 32768.0;

        split_5_bands(sample_l, &state_l, bands_l);
        split_5_bands(sample_r, &state_r, bands_r);

        for (int b = 0; b < 5; b++) {
            double absolute_energy = fabs(bands_l[b]) + fabs(bands_r[b]);
            
            if (absolute_energy > 0.0005) {
                double computed_gain = target_energy[b] / absolute_energy;
                dynamic_gains[b] += (computed_gain - dynamic_gains[b]) * 0.00008;
            }
            
            if (dynamic_gains[b] > 4.0)  dynamic_gains[b] = 4.0;
            if (dynamic_gains[b] < 0.1)  dynamic_gains[b] = 0.1;
            
            bands_l[b] *= dynamic_gains[b];
            bands_r[b] *= dynamic_gains[b];
        }

        double final_l = bands_l[0] + bands_l[1] + bands_l[2] + bands_l[3] + bands_l[4];
        double final_r = bands_r[0] + bands_r[1] + bands_r[2] + bands_r[3] + bands_r[4];

        if (final_l > 0.98) final_l = 0.98; else if (final_l < -0.98) final_l = -0.98;
        if (final_r > 0.98) final_r = 0.98; else if (final_r < -0.98) final_r = -0.98;

        left[i]  = (int16_t)(final_l * 32767.0);
        right[i] = (int16_t)(final_r * 32767.0);
    }
}

/* =========================================================================
 * ICES0 RECONCILED INTERFACE IMPLEMENTATION
 * ========================================================================= */

static lame_global_flags *gfp = NULL;

void ices_reencode_initialize(void) {
    gfp = lame_init();
    if (gfp == NULL) {
        fprintf(stderr, "ERROR: Failed to initialize LAME encoder engine library context.\n");
        return;
    }

    lame_set_in_samplerate(gfp, 44100);
    lame_set_num_channels(gfp, 2);
    lame_set_out_samplerate(gfp, 44100);
    
    lame_set_brate(gfp, 128); 
    lame_set_quality(gfp, 2); 
    
    if (lame_init_params(gfp) < 0) {
        fprintf(stderr, "ERROR: LAME dynamic parameter parsing setup failed.\n");
        return;
    }

    printf("INFO: Hardware-Style 5-Band Dynamics Master Engine compiled & initialized successfully.\n");
}

void ices_reencode_reset(input_stream_t* source) {
    (void)source;
    memset(&state_l, 0, sizeof(FilterState));
    memset(&state_r, 0, sizeof(FilterState));
}

int ices_reencode_decode(unsigned char* buf, size_t blen, size_t olen, int16_t* left, int16_t* right) {
    (void)buf; (void)blen; (void)olen; (void)left; (void)right;
    return 0; 
}

int ices_reencode(ices_stream_t* stream, int nsamples, int16_t* left, int16_t* right, unsigned char* outbuf, int outbuf_sz) {
    if (gfp == NULL) {
        return -1;
    }
    (void)stream; 

    // Run our AGC / 5-band Equalizer
    process_audio_channels(left, right, nsamples);

    // CRITICAL FIX: Pass 'nsamples' directly without doubling it!
    int bytes_encoded = lame_encode_buffer(gfp, left, right, nsamples, outbuf, outbuf_sz);
    
    if (bytes_encoded < 0) {
        fprintf(stderr, "WARNING: LAME non-interleaved encoding execution anomaly detected: %d\n", bytes_encoded);
        // Fallback safety to prevent stream stalling or hard loops if a frame fails
        return 0; 
    }

    return bytes_encoded;
}

int ices_reencode_flush(ices_stream_t* stream, unsigned char *outbuf, int outbuf_sz) {
    if (gfp == NULL) {
        return -1;
    }
    (void)stream;
    return lame_encode_flush(gfp, outbuf, outbuf_sz);
}

void ices_reencode_shutdown(void) {
    if (gfp != NULL) {
        lame_close(gfp);
        gfp = NULL;
    }
}
