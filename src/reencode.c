/* ices0 - Icecast Source Client
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

/*#include "ices.h" */
#include "reencode.h"
#include "log.h"

/* =========================================================================
 * PROFESSIONAL 5-BAND MULTIBAND DYNAMIC BROADCAST PROCESSOR
 * ========================================================================= */

// Crossover cutoffs normalized to 44.1kHz sample rate
// Sub/Low Bass (~100Hz), Mid-Bass (~400Hz), Vocal Mids (~1.5kHz), Presence (~5kHz)
#define FREQ_BAND_1 (100.0 / 44100.0)
#define FREQ_BAND_2 (400.0 / 44100.0)
#define FREQ_BAND_3 (1500.0 / 44100.0)
#define FREQ_BAND_4 (5000.0 / 44100.0)

typedef struct {
    // IIR Low-pass state histories for cascading splits
    double p1_0, p1_1;
    double p2_0, p2_1;
    double p3_0, p3_1;
    double p4_0, p4_1;
} FilterState;

static FilterState state_l = {0}, state_r = {0};

// Target root-mean-square (RMS) energy densities for a loud, commercial master profile
static const double target_energy[5] = {
    0.16,  // Band 1: Low Bass (Body, Punch)
    0.18,  // Band 2: Mid Bass (Warmth, Thick)
    0.15,  // Band 3: Mid Range (Vocals, Guitars)
    0.12,  // Band 4: Mid Highs (Clarity, Snare snap)
    0.09   // Band 5: Highs (Air, Crispness, Cymbals)
};

// Dynamic state tracking for independent AGC faders
static double dynamic_gains[5] = {1.0, 1.0, 1.0, 1.0, 1.0};

/* 5-Band IIR Filter Splitting Engine */
static void split_5_bands(double input, FilterState *state, double *bands) {
    // Recursive 2nd order low pass filtering to isolate bands cleanly
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

    // Isolate chunks via phase subtraction methods
    bands[0] = lp1;                 // Low Bass
    bands[1] = lp2 - lp1;           // Mid Bass
    bands[2] = lp3 - lp2;           // Mid Range
    bands[3] = lp4 - lp3;           // Mid Highs
    bands[4] = input - lp4;         // Highs
}

/* Core Audio DSP Processing Routine */
static void process_audio_frame(short *buffer, int samples_read) {
    double bands_l[5], bands_r[5];

    for (int i = 0; i < samples_read; i += 2) {
        // Convert to normalized floating point domains
        double sample_l = (double)buffer[i]   / 32768.0;
        double sample_r = (double)buffer[i+1] / 32768.0;

        // Split audio spectrum into 5 independent structural vectors
        split_5_bands(sample_l, &state_l, bands_l);
        split_5_bands(sample_r, &state_r, bands_r);

        // Calculate momentary aggregate energy inside the 5 bands
        for (int b = 0; b < 5; b++) {
            double absolute_energy = fabs(bands_l[b]) + fabs(bands_r[b]);
            
            if (absolute_energy > 0.0005) {
                double computed_gain = target_energy[b] / absolute_energy;
                // High-inertia leveling coefficients (avoids sound pumping artifacting)
                dynamic_gains[b] += (computed_gain - dynamic_gains[b]) * 0.00008;
            }
            
            // Apply safety bounds on hardware AGC expansion to protect silence
            if (dynamic_gains[b] > 4.0)  dynamic_gains[b] = 4.0;
            if (dynamic_gains[b] < 0.1)  dynamic_gains[b] = 0.1;
            
            // Re-level the decoupled structural audio segments
            bands_l[b] *= dynamic_gains[b];
            bands_r[b] *= dynamic_gains[b];
        }

        // Remux the vectors into a uniform stereo matrix
        double final_l = bands_l[0] + bands_l[1] + bands_l[2] + bands_l[3] + bands_l[4];
        double final_r = bands_r[0] + bands_r[1] + bands_r[2] + bands_r[3] + bands_r[4];

        // Broadcast Hard-Clipping Peak Limiter to eliminate digital overs
        if (final_l > 0.98) final_l = 0.98; else if (final_l < -0.98) final_l = -0.98;
        if (final_r > 0.98) final_r = 0.98; else if (final_r < -0.98) final_r = -0.98;

        // Flatten back into Interleaved standard signed 16-Bit format
        buffer[i]   = (short)(final_l * 32767.0);
        buffer[i+1] = (short)(final_r * 32767.0);
    }
}

/* =========================================================================
 * ICES0 SOURCE CORE ENGINE INTERFACING
 * ========================================================================= */

static lame_global_flags *gfp = NULL;

int reencode_init(void) {
    gfp = lame_init();
    if (gfp == NULL) {
        LOG_ERROR0("Failed to initialize LAME encoder engine library context.");
        return -1;
    }

    // Inherit configurations mapped out via icecast XML configuration formats
    lame_set_in_samplerate(gfp, 44100);
    lame_set_num_channels(gfp, 2);
    lame_set_out_samplerate(gfp, 44100);
    
    // Configure target bitrate defaults matching ices requirements
    lame_set_brate(gfp, 128); 
    lame_set_quality(gfp, 2); // High quality audio compilation configuration
    
    if (lame_init_params(gfp) < 0) {
        LOG_ERROR0("LAME dynamic parameter parsing setup failed.");
        return -1;
    }

    LOG_INFO0("Hardware-Style 5-Band Dynamics Master Engine compiled & initialized successfully.");
    return 0;
}

int reencode_data(short *pcm_buf, int samples, unsigned char *mp3_buf, int mp3_buf_sz) {
    if (gfp == NULL) {
        return -1;
    }

    // Intercept decoded PCM buffers and execute the 5-band leveling sequence
    process_audio_frame(pcm_buf, samples * 2);

    // Encode the perfectly processed master straight to the Icecast target pipeline
    int bytes_encoded = lame_encode_buffer_interleaved(gfp, pcm_buf, samples, mp3_buf, mp3_buf_sz);
    if (bytes_encoded < 0) {
        LOG_ERROR1("LAME system encoding execution anomaly detected: %d", bytes_encoded);
    }

    return bytes_encoded;
}

int reencode_flush(unsigned char *mp3_buf, int mp3_buf_sz) {
    if (gfp == NULL) {
        return -1;
    }
    return lame_encode_flush(gfp, mp3_buf, mp3_buf_sz);
}

void reencode_close(void) {
    if (gfp != NULL) {
        lame_close(gfp);
        gfp = NULL;
    }
}
