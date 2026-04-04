/*  vecstream_adaptive.h — Adaptive delta tolerance and keyframe scheduling
 *
 *  Self-tuning codec parameters based on observed delta statistics.
 *  Include after vecstream.h.
 *
 *  Usage:
 *    VecAdaptive va;
 *    vec_adaptive_init(&va, initial_tol, target_delta_fraction);
 *
 *    // After each P-frame:
 *    vec_adaptive_update(&va, n_patches, n_deltas, max_delta);
 *    float new_tol = vec_adaptive_tolerance(&va);
 *    int should_iframe = vec_adaptive_needs_iframe(&va);
 *    int should_kframe = vec_adaptive_needs_kframe(&va);
 */

#ifndef VECSTREAM_ADAPTIVE_H
#define VECSTREAM_ADAPTIVE_H

#include <math.h>

typedef struct {
    /* Configuration */
    float target_delta_fraction;  /* target: what fraction of patches have deltas (0.3-0.7) */
    float min_tol;                /* floor — never go below this */
    float max_tol;                /* ceiling — never go above this */
    float error_budget;           /* max accumulated drift before forcing I-frame */

    /* Current state */
    float current_tol;
    float smoothed_fraction;      /* EMA of actual delta fraction */
    float accumulated_drift;      /* estimated error from chained P-frames */
    int frames_since_iframe;
    int frames_since_kframe;

    /* Statistics */
    float max_delta_seen;         /* largest single-frame max_delta */
    float avg_delta_fraction;     /* running average */
    int total_iframes;
    int total_pframes;
    int total_kframes;
    int auto_iframes;             /* I-frames triggered by drift */
    int auto_kframes;             /* K-frames triggered by error */
} VecAdaptive;

static void vec_adaptive_init(VecAdaptive *va, float initial_tol, float target_fraction) {
    memset(va, 0, sizeof(VecAdaptive));
    va->current_tol = initial_tol;
    va->target_delta_fraction = (target_fraction > 0) ? target_fraction : 0.5f;
    va->min_tol = initial_tol * 0.01f;
    va->max_tol = initial_tol * 100.0f;
    va->error_budget = initial_tol * 50.0f;  /* allow ~50 P-frames of drift before I-frame */
    va->smoothed_fraction = target_fraction;
}

/* Call after each P-frame with the delta statistics */
static void vec_adaptive_update(VecAdaptive *va, int n_patches, int n_deltas, float max_delta) {
    float fraction = (n_patches > 0) ? (float)n_deltas / n_patches : 0;

    /* EMA smoothing (α = 0.1) */
    va->smoothed_fraction = 0.9f * va->smoothed_fraction + 0.1f * fraction;

    /* Track max delta */
    if (max_delta > va->max_delta_seen)
        va->max_delta_seen = max_delta;

    /* Accumulate drift estimate: each P-frame adds ~tol of potential error */
    va->accumulated_drift += va->current_tol;
    va->frames_since_iframe++;
    va->frames_since_kframe++;
    va->total_pframes++;

    /* Adjust tolerance toward target fraction */
    if (va->smoothed_fraction > va->target_delta_fraction * 1.2f) {
        /* Too many deltas — relax tolerance */
        va->current_tol *= 1.05f;
    } else if (va->smoothed_fraction < va->target_delta_fraction * 0.8f) {
        /* Too few deltas — tighten tolerance for better quality */
        va->current_tol *= 0.95f;
    }

    /* Clamp */
    if (va->current_tol < va->min_tol) va->current_tol = va->min_tol;
    if (va->current_tol > va->max_tol) va->current_tol = va->max_tol;
}

/* Get the current adaptive tolerance */
static float vec_adaptive_tolerance(const VecAdaptive *va) {
    return va->current_tol;
}

/* Should we force an I-frame? (drift exceeded budget) */
static int vec_adaptive_needs_iframe(VecAdaptive *va) {
    if (va->accumulated_drift >= va->error_budget) {
        va->accumulated_drift = 0;
        va->frames_since_iframe = 0;
        va->auto_iframes++;
        va->total_iframes++;
        return 1;
    }
    return 0;
}

/* Should we force a K-frame for verification? */
static int vec_adaptive_needs_kframe(VecAdaptive *va) {
    /* K-frame every ~10 I-frames, or if drift is high */
    if (va->frames_since_kframe > 0 &&
        va->total_iframes > 0 &&
        va->total_iframes % 10 == 0) {
        va->frames_since_kframe = 0;
        va->auto_kframes++;
        va->total_kframes++;
        return 1;
    }
    return 0;
}

/* Record that an I-frame was written (resets drift) */
static void vec_adaptive_iframe_written(VecAdaptive *va) {
    va->accumulated_drift = 0;
    va->frames_since_iframe = 0;
    va->total_iframes++;
}

/* Record that a K-frame was written */
static void vec_adaptive_kframe_written(VecAdaptive *va) {
    va->frames_since_kframe = 0;
    va->total_kframes++;
}

/* Print summary */
static void vec_adaptive_summary(const VecAdaptive *va) {
    fprintf(stderr, "VecAdaptive: tol=%.6f (range %.6f-%.6f)\n",
            va->current_tol, va->min_tol, va->max_tol);
    fprintf(stderr, "  delta fraction: %.1f%% (target %.1f%%)\n",
            va->smoothed_fraction * 100, va->target_delta_fraction * 100);
    fprintf(stderr, "  frames: %d I (%d auto), %d P, %d K (%d auto)\n",
            va->total_iframes, va->auto_iframes,
            va->total_pframes,
            va->total_kframes, va->auto_kframes);
    fprintf(stderr, "  max delta seen: %.6f\n", va->max_delta_seen);
}

#endif /* VECSTREAM_ADAPTIVE_H */
