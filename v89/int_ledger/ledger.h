/* Integer ledger helpers for cellfabi.c — see v89/INT_LEDGER.md
 *
 * ledger_mode:
 *   0 — off (FP64 only; should match cellfab when ledger_mode absent)
 *   1 — shadow: FP64 dynamics + parallel int64 accounts (quantized Δ)
 *   2 — dense+space truth: Es/Em integer; FP synced for kinematics
 *   3 — full: + lem + conversion quanta + field energy inject bridge
 */
#ifndef CELLFABI_LEDGER_H
#define CELLFABI_LEDGER_H

#include <stdint.h>
#include <math.h>

/* residual + quantize: deterministic floor toward -inf of (x + resid) */
static inline int64_t led_qdelta(double x_over_u, double *resid)
{
    double raw = x_over_u + *resid;
    double fl = floor(raw);
    /* guard huge values */
    if (fl > (double)INT64_MAX) fl = (double)INT64_MAX;
    if (fl < (double)INT64_MIN) fl = (double)INT64_MIN;
    int64_t dn = (int64_t)fl;
    *resid = raw - fl;
    return dn;
}

static inline int64_t led_from_fp(double e, double u)
{
    if (u <= 0) return 0;
    double q = floor(e / u + 0.5);
    if (q > (double)INT64_MAX) return INT64_MAX;
    if (q < (double)INT64_MIN) return INT64_MIN;
    return (int64_t)q;
}

static inline double led_to_fp(int64_t n, double u)
{
    return (double)n * u;
}

#endif
