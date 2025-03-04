/*
===========================================================================
Copyright (C) 1997-2006 Id Software, Inc.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
===========================================================================
*/

#ifndef __MATHLIB__
#define __MATHLIB__

#define MATH_INLINE

// mathlib.h

#include <math.h>

// big enough without being so big we get major floating point errors.
#define BOGUS_RANGE (1 << 20) //qb: to comply with 24-bit 20.3 network max.  1048576

#ifdef BSP          // only do this for qbsp, leads to stack overflows on qrad3.
#define DOUBLEVEC_T // jit - might as well be more accurate, and sometimes doubles are even faster on modern hardware, anyway...
#endif

#ifdef DOUBLEVEC_T
typedef double vec_t;
#else
typedef float vec_t;
#endif
typedef vec_t vec3_t[3];

#define SIDE_FRONT 0
#define SIDE_ON    2
#define SIDE_BACK  1
#define SIDE_CROSS -2

#define Q_PI       3.14159265358979323846

// angle indexes
#define PITCH      0 // up / down
#define YAW        1 // left / right
#define ROLL       2 // fall over

extern vec3_t vec3_origin;

#define EQUAL_EPSILON               0.001
#define NORMAL_EPSILON              0.00001

#define DotProduct(x, y)            (x[0] * y[0] + x[1] * y[1] + x[2] * y[2])
#define CrossProduct(v1, v2, cross) ((cross)[0] = (v1)[1] * (v2)[2] - (v1)[2] * (v2)[1], (cross)[1] = (v1)[2] * (v2)[0] - (v1)[0] * (v2)[2], (cross)[2] = (v1)[0] * (v2)[1] - (v1)[1] * (v2)[0])

#define VectorSubtract(a, b, c)     (c[0] = a[0] - b[0], c[1] = a[1] - b[1], c[2] = a[2] - b[2])
#define VectorAdd(a, b, c)          (c[0] = a[0] + b[0], c[1] = a[1] + b[1], c[2] = a[2] + b[2])
#define VectorCopy(a, b)            (b[0] = a[0], b[1] = a[1], b[2] = a[2])
#define VectorScale(a, b, c)        (c[0] = b * (a)[0], c[1] = b * (a)[1], c[2] = b * (a)[2])

#define VectorClear(x)              (x[0] = x[1] = x[2] = 0)
#define VectorNegate(x)             (x[0] = -x[0], x[1] = -x[1], x[2] = -x[2])
#define VectorAverage(a, b, c)      (c[0] = (a[0] + b[0]) / 2, c[1] = (a[1] + b[1]) / 2, c[2] = (a[2] + b[2]) / 2)

#define VectorLength(v)             (sqrt(DotProduct((v), (v))))
#define VectorInverse(v)            ((v)[0] = -(v)[0], (v)[1] = -(v)[1], (v)[2] = -(v)[2])
#define VectorMA(a, b, c, d)        ((d)[0] = (a)[0] + (b) * (c)[0], (d)[1] = (a)[1] + (b) * (c)[1], (d)[2] = (a)[2] + (b) * (c)[2])

#define Q_rint(in)                  (floor((in) + 0.5))

#define MIN(X, Y)                   (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y)                   (((X) > (Y)) ? (X) : (Y))
#define BOUND(a, b, c)              ((a) >= (c) ? (a) : (b) < (a) ? (a) \
                                       : (b) > (c)                ? (c) \
                                                                  : (b))

#ifndef MATH_INLINE
qboolean VectorCompare(vec3_t v1, vec3_t v2);

vec_t _DotProduct(vec3_t v1, vec3_t v2);
void _VectorSubtract(vec3_t va, vec3_t vb, vec3_t out);
void _VectorAdd(vec3_t va, vec3_t vb, vec3_t out);
void _VectorCopy(vec3_t in, vec3_t out);
void _VectorScale(vec3_t v, vec_t scale, vec3_t out);

vec_t VectorNormalize(vec3_t in, vec3_t out);

vec_t ColorNormalize(vec3_t in, vec3_t out);

void ClearBounds(vec3_t mins, vec3_t maxs);
void AddPointToBounds(vec3_t v, vec3_t mins, vec3_t maxs);

qboolean RayPlaneIntersect(vec3_t p_n, vec_t p_d, vec3_t l_o, vec3_t l_n,
                           vec3_t res);

#else

static inline double VectorLengthSq(vec3_t v) // jit
{
    return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
}

static inline qboolean VectorCompare(vec3_t v1, vec3_t v2) {
    int32_t i;

    for (i = 0; i < 3; i++)
        if (fabs(v1[i] - v2[i]) > EQUAL_EPSILON)
            return false;

    return true;
}

static inline vec_t _DotProduct(vec3_t v1, vec3_t v2) {
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

static inline void _VectorSubtract(vec3_t va, vec3_t vb, vec3_t out) {
    out[0] = va[0] - vb[0];
    out[1] = va[1] - vb[1];
    out[2] = va[2] - vb[2];
}

static inline void _VectorAdd(vec3_t va, vec3_t vb, vec3_t out) {
    out[0] = va[0] + vb[0];
    out[1] = va[1] + vb[1];
    out[2] = va[2] + vb[2];
}

static inline void _VectorCopy(vec3_t in, vec3_t out) {
    out[0] = in[0];
    out[1] = in[1];
    out[2] = in[2];
}

static inline void _VectorScale(vec3_t v, vec_t scale, vec3_t out) {
    out[0] = v[0] * scale;
    out[1] = v[1] * scale;
    out[2] = v[2] * scale;
}

static inline vec_t VectorNormalize(vec3_t in, vec3_t out) {
    vec_t length, ilength;

    length = sqrt(in[0] * in[0] + in[1] * in[1] + in[2] * in[2]);
    if (length == 0) {
        VectorClear(out);
        return 0;
    }

    ilength = 1.0 / length;
    out[0]  = in[0] * ilength;
    out[1]  = in[1] * ilength;
    out[2]  = in[2] * ilength;

    return length;
}

static inline vec_t ColorNormalize(vec3_t in, vec3_t out) {
    float max, scale;

    max = in[0];
    if (in[1] > max)
        max = in[1];
    if (in[2] > max)
        max = in[2];

    if (max == 0)
        return 0;

    scale = 1.0 / max;

    VectorScale(in, scale, out);

    return max;
}

static inline void ClearBounds(vec3_t mins, vec3_t maxs) {
    mins[0] = mins[1] = mins[2] = BOGUS_RANGE;
    maxs[0] = maxs[1] = maxs[2] = -BOGUS_RANGE;
}

static inline void AddPointToBounds(vec3_t v, vec3_t mins, vec3_t maxs) {
    int32_t i;
    vec_t val;

    for (i = 0; i < 3; i++) {
        val = v[i];
        if (val < mins[i])
            mins[i] = val;
        if (val > maxs[i])
            maxs[i] = val;
    }
}

static inline qboolean RayPlaneIntersect(vec3_t p_n, vec_t p_d, vec3_t l_o, vec3_t l_n,
                                         vec3_t res) {
    float dot, t;

    dot = DotProduct(p_n, l_n);

    if (dot > -0.001)
        return false;

    t      = (p_d - (l_o[0] * p_n[0]) - (l_o[1] * p_n[1]) - (l_o[2] * p_n[2])) / dot;

    res[0] = l_o[0] + (t * l_n[0]);
    res[1] = l_o[1] + (t * l_n[1]);
    res[2] = l_o[2] + (t * l_n[2]);

    return true;
}

static inline void AngleVector(vec3_t angles, vec3_t out) // jit
{
    float angle;
    static float sp, sy, cp, cy;
    // static to help MS compiler fp bugs

    angle = angles[YAW] * ((float)Q_PI * 2.0f / 360.0f);
    sy    = sin(angle);
    cy    = cos(angle);
    angle = angles[PITCH] * ((float)Q_PI * 2.0f / 360.0f);
    sp    = sin(angle);
    cp    = cos(angle);

    if (out) {
        out[0] = cp * cy;
        out[1] = cp * sy;
        out[2] = -sp;
    }
}
#endif

#endif
