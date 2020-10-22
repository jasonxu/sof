// SPDX-License-Identifier: BSD-3-Clause
//
// Copyright(c) 2020 Google LLC. All rights reserved.
//
// Author: Pin-chih Lin <johnylin@google.com>

#include <sof/audio/drc/drc_math.h>
#include <sof/audio/format.h>
#include <sof/math/decibels.h>
#include <sof/math/numbers.h>
#include <math.h>

inline int32_t decibels_to_linear(int32_t decibels)
{
/*
 * Input is Q8.24: max 128.0
 * Output is Q12.20: max 2048.0
 */
	return db2lin_fixed(decibels);
}

static inline int32_t warp_rexp(int32_t x, int precision_x, int *e)
{
/*
 * Input depends on precision_x
 * Output range [0.5, 1); regulated to Q2.30
 */
	int bit;
	int32_t bitmask = 0x40000000;

	for (bit = 31; bit > 0; bit--) {
		if (x & bitmask)
			break;
		bitmask >>= 1;
	}

	*e = bit - precision_x;

	if (bit > 30)
		return Q_SHIFT_RND(x, bit, 30);
	if (bit < 30)
		return Q_SHIFT_LEFT(x, bit, 30);
	return x;
}

static inline int32_t warp_log10(int32_t x)
{
/*
 * Input is Q6.26: max 32.0
 * Output range ~ (-inf, 1.505); regulated to Q6.26: (-32.0, 32.0)
 */
#define q_v 26
#define q_mult(a, b, qa, qb, qy) ((int32_t)Q_MULTSR_32X32((int64_t)a, b, qa, qb, qy))
	/* Coefficients obtained from:
	 * fpminimax(log10(x), 5, [|SG...|], [1/2;sqrt(2)/2], absolute);
	 * max err ~= 6.088e-8
	 */
	const int32_t ONE_OVER_SQRT2 = Q_CONVERT_FLOAT(0.70710678118654752f, 30); /* 1/sqrt(2) */
	const int32_t A5 = Q_CONVERT_FLOAT(1.131880283355712890625f, q_v);
	const int32_t A4 = Q_CONVERT_FLOAT(-4.258677959442138671875f, q_v);
	const int32_t A3 = Q_CONVERT_FLOAT(6.81631565093994140625f, q_v);
	const int32_t A2 = Q_CONVERT_FLOAT(-6.1185703277587890625f, q_v);
	const int32_t A1 = Q_CONVERT_FLOAT(3.6505267620086669921875f, q_v);
	const int32_t A0 = Q_CONVERT_FLOAT(-1.217894077301025390625f, q_v);
	const int32_t LOG10_2 = Q_CONVERT_FLOAT(0.301029995663981195214f, q_v);
	int e;
	int32_t exp; /* Q31.1 */
	int32_t x2, x4; /* Q2.30 */
	int32_t A5Xx, A3Xx;

	x = warp_rexp(x, 26, &e); /* Q2.30 */
	exp = (int32_t)e << 1; /* Q_CONVERT_FLOAT(e, 1) */

	if (x > ONE_OVER_SQRT2) {
		x = q_mult(x, ONE_OVER_SQRT2, 30, 30, 30);
		exp += 1; /* Q_CONVERT_FLOAT(0.5, 1) */
	}

	x2 = q_mult(x, x, 30, 30, 30);
	x4 = q_mult(x2, x2, 30, 30, 30);
	A5Xx = q_mult(A5, x, q_v, 30, q_v);
	A3Xx = q_mult(A3, x, q_v, 30, q_v);
	return q_mult((A5Xx + A4), x4, q_v, 30, q_v) + q_mult((A3Xx + A2), x2, q_v, 30, q_v)
		+ q_mult(A1, x, q_v, 30, q_v) + A0 + q_mult(exp, LOG10_2, 1, q_v, q_v);
#undef q_mult
#undef q_v
}

inline int32_t linear_to_decibels(int32_t linear)
{
/*
 * Input is Q6.26: max 32.0
 * Output range ~ (-inf, 30.1030); regulated to Q11.21: (-1024.0, 1024.0)
 */
	/* For negative or zero, just return a very small dB value. */
	if (linear <= 0)
		return Q_CONVERT_FLOAT(-1000.0f, 21);

	int32_t log10_linear = warp_log10(linear); /* Q6.26 */
	return Q_MULTSR_32X32((int64_t)20, log10_linear, 0, 26, 21);
}

inline int32_t warp_log(int32_t x)
{
/*
 * Input is Q6.26: max 32.0
 * Output range ~ (-inf, 3.4657); regulated to Q6.26: (-32.0, 32.0)
 */
	if (x <= 0)
		return Q_CONVERT_FLOAT(-30.0f, 26);

	/* log(x) = log(10) * log10(x) */
	const int32_t LOG10 = Q_CONVERT_FLOAT(2.3025850929940457f, 29);
	int32_t log10_x = warp_log10(x); /* Q6.26 */
	return Q_MULTSR_32X32((int64_t)LOG10, log10_x, 29, 26, 26);
}

inline int32_t warp_sin(int32_t x)
{
/*
 * Input is Q2.30: (-2.0, 2.0)
 * Output range: [-1.0, 1.0]; regulated to Q2.30: (-2.0, 2.0)
 */
#define q_v 30
#define q_multv(a, b) ((int32_t)Q_MULTSR_32X32((int64_t)a, b, q_v, q_v, q_v))
	/* Coefficients obtained from:
	 * fpminimax(sin(x*pi/2), [|1,3,5,7|], [|SG...|], [-1e-30;1], absolute)
	 * max err ~= 5.901e-7
	 */
	const int32_t A7 = Q_CONVERT_FLOAT(-4.3330336920917034149169921875e-3f, q_v);
	const int32_t A5 = Q_CONVERT_FLOAT(7.9434238374233245849609375e-2f, q_v);
	const int32_t A3 = Q_CONVERT_FLOAT(-0.645892798900604248046875f, q_v);
	const int32_t A1 = Q_CONVERT_FLOAT(1.5707910060882568359375f, q_v);
	int32_t x2 = q_multv(x, x);
	int32_t x4 = q_multv(x2, x2);

	int32_t A3Xx2 = q_multv(A3, x2);
	int32_t A7Xx2 = q_multv(A7, x2);

	return q_multv(x, (q_multv(x4, (A7Xx2 + A5)) + A3Xx2 + A1));
#undef q_multv
#undef q_v
}

inline int32_t warp_asin(int32_t x)
{
/*
 * Input is Q2.30: (-2.0, 2.0)
 * Output range: [-1.0, 1.0]; regulated to Q2.30: (-2.0, 2.0)
 */
#define q_vl 30
#define q_vh 26
#define q_multv(a, b, q) ((int32_t)Q_MULTSR_32X32((int64_t)a, b, q, q, q))
	/* Coefficients obtained from:
	 * If x <= 1/sqrt(2), then
	 *   fpminimax(asin(x), [|1,3,5,7|], [|SG...|], [-1e-30;1/sqrt(2)], absolute)
	 *   max err ~= 1.89936e-5
	 * Else then
	 *   fpminimax(asin(x), [|1,3,5,7|], [|SG...|], [1/sqrt(2);1], absolute)
	 *   max err ~= 3.085226e-2
	 */
	const int32_t TWO_OVER_PI = Q_CONVERT_FLOAT(0.63661977236758134f, q_vl); /* 2/pi */
	const int32_t ONE_OVER_SQRT2 = Q_CONVERT_FLOAT(0.70710678118654752f, q_vl); /* 1/sqrt(2) */
	const int32_t A7L = Q_CONVERT_FLOAT(0.1181826665997505187988281f, q_vl);
	const int32_t A5L = Q_CONVERT_FLOAT(4.0224377065896987915039062e-2f, q_vl);
	const int32_t A3L = Q_CONVERT_FLOAT(0.1721895635128021240234375f, q_vl);
	const int32_t A1L = Q_CONVERT_FLOAT(0.99977016448974609375f, q_vl);

	const int32_t A7H = Q_CONVERT_FLOAT(14.12774658203125f, q_vh);
	const int32_t A5H = Q_CONVERT_FLOAT(-30.1692714691162109375f, q_vh);
	const int32_t A3H = Q_CONVERT_FLOAT(21.4760608673095703125f, q_vh);
	const int32_t A1H = Q_CONVERT_FLOAT(-3.894591808319091796875f, q_vh);

	int32_t A7, A5, A3, A1, q_v;
	int32_t x2, x4;
	int32_t A3Xx2, A7Xx2, asinx;

	if (ABS(x) <= ONE_OVER_SQRT2) {
		A7 = A7L;
		A5 = A5L;
		A3 = A3L;
		A1 = A1L;
		q_v = q_vl;
	} else {
		A7 = A7H;
		A5 = A5H;
		A3 = A3H;
		A1 = A1H;
		q_v = q_vh;
		x = Q_SHIFT_RND(x, q_vl, q_vh); /* Q6.26 */
	}

	x2 = q_multv(x, x, q_v);
	x4 = q_multv(x2, x2, q_v);

	A3Xx2 = q_multv(A3, x2, q_v);
	A7Xx2 = q_multv(A7, x2, q_v);

	asinx = q_multv(x, (q_multv(x4, (A7Xx2 + A5), q_v) + A3Xx2 + A1), q_v);
	return Q_MULTSR_32X32((int64_t)asinx, TWO_OVER_PI, q_v, q_vl, 30);
#undef q_multv
#undef q_vh
#undef q_vl
}

inline int32_t warp_pow(int32_t x, int32_t y)
{
/*
 * Input x is Q6.26: (-32.0, 32.0)
 *       y is Q2.30: (-2.0, 2.0)
 * Output is Q12.20: max 2048.0
 */
	/* x^y = expf(y * log(x)) */
	return exp_fixed(Q_MULTSR_32X32((int64_t)y, warp_log(x), 30, 26, 27));
}

inline int32_t warp_inv(int32_t x, int precision_x, int precision_y)
{
#define q_v 25
#define q_mult(a, b, qa, qb, qy) ((int32_t)Q_MULTSR_32X32((int64_t)a, b, qa, qb, qy))
	/* Coefficients obtained from:
	 * fpminimax(1/x, 5, [|SG...|], [sqrt(2)/2;1], absolute);
	 * max err ~= 1.00388e-6
	 */
	const int32_t ONE_OVER_SQRT2 = Q_CONVERT_FLOAT(0.70710678118654752f, 30); /* 1/sqrt(2) */
	const int32_t SQRT2 = Q_CONVERT_FLOAT(1.4142135623730950488f, 30); /* sqrt(2) */
	const int32_t A5 = Q_CONVERT_FLOAT(-2.742647647857666015625f, q_v);
	const int32_t A4 = Q_CONVERT_FLOAT(14.01327800750732421875f, q_v);
	const int32_t A3 = Q_CONVERT_FLOAT(-29.74465179443359375f, q_v);
	const int32_t A2 = Q_CONVERT_FLOAT(33.57208251953125f, q_v);
	const int32_t A1 = Q_CONVERT_FLOAT(-21.25031280517578125f, q_v);
	const int32_t A0 = Q_CONVERT_FLOAT(7.152250766754150390625f, q_v);
	int e;
	int sqrt2_extracted = 0;
	int32_t x2, x4; /* Q2.30 */
	int32_t A5Xx, A3Xx;
	int32_t inv;

	x = warp_rexp(x, precision_x, &e); /* Q2.30 */

	if (x < ONE_OVER_SQRT2) {
		x = q_mult(x, SQRT2, 30, 30, 30);
		sqrt2_extracted = 1;
	}

	x2 = q_mult(x, x, 30, 30, 30);
	x4 = q_mult(x2, x2, 30, 30, 30);
	A5Xx = q_mult(A5, x, q_v, 30, q_v);
	A3Xx = q_mult(A3, x, q_v, 30, q_v);
	inv = q_mult((A5Xx + A4), x4, q_v, 30, q_v) + q_mult((A3Xx + A2), x2, q_v, 30, q_v)
		+ q_mult(A1, x, q_v, 30, q_v) + A0;

	if (sqrt2_extracted)
		inv = q_mult(inv, SQRT2, q_v, 30, q_v);

	e += q_v;
	if (e > precision_y)
		return Q_SHIFT_RND(inv, e, precision_y);
	if (e < precision_y)
		return Q_SHIFT_LEFT(inv, e, precision_y);
	return inv;
#undef q_mult
#undef q_v
}

inline int32_t knee_exp(int32_t input)
{
/*
 * Input is Q5.27: max 16.0
 * Output is Q12.20: max 2048.0
 */
	return exp_fixed(input);
}
