/* SPDX-License-Identifier: BSD-3-Clause
 *
 * Copyright(c) 2020 Google LLC. All rights reserved.
 *
 * Author: Pin-chih Lin <johnylin@google.com>
 */
#ifndef __SOF_AUDIO_DRC_DRC_MATH_H__
#define __SOF_AUDIO_DRC_DRC_MATH_H__

#include <stddef.h>
#include <stdint.h>

#define DRC_NEG_TWO_DB 0.7943282347242815f /* -2dB = 10^(-2/20) */

int32_t decibels_to_linear(int32_t decibels); /* Input:Q8.24 Output:Q12.20 */
int32_t linear_to_decibels(int32_t linear); /* Input:Q6.26 Output:Q11.21 */
int32_t warp_log(int32_t x); /* Input:Q6.26 Output:Q6.26 */
int32_t warp_sin(int32_t x); /* Input:Q2.30 Output:Q2.30 */
int32_t warp_asin(int32_t x); /* Input:Q2.30 Output:Q2.30 */
int32_t warp_pow(int32_t x, int32_t y); /* Input:Q6.26, Q2.30 Output:Q12.20 */
int32_t warp_inv(int32_t x, int precision_x, int precision_y);
int32_t knee_exp(int32_t input); /* Input:Q5.27 Output:Q12.20 */

#endif //  __SOF_AUDIO_DRC_DRC_MATH_H__
