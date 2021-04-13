/*

  Copyright (C) 2018 Gonzalo José Carracedo Carballal

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as
  published by the Free Software Foundation, either version 3 of the
  License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this program.  If not, see
  <http: *www.gnu.org/licenses/>

*/

#ifndef GRAVES_GRAVES_H
#define GRAVES_GRAVES_H

#include <util/util.h>

#include <sigutils/iir.h>
#include <sigutils/ncqo.h>
#include <sigutils/log.h>
#include <sigutils/sampling.h>

#ifdef __cplusplus
extern "C" {
#endif

#define GRAVES_CENTER_FREQ SU_ADDSFX(143050000.)
#define SPEED_OF_LIGHT     SU_ADDSFX(299792458.)

#if 0
#define LPF1_CUTOFF SU_ADDSFX(300.) /* In Hz */
#define LPF2_CUTOFF SU_ADDSFX(50.)  /* In Hz */
#define GRAVES_POWER_RATIO (LPF2_CUTOFF / LPF1_CUTOFF)
#define GRAVES_POWER_THRESHOLD (SU_ASFLOAT(2.) * GRAVES_POWER_RATIO)
#endif

#define GRAVES_MIN_LPF_CUTOFF SU_ABS2NORM_FREQ(8000, 50)

#define MIN_CHIRP_DURATION SU_ADDSFX(0.07)

struct graves_chirp_info {
  SUSCOUNT t0;  /* Start time */
  SUFLOAT t0f;  /* Decimal part of the start time */

  SUSCOUNT fs;
  SUFLOAT  rbw; /* Bandwidth ratio */

  /* Unsigned int length */
  unsigned int length;

  /* Chirp data */
  const SUCOMPLEX *x;

  /* Quotient data */
  const SUFLOAT   *q;

  /* Narrow channel power data */
  const SUFLOAT   *p_n;

  /* Wide channel power data */
  const SUFLOAT   *p_w;
};

typedef SUBOOL (*graves_chirp_cb_t) (
    void *privdata,
    const struct graves_chirp_info *info);

struct graves_det_params {
  SUSCOUNT fs;
  SUFLOAT  fc;
  SUFLOAT  lpf1;
  SUFLOAT  lpf2;
  SUFLOAT  threshold;
};

#define graves_det_params_INITIALIZER \
{                                     \
  8000,             /* fs */          \
  SU_ADDSFX(1000.), /* fc */          \
  SU_ADDSFX(300.),  /* lpf1 */        \
  SU_ADDSFX(50.),   /* lpf2 */        \
  SU_ADDSFX(2.),    /* threshoid */   \
}

struct graves_det {
  struct graves_det_params params;
  SUFLOAT ratio;
  SUSCOUNT n;          /* Samples consumed */
  su_iir_filt_t lpf1; /* Low pass filter 1. Used to detect noise power */
  su_iir_filt_t lpf2; /* Low pass filter 2. Used to isolate chirps */
  su_ncqo_t lo;
  SUFLOAT alpha; /* Slow decay, used to detect chirps */
  SUFLOAT last_good_q;
  SUFLOAT p_w; /* Wide channel power */
  SUFLOAT p_n; /* Narrow channel power */

  SUSCOUNT hist_len;
  SUSCOUNT p;
  SUFLOAT   *p_n_hist;
  SUFLOAT   *p_w_hist;
  SUFLOAT   *q_hist;
  SUCOMPLEX *samp_hist;

  SUFLOAT   energy_thres;
  SUBOOL    in_chirp;

  grow_buf_t chirp;
  grow_buf_t q;
  grow_buf_t p_n_buf;
  grow_buf_t p_w_buf;

  void *privdata;

  graves_chirp_cb_t on_chirp;
};

typedef struct graves_det graves_det_t;

/*
 * Q is actually the quotient of the averaged power at the output of
 * two filters: P_n / P_w
 *
 * When a signal is present, the power can be expressed as the sum of 2
 * contributions:
 *
 *   P_n = W_n * N + S
 *   P_w = W_w * N + S
 *
 * Where N is the PSD of the noise and W_n, W_w are the bandwidths of the
 * narrow and wide low pass filters respectively. We want to deduce S / (W_w * N),
 * therefore:
 *
 * P_n - W_n * N = P_w - W_w * N -> (P_n - P_w) / (W_n - W_w) = N
 *
 * And replacing:
 *
 * S = P_n - W_n * (P_n - P_w) / (W_n - W_w)
 *
 * And finally :
 *
 * S / N = (P_n - W_n * (P_n - P_w) / (W_n - W_w))  / ((P_n - P_w) / (W_n - W_w))
 *
 */


SUINLINE SUFLOAT
graves_det_q_to_snr(SUFLOAT ratio, SUFLOAT q)
{
  return (q - ratio) / (SU_ADDSFX(1.) - q);
}

SUINLINE SUFLOAT
graves_det_get_N0(SUFLOAT ratio, SUFLOAT p_n, SUFLOAT snr)
{
  return p_n / (ratio + snr);
}

SUINLINE SUFLOAT
graves_det_get_ratio(const graves_det_t *det)
{
  return det->ratio;
}

SUINLINE const struct graves_det_params *
graves_det_get_params(const graves_det_t *det)
{
  return &det->params;
}

void graves_det_destroy(graves_det_t *detect);

void graves_det_set_center_freq(graves_det_t *md, SUFLOAT fc);

SUBOOL graves_det_feed(graves_det_t *md, SUCOMPLEX x);

graves_det_t *
graves_det_new(
    const struct graves_det_params *params,
    graves_chirp_cb_t chrp_fn,
    void *privdata);

#ifdef __cplusplus
}
#endif


#endif /* _GRAVES_H */

