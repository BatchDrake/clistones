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
  <http://www.gnu.org/licenses/>

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef FILENAME
#  define FILENAME __FILENAME__
#endif /* FILENAME */

#include <graves.h>

void
graves_det_destroy(graves_det_t *detect)
{
  if (detect->q_hist != NULL)
    free(detect->q_hist);

  if (detect->p_n_hist != NULL)
    free(detect->p_n_hist);

  if (detect->p_w_hist != NULL)
    free(detect->p_w_hist);

  if (detect->samp_hist != NULL)
    free(detect->samp_hist);

  su_iir_filt_finalize(&detect->lpf1);
  su_iir_filt_finalize(&detect->lpf2);

  grow_buf_clear(&detect->chirp);
  grow_buf_clear(&detect->q);
  grow_buf_clear(&detect->p_n_buf);
  grow_buf_clear(&detect->p_w_buf);

  free(detect);
}

SUPRIVATE void
graves_det_filt_back(graves_det_t *md)
{
  int i;
  int shift = (int) md->hist_len;
  int len = (int) (grow_buf_get_size(&md->p_n_buf) / sizeof (SUFLOAT));
  SUFLOAT *p_n_ptr = grow_buf_get_buffer(&md->p_n_buf);
  SUFLOAT *p_w_ptr = grow_buf_get_buffer(&md->p_w_buf);
  SUFLOAT *q_ptr;
  SUFLOAT  p_n = md->p_n;
  SUFLOAT  p_w = md->p_w;

  grow_buf_alloc(&md->q, (unsigned) len * sizeof (SUFLOAT));
  q_ptr = grow_buf_get_buffer(&md->q);

  memset(q_ptr, 0, (unsigned) len * sizeof (SUFLOAT));

  /* Apply filters in reverse order. */
  for (i = len - 1; i >= 0; --i) {
    p_w += md->alpha * (p_w_ptr[i] - p_w);
    p_n += md->alpha * (p_n_ptr[i] - p_n);

    p_n_ptr[i] = p_n;
    p_w_ptr[i] = p_w;
  }

  /* Shift */
  for (i = 0; i < len - shift; ++i) {
    p_n_ptr[i] = p_n_ptr[i + shift];
    p_w_ptr[i] = p_w_ptr[i + shift];
    q_ptr[i]   = p_n_ptr[i] / p_w_ptr[i];
  }
}

SUBOOL
graves_det_feed(graves_det_t *md, SUCOMPLEX x)
{
  SUCOMPLEX y;
  SUFLOAT   Q;
  SUFLOAT   energy;
  struct graves_chirp_info info;
  unsigned int i;

  x *= SU_C_CONJ(su_ncqo_read(&md->lo));

  y = su_iir_filt_feed(&md->lpf1, x);
  md->p_w += md->alpha * (SU_C_REAL(y * SU_C_CONJ(y)) - md->p_w);

  y = su_iir_filt_feed(&md->lpf2, x);
  md->p_n += md->alpha * (SU_C_REAL(y * SU_C_CONJ(y)) - md->p_n);

  /* Compute power quotient */
  Q = md->p_n / md->p_w;

  if (Q >= 1 || Q < md->ratio)
    Q = md->last_good_q;
  else
    md->last_good_q = Q;

  /* Update histories */
  md->p_n_hist[md->p]  = md->p_n;
  md->p_w_hist[md->p]  = md->p_w;
  md->q_hist[md->p]    = Q;
  md->samp_hist[md->p] = y;

  if (++md->p == md->hist_len)
    md->p = 0;

  /* md->p now points to the OLDEST sample */

  /* Compute cross-correlation */
  energy = 0;
  for (i = 0; i < md->hist_len; ++i)
    energy += md->q_hist[i];

  /* Detect chirp limits */
  if (md->in_chirp) {
    if (energy < md->energy_thres) {
      /* DETECTED: CHIRP END */
      md->in_chirp = SU_FALSE;

      graves_det_filt_back(md);

      info.length = (unsigned int) grow_buf_get_size(&md->chirp) / sizeof(SUCOMPLEX);
      info.length -= md->hist_len;

      if (info.length > 0) {
        info.t0     = (md->n - info.length) / md->params.fs;
        info.t0f    = SU_ASFLOAT((md->n - info.length) % md->params.fs) / md->params.fs;
        info.x      = (const SUCOMPLEX *) grow_buf_get_buffer(&md->chirp);
        info.q      = (const SUFLOAT *) grow_buf_get_buffer(&md->q);
        info.p_n    = (const SUFLOAT *) grow_buf_get_buffer(&md->p_n_buf);
        info.p_w    = (const SUFLOAT *) grow_buf_get_buffer(&md->p_w_buf);

        info.fs     = md->params.fs;
        info.rbw    = md->ratio;

        SU_TRYCATCH((md->on_chirp) (md->privdata, &info), return SU_FALSE);
      }
#ifdef DEBUG
      printf(
          "Chirp of length %5d detected (at %02d:%02d:%02d)\n",
          grow_buf_get_size(&md->chirp) / sizeof(SUCOMPLEX),
          start / 3600,
          (start / 60) % 60,
          start % 60);
#endif
    } else {
      /* Sample belongs to chirp. Save it for later processing */
      SU_TRYCATCH(
          grow_buf_append(&md->chirp, &y, sizeof(SUCOMPLEX)) != -1,
          return SU_FALSE)
      SU_TRYCATCH(
          grow_buf_append(&md->p_n_buf, &md->p_n, sizeof(SUFLOAT)) != -1,
          return SU_FALSE)
      SU_TRYCATCH(
          grow_buf_append(&md->p_w_buf, &md->p_w, sizeof(SUFLOAT)) != -1,
          return SU_FALSE)
    }
  } else {
    if (energy >= md->energy_thres) {
      /* DETECTED: CHIRP START */
      md->in_chirp = SU_TRUE;

      /* Save all samples in the delay line to the grow buffer */
      grow_buf_shrink(&md->chirp);
      grow_buf_shrink(&md->q);
      grow_buf_shrink(&md->p_n_buf);
      grow_buf_shrink(&md->p_w_buf);

      for (i = 0; i < md->hist_len; ++i) {
        SU_TRYCATCH(
            grow_buf_append(
                &md->chirp,
                md->samp_hist + (i + md->p) % md->hist_len,
                sizeof(SUCOMPLEX)) != -1,
            return SU_FALSE)
        SU_TRYCATCH(
            grow_buf_append(
                &md->p_n_buf,
                md->p_n_hist + (i + md->p) % md->hist_len,
                sizeof(SUCOMPLEX)) != -1,
            return SU_FALSE)
        SU_TRYCATCH(
            grow_buf_append(
                &md->p_w_buf,
                md->p_w_hist + (i + md->p) % md->hist_len,
                sizeof(SUCOMPLEX)) != -1,
            return SU_FALSE)
      }
    }
  }

  ++md->n;
  return SU_TRUE;
}

void
graves_det_set_center_freq(graves_det_t *md, SUFLOAT fc)
{
  su_ncqo_set_freq(
        &md->lo,
        SU_ABS2NORM_FREQ(md->params.fs, fc));
}

SUPRIVATE SUBOOL
graves_det_check_params(const struct graves_det_params *params)
{
  if (params->lpf1 <= params->lpf2) {
    SU_ERROR("Illegal filter cutoff frequencies (lpf1 < lpf2)\n");
    return SU_FALSE;
  }

  if (SU_ABS2NORM_FREQ(
        params->fs,
        params->lpf1) < GRAVES_MIN_LPF_CUTOFF) {
    SU_ERROR(
          "LPF1 is too narrow (safe minimum is %g Hz)",
          SU_NORM2ABS_FREQ(params->fs, GRAVES_MIN_LPF_CUTOFF));
    return SU_FALSE;
  }

  if (SU_ABS2NORM_FREQ(
        params->fs,
        params->lpf2) < GRAVES_MIN_LPF_CUTOFF) {
    SU_ERROR(
          "LPF2 is too narrow (safe minimum is %g Hz)",
          SU_NORM2ABS_FREQ(params->fs, GRAVES_MIN_LPF_CUTOFF));
    return SU_FALSE;
  }

  return SU_TRUE;
}

graves_det_t *
graves_det_new(
    const struct graves_det_params *params,
    graves_chirp_cb_t chrp_fn,
    void *privdata)
{
  graves_det_t *new = NULL;

  if (!graves_det_check_params(params))
    return NULL;

  SU_TRYCATCH(new = calloc(1, sizeof (graves_det_t)), goto fail)

  new->params = *params;
  new->ratio  = params->lpf2 / params->lpf1;
  new->alpha = 1 - SU_EXP(-SU_ADDSFX(1.) / (params->fs * MIN_CHIRP_DURATION));
  new->on_chirp = chrp_fn;
  new->privdata = privdata;

  su_ncqo_init(&new->lo, SU_ABS2NORM_FREQ(params->fs, params->fc));

  SU_TRYCATCH(
      su_iir_bwlpf_init(
          &new->lpf1,
          4,
          SU_ABS2NORM_FREQ(params->fs, params->lpf1)),
      goto fail)

  SU_TRYCATCH(
      su_iir_bwlpf_init(
          &new->lpf2,
          4,
          SU_ABS2NORM_FREQ(params->fs, params->lpf2)),
      goto fail)

# if 0
  int i = 0;
  SUCOMPLEX y1, y2;
  SUCOMPLEX e1, e2;

  y1 = su_iir_filt_feed(&new->lpf1, 1);
  y2 = su_iir_filt_feed(&new->lpf2, 1);

  e1 = SU_C_CONJ(y1) * y1;
  e2 = SU_C_CONJ(y2) * y2;

  for (i = 0; i < 1000; ++i) {
    printf("Ratio: %g, measured: %g\n", new->ratio, SU_C_REAL(e2 / e1));

    y1 = su_iir_filt_feed(&new->lpf1, 0);
    y2 = su_iir_filt_feed(&new->lpf2, 0);

    e1 += SU_C_CONJ(y1) * y1;
    e2 += SU_C_CONJ(y2) * y2;
  }
#endif

  new->hist_len = (SUSCOUNT) (SU_CEIL(params->fs * MIN_CHIRP_DURATION));
  new->energy_thres = params->threshold * new->ratio * new->hist_len;

  SU_TRYCATCH(
      new->q_hist   = calloc(sizeof(SUFLOAT), new->hist_len),
      goto fail)

  SU_TRYCATCH(
      new->samp_hist = calloc(sizeof(SUCOMPLEX), new->hist_len),
      goto fail)

  SU_TRYCATCH(
      new->p_n_hist = calloc(sizeof(SUFLOAT), new->hist_len),
      goto fail)

  SU_TRYCATCH(
      new->p_w_hist = calloc(sizeof(SUFLOAT), new->hist_len),
      goto fail)
  return new;

fail:
  if (new != NULL)
    graves_det_destroy(new);

  return new;
}
