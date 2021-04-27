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

SUPRIVATE void
graves_det_element_finalize(struct graves_det_element *self)
{
  if (self->q_hist != NULL)
    free(self->q_hist);

  if (self->p_n_hist != NULL)
    free(self->p_n_hist);

  if (self->p_w_hist != NULL)
    free(self->p_w_hist);

  if (self->samp_hist != NULL)
    free(self->samp_hist);

  if (self->pres_hist != NULL)
    free(self->pres_hist);

  su_iir_filt_finalize(&self->lpf1);
  su_iir_filt_finalize(&self->lpf2);
}


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
SUPRIVATE SUBOOL
graves_det_element_feed(struct graves_det_element *self, SUCOMPLEX x)
{
  SUCOMPLEX y;
  SUFLOAT   Q;
  SUFLOAT   energy;
  struct graves_chirp_info info;
  unsigned int i;

  y = su_iir_filt_feed(&self->lpf1, x);
  self->p_w += self->alpha * (SU_C_REAL(y * SU_C_CONJ(y)) - self->p_w);

  y = su_iir_filt_feed(&self->lpf2, x);
  self->p_n += self->alpha * (SU_C_REAL(y * SU_C_CONJ(y)) - self->p_n);

  /* Compute power quotient */
  Q = self->p_n / self->p_w;

  if (Q >= 1 || Q < self->ratio)
    Q = self->last_good_q;
  else
    self->last_good_q = Q;

  /* Update histories */
  self->p_n_hist[self->p]  = self->p_n;
  self->p_w_hist[self->p]  = self->p_w;
  self->q_hist[self->p]    = Q;

  /* Compute cross-correlation */
  energy = 0;
  for (i = 0; i < self->hist_len; ++i)
    energy += self->q_hist[i];

  self->present = energy >= self->energy_thres;
  self->y = y;

  self->samp_hist[self->p] = self->y;
  self->pres_hist[self->p] = self->present;

  if (++self->p == self->hist_len)
    self->p = 0;

  return self->present;
}

SUPRIVATE SUBOOL
graves_det_element_init(
    struct graves_det_element *self,
    const struct graves_det_params *params)
{
  memset(self, 0, sizeof (struct graves_det_element));

  SU_TRYCATCH(
      su_iir_bwlpf_init(
          &self->lpf1,
          4,
          SU_ABS2NORM_FREQ(params->fs, params->lpf1)),
      goto fail)

  SU_TRYCATCH(
      su_iir_bwlpf_init(
          &self->lpf2,
          4,
          SU_ABS2NORM_FREQ(params->fs, params->lpf2)),
      goto fail)

  self->ratio  = params->lpf2 / params->lpf1;
  self->hist_len = (SUSCOUNT) (SU_CEIL(params->fs * MIN_CHIRP_DURATION));
  self->energy_thres = params->threshold * self->ratio * self->hist_len;
  self->alpha = 1 - SU_EXP(-SU_ADDSFX(1.) / (params->fs * MIN_CHIRP_DURATION));

  SU_TRYCATCH(
      self->q_hist   = calloc(sizeof(SUFLOAT), self->hist_len),
      goto fail)

  SU_TRYCATCH(
      self->samp_hist = calloc(sizeof(SUCOMPLEX), self->hist_len),
      goto fail)

  SU_TRYCATCH(
      self->p_n_hist = calloc(sizeof(SUFLOAT), self->hist_len),
      goto fail)

  SU_TRYCATCH(
      self->p_w_hist = calloc(sizeof(SUFLOAT), self->hist_len),
      goto fail)

  SU_TRYCATCH(
      self->p_w_hist = calloc(sizeof(SUFLOAT), self->hist_len),
      goto fail)

  SU_TRYCATCH(
      self->pres_hist = calloc(sizeof(SUBOOL), self->hist_len),
      goto fail)

  return SU_TRUE;

fail:
  graves_det_element_finalize(self);

  return SU_FALSE;
}

void
graves_det_destroy(graves_det_t *self)
{
  int i;

  if (self->det_bank != NULL) {
    for (i = 0; i < self->params.multiplicity; ++i)
      graves_det_element_finalize(self->det_bank + i);

    free(self->det_bank);
  }

  if (self->mixer_hist)
    free(self->mixer_hist);

  grow_buf_clear(&self->chirp);
  grow_buf_clear(&self->snr_buf);
  grow_buf_clear(&self->S_buf);
  grow_buf_clear(&self->N_buf);

  free(self);
}

SUINLINE void
graves_det_effective_signal_and_noise_at(
    const graves_det_t *self,
    int at,
    SUFLOAT *S,
    SUFLOAT *N)
{
  SUFLOAT S_power = 0;
  SUFLOAT N_power = 0;
  SUFLOAT curr_S;

  SUFLOAT w_w = self->params.lpf1;
  SUFLOAT w_n = self->params.lpf2;
  SUFLOAT P_w, P_n;
  int n;

  for (n = 0; n < self->params.multiplicity; ++n) {
    P_n = self->det_bank[n].p_n_hist[at];
    P_w = self->det_bank[n].p_w_hist[at];

    curr_S = P_n - w_n * (P_w - P_n) / (w_w - w_n);

    if (self->det_bank[n].pres_hist[at])
      S_power += curr_S;

    N_power += P_n - curr_S;
  }

  /* The effective noise level is measured in a single filter */

  *S = S_power;
  *N = N_power / self->params.multiplicity;
}

SUINLINE void
graves_det_effective_signal_and_noise(
    const graves_det_t *self,
    SUFLOAT *S,
    SUFLOAT *N)
{
  int hist_len = self->det_bank[0].hist_len;
  int p = (hist_len + self->det_bank[0].p - 1) % hist_len;

  return graves_det_effective_signal_and_noise_at(self, p, S, N);
}

/*
 * graves_det_filt_back computes the narrow and wide power outputs according
 * to an effective SNR based on the smart signal isolation using multiple
 * filters.
 */
SUPRIVATE void
graves_det_filt_back(graves_det_t *self)
{
  /* We can do this because all history lengths and alphas are the same */
  int i;
  int shift = (int) self->det_bank[0].hist_len;
  SUFLOAT alpha = self->det_bank[0].alpha;
  SUFLOAT ratio = self->det_bank[0].ratio;
  int len = (int) (grow_buf_get_size(&self->S_buf) / sizeof (SUFLOAT));

  SUFLOAT *S_ptr = grow_buf_get_buffer(&self->S_buf);
  SUFLOAT *N_ptr = grow_buf_get_buffer(&self->N_buf);
  SUFLOAT *snr_ptr;
  SUFLOAT  S = 0;
  SUFLOAT  N = 0;

  graves_det_effective_signal_and_noise(self, &S, &N);

  grow_buf_alloc(&self->snr_buf, (unsigned) len * sizeof (SUFLOAT));
  snr_ptr = grow_buf_get_buffer(&self->snr_buf);

  memset(snr_ptr, 0, (unsigned) len * sizeof (SUFLOAT));

  /* Apply filters in reverse order. */
  for (i = len - 1; i >= 0; --i) {
    S += alpha * (S_ptr[i] - S);
    N += alpha * (N_ptr[i] - N);

    S_ptr[i] = S;
    N_ptr[i] = N;
  }

  /* Shift */
  for (i = 0; i < len - shift; ++i) {
    S_ptr[i] = S_ptr[i + shift];
    N_ptr[i] = N_ptr[i + shift];

    snr_ptr[i]   = S_ptr[i] / N_ptr[i];
  }
}

SUBOOL
graves_det_feed(graves_det_t *self, SUCOMPLEX x)
{
  SUCOMPLEX y;
  SUCOMPLEX m, curr_m, c;
  SUFLOAT S, N;
  SUBOOL any_present = SU_FALSE;
  SUBOOL detected;

  int hist_len = self->det_bank[0].hist_len;
  SUFLOAT ratio = self->det_bank[0].ratio;

  struct graves_chirp_info info;
  unsigned int i, n;
  int prev_i;

  m = su_ncqo_read(&self->mixer); /* Mixer sample */
  x *= SU_C_CONJ(su_ncqo_read(&self->lo));

  for (n = 0; n < self->params.multiplicity; ++n) {
    detected = graves_det_element_feed(&self->det_bank[n], x);
    any_present = any_present || detected;
    x *= SU_C_CONJ(m);
  }

  /* Detect chirp limits */
  if (self->in_chirp) {
    if (!any_present) {
      /* DETECTED: CHIRP END */
      self->in_chirp = SU_FALSE;

      graves_det_filt_back(self);

      info.length = (unsigned int) grow_buf_get_size(&self->chirp) / sizeof(SUCOMPLEX);
      info.length -= hist_len;

      if (info.length > 0) {
        info.t0     = (self->n - info.length) / self->params.fs;
        info.t0f    = SU_ASFLOAT((self->n - info.length) % self->params.fs) / self->params.fs;
        info.x      = (const SUCOMPLEX *) grow_buf_get_buffer(&self->chirp);
        info.snr    = (const SUFLOAT *) grow_buf_get_buffer(&self->snr_buf);
        info.S      = (const SUFLOAT *) grow_buf_get_buffer(&self->S_buf);
        info.N      = (const SUFLOAT *) grow_buf_get_buffer(&self->N_buf);

        info.fs     = self->params.fs;
        info.rbw    = ratio;

        SU_TRYCATCH((self->on_chirp) (self->privdata, &info), return SU_FALSE);
      }
#ifdef DEBUG
      printf(
          "Chirp of length %5d detected\n",
          grow_buf_get_size(&self->chirp) / sizeof(SUCOMPLEX));
#endif
    } else {

      /* Sample belongs to chirp. Save it for later processing */
      /* We include multiple chirps at a time according to their presence */
      graves_det_effective_signal_and_noise(self, &S, &N);

      y = 0;
      curr_m = SU_C_CONJ(su_ncqo_read(&self->center));

      for (n = 0; n < self->params.multiplicity; ++n) {
        if (self->det_bank[n].present)
          y += curr_m * self->det_bank[n].y;

        /* Increase frequency of this oscillator in m */
        curr_m *= m;
      }

      SU_TRYCATCH(
          grow_buf_append(&self->chirp, &y, sizeof(SUCOMPLEX)) != -1,
          return SU_FALSE)
      SU_TRYCATCH(
          grow_buf_append(&self->S_buf, &S, sizeof(SUFLOAT)) != -1,
          return SU_FALSE)
      SU_TRYCATCH(
          grow_buf_append(&self->N_buf, &N, sizeof(SUFLOAT)) != -1,
          return SU_FALSE)
    }
  } else {
    if (any_present) {
      /* DETECTED: CHIRP START */
      self->in_chirp = SU_TRUE;

      /* Save all samples in the delay line to the grow buffer */
      grow_buf_shrink(&self->chirp);
      grow_buf_shrink(&self->snr_buf);
      grow_buf_shrink(&self->S_buf);
      grow_buf_shrink(&self->N_buf);

      for (i = 0; i < hist_len; ++i) {
        prev_i = (i + self->det_bank[0].p) % hist_len;

        graves_det_effective_signal_and_noise_at(self, prev_i, &S, &N);

        curr_m = SU_C_CONJ(su_ncqo_read(&self->center));
        m = self->mixer_hist[prev_i];
        y = 0;
        for (n = 0; n < self->params.multiplicity; ++n) {
          if (self->det_bank[n].pres_hist[prev_i])
            y += curr_m * self->det_bank[n].samp_hist[prev_i];

          /* Increase frequency of this oscillator in m */
          curr_m *= m;
        }

        SU_TRYCATCH(
            grow_buf_append(&self->chirp, &y, sizeof(SUCOMPLEX)) != -1,
            return SU_FALSE)
        SU_TRYCATCH(
            grow_buf_append(&self->S_buf, &S, sizeof(SUCOMPLEX)) != -1,
            return SU_FALSE)
        SU_TRYCATCH(
            grow_buf_append(&self->N_buf, &N, sizeof(SUCOMPLEX)) != -1,
            return SU_FALSE)
      }
    }
  }

  ++self->n;
  return SU_TRUE;
}

void
graves_det_set_center_freq(graves_det_t *md, SUFLOAT fc)
{
  su_ncqo_set_freq(&md->lo, SU_ABS2NORM_FREQ(md->params.fs, fc));
}

SUPRIVATE SUBOOL
graves_det_check_params(const struct graves_det_params *params)
{
  if (params->multiplicity < 1) {
    SU_ERROR("Illegal multiplicity (at least one channel is required)\n");
    return SU_FALSE;
  }

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
  int i;
  SUFLOAT f0;

  if (!graves_det_check_params(params))
    return NULL;

  SU_TRYCATCH(new = calloc(1, sizeof (graves_det_t)), goto fail)

  new->params = *params;
  new->on_chirp = chrp_fn;
  new->privdata = privdata;

  f0 = params->fc - (params->multiplicity - 1) * params->lpf2;

  su_ncqo_init(
      &new->lo, SU_ABS2NORM_FREQ(params->fs, f0));

  su_ncqo_init(
      &new->mixer,
      SU_ABS2NORM_FREQ(params->fs, 2 * params->lpf2));

  su_ncqo_init(
      &new->center,
      SU_ABS2NORM_FREQ(params->fs, (params->multiplicity - 1) * params->lpf2));

  SU_TRYCATCH(
      new->det_bank = calloc(
          params->multiplicity,
          sizeof(struct graves_det_element)),
      goto fail);

  for (i = 0; i < params->multiplicity; ++i)
    SU_TRYCATCH(
        graves_det_element_init(&new->det_bank[i], params),
        goto fail);

  SU_TRYCATCH(
      new->mixer_hist = calloc(new->det_bank[0].hist_len, sizeof(SUCOMPLEX)),
      goto fail);

  return new;

fail:
  if (new != NULL)
    graves_det_destroy(new);

  return new;
}
