/*

  Copyright (C) 2021 Gonzalo José Carracedo Carballal

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

#ifndef _CLISTONES_CLISTONES_H
#define _CLISTONES_CLISTONES_H

#include <graves.h>
#include <alsa/asoundlib.h>
#include <stdint.h>

#define CLISTONES_SAMP_RATE 8000
#define CLISTONES_READ_SIZE  128

struct clistones_params {
  const char *output_dir;
  const char *device;
  SUFLOAT freq_offset;
  SUFLOAT snr_threshold;
  SUFLOAT duration_threshold;
  unsigned int cycle_len;
};

#define clistones_params_INITIALIZER    \
{                                       \
  NULL,      /* output_dir */           \
  "default", /* device */               \
  1000.,     /* freq_offset */          \
  1,         /* snr_threshold */        \
  0.25,      /* duration_threshold */   \
  10         /* cycle_len */            \
}

struct clistones_chirp_summary {
  unsigned int index;
  struct timeval tv;
  SUFLOAT duration;
  SUFLOAT mean_snr;
  SUFLOAT max_snr;
  SUFLOAT mean_vel;
  SUBOOL  weak;
};

struct clistones {
  struct clistones_params params;
  struct graves_det_params det_params;
  char *directory;
  graves_det_t *detector;
  snd_pcm_t *pcm;
  FILE *logfp;

  unsigned int event_count;

  uint16_t *buffer;
  SUBOOL cancelled;
  struct timeval first;
};

typedef struct clistones clistones_t;

SUINLINE const char *
clistones_data_directory(const clistones_t *self)
{
  return self->directory;
}

clistones_t *clistones_new(const struct clistones_params *params);
void clistones_cancel(clistones_t *self);
SUBOOL clistones_loop(clistones_t *self);
void clistones_destroy(clistones_t *self);

#endif /* _CLISTONES_CLISTONES_H */
