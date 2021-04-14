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

#ifndef FILENAME
#  define FILENAME __FILENAME__
#endif /* FILENAME */

#include <stdio.h>
#include <clistones.h>
#include <errno.h>
#include <sigutils/sigutils.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <getopt.h>

SUPRIVATE SUBOOL
clistones_register_chirp(
    clistones_t *self,
    struct clistones_chirp_summary *summary,
    struct timeval tv,
    const struct graves_chirp_info *chirp)
{
  FILE *fp = NULL;
  unsigned int i;
  char *path = NULL;
  SUFLOAT K, offset, snr, doppler;
  SUFLOAT cum_doppler = 0;
  SUFLOAT cum_snr = 0;
  SUFLOAT max_snr = 0;
  SUFLOAT val;
  SUCOMPLEX prev = 0;
  SUBOOL ok = SU_FALSE;

  SU_TRYCATCH(
      path = strbuild("%s/event_%06d.dat", self->directory, self->event_count),
      goto done);

  if ((fp = fopen(path, "wb")) == NULL) {
    SU_ERROR("Failed to open `%s' for writing: %s\n", path, strerror(errno));
    goto done;
  }

  /* Save data */
  SU_TRYCATCH(
      fprintf(fp, "EVENT_INDEX     =%15d", self->event_count) > 0,
      goto done);

  SU_TRYCATCH(
      fprintf(fp, "TIMESTAMP_SEC   =%15d", tv.tv_sec) > 0,
      goto done);

  SU_TRYCATCH(
      fprintf(fp, "TIMESTAMP_USEC  =%15d", tv.tv_usec) > 0,
      goto done);

  SU_TRYCATCH(
      fprintf(fp, "SAMPLE_RATE     =%15d", self->det_params.fs) > 0,
      goto done);

  SU_TRYCATCH(
      fprintf(fp, "CAPTURE_LEN     =%15d", chirp->length) > 0,
      goto done);

  SU_TRYCATCH(fprintf(fp, "DATA SECTION START              ") > 0, goto done);

  /* Save I/Q block */
  SU_TRYCATCH(
      fwrite(chirp->x, chirp->length * sizeof(SUCOMPLEX), 1, fp) == 1,
      goto done);

  /* Save SNR block */
  for (i = 0; i < chirp->length; ++i) {
    val = graves_det_q_to_snr(
            self->det_params.lpf2 / self->det_params.lpf1,
            chirp->q[i]);
    SU_TRYCATCH(fwrite(&val, sizeof(SUFLOAT), 1, fp) == 1, goto done);
  }

  /* Do some post processing on the chirp data */
  K = self->det_params.fs * SU_ADDSFX(.25) * SPEED_OF_LIGHT /
      (GRAVES_CENTER_FREQ * SU_ADDSFX(M_PI));

  /* Save Doppler block */
  max_snr = 0;
  for (i = 0; i < chirp->length; ++i) {
    offset = SU_C_ARG(chirp->x[i] * SU_C_CONJ(prev));
    prev = chirp->x[i];
    doppler = K * offset;

    snr = graves_det_q_to_snr(
            self->det_params.lpf2 / self->det_params.lpf1,
            chirp->q[i]);
    cum_snr += snr;
    if (snr > max_snr)
      max_snr = snr;

    cum_doppler += doppler * snr;
    SU_TRYCATCH(fwrite(&doppler, sizeof(SUFLOAT), 1, fp) == 1, goto done);
  }

  summary->index    = self->event_count;
  summary->tv       = tv;
  summary->duration = chirp->length / SU_ASFLOAT(self->det_params.fs);
  summary->mean_snr = cum_snr / chirp->length;
  summary->max_snr  = max_snr;
  summary->mean_vel = cum_doppler / cum_snr;

  summary->weak     = summary->max_snr < self->params.snr_threshold ||
      summary->duration < self->params.duration_threshold;

  if (summary->weak)
    (void) unlink(path);

  ok = SU_TRUE;

done:
  if (fp != NULL)
    fclose(fp);

  if (path != NULL)
    free(path);

  return ok;
}

void
clistones_cancel(clistones_t *self)
{
  self->cancelled = SU_TRUE;
}

SUPRIVATE SUBOOL
clistones_on_chirp(void *privdata, const struct graves_chirp_info *chirp)
{
  struct clistones_chirp_summary summary;
  clistones_t *self = (clistones_t *) privdata;
  SUBOOL ok = SU_FALSE;
  SUFLOAT snr, delta_t;
  unsigned int ticks, i;

  struct timeval now, prev, sub;
  struct tm *tm;

  gettimeofday(&now, NULL);

  SU_TRYCATCH(clistones_register_chirp(self, &summary, now, chirp), goto done);

  /* We ignore weak chirps */
  if (!summary.weak) {
    tm = gmtime(&now.tv_sec);

    printf(
        "[%04d/%02d/%02d - %02d:%02d:%02d U] ",
        tm->tm_year + 1900,
        tm->tm_mon  + 1,
        tm->tm_mday,
        tm->tm_hour,
        tm->tm_min,
        tm->tm_sec);

        snr = SU_POWER_DB(summary.mean_snr);

        ticks = snr < 1 ? 1 : floor(snr);

        printf(
            "STONE EVENT %07d %6.2f s (%+6.2f m/s) SNR: %+6.2f dB (max %+6.2f dB) [",
            self->event_count + 1,
            summary.duration,
            summary.mean_vel,
            snr,
            SU_POWER_DB(summary.max_snr));

        if (ticks >= 10)
          printf("\033[1;31m");
        else if (ticks >= 5)
          printf("\033[1;33m");
        else
          printf("\033[1;32m");

        if (ticks >= 16)
          ticks = 16;

        for (i = 0; i < ticks; ++i)
          putchar('|');

        printf("\033[0m");

        if (ticks == 16) {
          --ticks;
          putchar('+');
        }

        for (i = 0; i < 16 - ticks; ++i)
          putchar(' ');
        putchar(']');
        printf("\n");

    SU_TRYCATCH(
        fprintf(
            self->logfp,
            "%d,%d.%d,%e.10,%e.10,%e.10,%e.10,%e.10\n",
            summary.index,
            summary.tv.tv_sec,
            summary.tv.tv_usec,
            summary.duration,
            summary.mean_snr,
            summary.max_snr,
            summary.mean_vel) > 0,
        goto done);

    fflush(self->logfp);

    ++self->event_count;

    /* Show ZHR notice */
    if (self->params.cycle_len > 0) {
      if ((self->event_count % self->params.cycle_len) == 0) {
        if (self->event_count > 0) {
          timersub(&now, &self->first, &sub);

          delta_t = (sub.tv_sec + 1e-6 * sub.tv_usec);
          printf(
              "[%04d/%02d/%02d - %02d:%02d:%02d U] ",
              tm->tm_year + 1900,
              tm->tm_mon  + 1,
              tm->tm_mday,
              tm->tm_hour,
              tm->tm_min,
              tm->tm_sec);
          printf(
              "ZHR report update: %g events / hour\n",
              3600. * self->params.cycle_len / delta_t);
        }

        self->first = now;
      }
    }
  }

  ok = SU_TRUE;

done:
  return ok;
}

SUBOOL
clistones_loop(clistones_t *self)
{
  int err;
  unsigned int i;
  SUBOOL ok = SU_FALSE;

  while (!self->cancelled) {
    /* Read samples from soundcard */
    if ((err = snd_pcm_readi(self->pcm, self->buffer, CLISTONES_READ_SIZE))
        != CLISTONES_READ_SIZE) {
      SU_ERROR(
          "Error %d while capturing samples: %s\n",
          err,
          snd_strerror (err));
      goto done;
    }

    /* Forward them to meteorite detector */
    for (i = 0; i < CLISTONES_READ_SIZE; ++i)
      SU_TRYCATCH(
          graves_det_feed(self->detector, self->buffer[i] / 65535.),
          goto done);
  }

  ok = SU_TRUE;

done:
  self->cancelled = SU_TRUE;

  return ok;
}

SUPRIVATE snd_pcm_t *
clistones_open_audio(const struct clistones_params *params)
{
  int i, err;
  unsigned int rate = CLISTONES_SAMP_RATE;
  snd_pcm_t *capture_handle = NULL;
  snd_pcm_hw_params_t *hw_params;
  snd_pcm_format_t format = SND_PCM_FORMAT_S16_LE;
  SUBOOL ok = SU_FALSE;

  if ((err = snd_pcm_open(
      &capture_handle,
      params->device,
      SND_PCM_STREAM_CAPTURE,
      0)) < 0) {
    SU_ERROR(
        "Cannot open audio device `%s' (%s)\n",
        params->device,
        snd_strerror(err));
    goto done;
  }

  if ((err = snd_pcm_hw_params_malloc(&hw_params)) < 0) {
    SU_ERROR(
        "Cannot allocate hardware parameter structure (%s)\n",
        snd_strerror(err));
    goto done;
  }

  if ((err = snd_pcm_hw_params_any(capture_handle, hw_params)) < 0) {
    SU_ERROR(
        "Cannot initialize hardware parameter structure (%s)\n",
        snd_strerror(err));
    goto done;
  }

  if ((err = snd_pcm_hw_params_set_access(
      capture_handle,
      hw_params,
      SND_PCM_ACCESS_RW_INTERLEAVED)) < 0) {
    SU_ERROR("Cannot set access type (%s)\n", snd_strerror (err));
    goto done;
  }

  if ((err = snd_pcm_hw_params_set_format(
      capture_handle,
      hw_params,
      format)) < 0) {
    SU_ERROR("Cannot set sample format (%s)\n", snd_strerror (err));
    goto done;
  }

  if ((err = snd_pcm_hw_params_set_rate_near(
      capture_handle,
      hw_params,
      &rate,
      0)) < 0) {
    SU_ERROR("Cannot set sample rate (%s)\n", snd_strerror (err));
    goto done;
  }

  if (rate != CLISTONES_SAMP_RATE) {
    SU_ERROR(
        "Sample rate %d Hz not supported (offered %d instead)\n",
        CLISTONES_SAMP_RATE,
        rate);
    goto done;
  }

  if ((err = snd_pcm_hw_params_set_channels(
      capture_handle,
      hw_params,
      1)) < 0) {
    SU_ERROR("Cannot set channel count (%s)\n", snd_strerror(err));
    goto done;
  }

  if ((err = snd_pcm_hw_params(capture_handle, hw_params)) < 0) {
    SU_ERROR("Cannot set parameters (%s)\n", snd_strerror(err));
    goto done;
  }

  if ((err = snd_pcm_prepare (capture_handle)) < 0) {
    SU_ERROR(
        "Cannot prepare audio interface for use (%s)\n",
        snd_strerror(err));
    goto done;
  }

  ok = SU_TRUE;

done:
  if (hw_params != NULL)
    snd_pcm_hw_params_free(hw_params);

  if (!ok && capture_handle != NULL) {
    snd_pcm_close(capture_handle);
    capture_handle = NULL;
  }

  return capture_handle;
}

clistones_t *
clistones_new(const struct clistones_params *params)
{
  struct graves_det_params det_params = graves_det_params_INITIALIZER;
  char *default_directory = NULL;
  const char *directory;
  char *path = NULL;
  clistones_t *new = NULL;
  time_t t;
  struct tm *tm;

  SU_TRYCATCH(new = calloc(1, sizeof (clistones_t)), goto fail);
  new->params = *params;

  if (params->output_dir == NULL) {
    time(&t);
    tm = gmtime(&t);
    SU_TRYCATCH(
        new->directory = strbuild(
            "clistones_%04d%02d%02d_%02d%02d%02d",
            tm->tm_year + 1900,
            tm->tm_mon + 1,
            tm->tm_mday,
            tm->tm_hour,
            tm->tm_min,
            tm->tm_sec),
        goto fail);
  } else {
    SU_TRYCATCH(new->directory = strdup(params->output_dir), goto fail);
  }

  directory = new->directory;

  SU_TRYCATCH(
      new->buffer = malloc(sizeof(uint16_t) * CLISTONES_READ_SIZE),
      goto fail);

  det_params.fs = CLISTONES_SAMP_RATE;
  new->det_params = det_params;

  SU_TRYCATCH(
      new->detector = graves_det_new(
          &new->det_params,
          clistones_on_chirp,
          new),
      goto fail);

  SU_TRYCATCH(new->pcm = clistones_open_audio(params), goto fail);

  if (strcmp(directory, ".") != 0 && access(directory, F_OK) == -1) {
    if (mkdir(directory, 0755) == -1 && errno != EEXIST) {
      SU_ERROR(
          "Failed to create output directory `%s': %s\n",
          directory,
          strerror(errno));
      goto fail;
    }
  }

  SU_TRYCATCH(path = strbuild("%s/events.csv", directory), goto fail);

  if ((new->logfp = fopen(path, "w")) == NULL) {
    SU_ERROR(
        "Failed to create event log file `%s': %s\n",
        path,
        strerror(errno));
    goto fail;
  }

  gettimeofday(&new->first, NULL);

  return new;

fail:
  if (path != NULL)
    free(path);

  if (default_directory != NULL)
    free(default_directory);

  if (new != NULL)
    clistones_destroy(new);

  return NULL;
}

void
clistones_destroy(clistones_t *self)
{
  if (self->pcm != NULL)
    snd_pcm_close(self->pcm);

  if (self->logfp != NULL)
    fclose(self->logfp);

  if (self->detector != NULL)
    graves_det_destroy(self->detector);

  if (self->buffer != NULL)
    free(self->buffer);

  free(self);
}

void
help(const char *a0)
{
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "  %s [OPTIONS]\n\n", a0);
  fprintf(stderr, "OPTIONS:\n");
  fprintf(stderr, "  -d, --device=DEV  Sets ALSA capture device to DEV\n");
  fprintf(stderr, "  -o, --dir=DIR     Sets the output data directory to DIR\n");
  fprintf(stderr, "  -s, --snr=SNR_DB  Sets the SNR threshold for detection (dB)\n");
  fprintf(stderr, "  -t, --duration=T  Sets the duration threshold in seconds\n");
  fprintf(stderr, "  -Z, --zhr=EVENTS  Sets the ZHR report update interval\n\n");
  fprintf(stderr, "  -h, --help        This help\n");
}

static struct option long_options[] =
{
  {"device",   required_argument, 0, 'd'},
  {"dir",      required_argument, 0, 'o'},
  {"snr",      required_argument, 0, 's'},
  {"duration", required_argument, 0, 't'},
  {"zhr",      required_argument, 0, 'Z'},
  {"help",     no_argument, 0, 'h'},
  {0, 0, 0, 0}
};

int
main(int argc, char **argv, char **envp)
{
  clistones_t *clistones = NULL;
  struct clistones_params params = clistones_params_INITIALIZER;
  int ret = EXIT_FAILURE;
  int option_index = 0;
  int c;

  if (!su_lib_init()) {
    fprintf(stderr, "%s: failed to initialize library\n", argv[0]);
    goto done;
  }

  for (;;) {
    c = getopt_long(argc, argv, "d:o:s:t:Z:h", long_options, &option_index);

    if (c == -1)
      break;

    switch (c) {
      case 'd':
        params.device = optarg;
        break;

      case 'o':
        params.output_dir = optarg;
        break;

      case 's':
        if (sscanf(optarg, "%g", &params.snr_threshold) < 1) {
          fprintf(stderr, "%s: invalid SNR\n\n", argv[0]);
          help(argv[0]);
          goto done;
        }

        params.snr_threshold = SU_POWER_MAG(params.snr_threshold);
        break;

      case 't':
        if (sscanf(optarg, "%g", &params.duration_threshold) < 1) {
          fprintf(stderr, "%s: invalid duration\n\n", argv[0]);
          help(argv[0]);
          goto done;
        }
        break;

      case 'Z':
        if (sscanf(optarg, "%u", &params.cycle_len) < 1) {
          fprintf(stderr, "%s: invalid ZHR update interval\n\n", argv[0]);
          help(argv[0]);
          goto done;
        }
        break;

      case 'h':
        help(argv[0]);
        ret = EXIT_SUCCESS;
        goto done;

      case '?':
      case ':':
        help(argv[0]);
        goto done;

      default:
        abort ();
    }
  }

  if ((clistones = clistones_new(&params)) == NULL) {
    fprintf(stderr, "%s: failed to create clistones object\n", argv[0]);
    goto done;
  }

  printf(
      "Welcome to...\n"
      "   _____ _ _  _____ _                        \n"
      "  / ____| (_)/ ____| |                       \n"
      " | |    | |_| (___ | |_ ___  _ __   ___  ___ \n"
      " | |    | | |\\___ \\| __/ _ \\| '_ \\ / _ \\/ __|\n"
      " | |____| | |____) | || (_) | | | |  __/\\__ \\\n"
      "  \\_____|_|_|_____/ \\__\\___/|_| |_|\\___||___/\n"
      "                                             \n"
      "      The automatic meteor echo detector\n");
  printf("\n");
  printf("Brought to you with love and kindness by Gonzalo J. Carracedo\n\n");
  printf("  Listening samples from audio device \"%s\"\n", params.device);
  printf("  Data directory:  %s\n", clistones->directory);
  printf("  SNR threshold:   %g dB\n", SU_POWER_DB(params.snr_threshold));
  printf("  Min duration:    %g seconds\n", params.duration_threshold);

  if (params.cycle_len != 0)
    printf("  ZHR report update every %d events\n", params.cycle_len);
  else
    printf("  ZHR reports disabled\n");

  printf("\n");

  SU_TRYCATCH(clistones_loop(clistones), goto done);

  ret = EXIT_SUCCESS;

done:
  if (clistones != NULL)
    clistones_destroy(clistones);

  return ret;
}
