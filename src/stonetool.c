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
  <http://www.gnu.org/licenses/>

*/

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <errno.h>
#include <stdint.h>
#include <fcntl.h>
#include <stdlib.h>
#include <getopt.h>
#include <sigutils/sigutils.h>
#include <string.h>

struct stonefile {
  int index;
  struct timeval tv;
  unsigned long samp_rate;
  unsigned int length;

  SUCOMPLEX *iq;
  SUFLOAT *snr;
  SUFLOAT *doppler;
};

typedef struct stonefile stonefile_t;

SUPRIVATE void
stonefile_destroy(stonefile_t *self)
{
  if (self->iq != NULL)
    free(self->iq);

  if (self->snr != NULL)
    free(self->snr);

  if (self->doppler != NULL)
    free(self->doppler);

  free(self);
}

SUPRIVATE SUBOOL
stonefile_parse_keys(stonefile_t *self, uint8_t *bytes, size_t size)
{
  char buffer[33];
  char *space;
  char *value;
  unsigned long longval;
  unsigned int p = 0;
  SUBOOL have_data = SU_FALSE;
  SUBOOL ok = SU_FALSE;
  buffer[32] = '\0';

  while ((p + 32) < size) {
    memcpy(buffer, bytes + p, 32);
    p += 32;

    if (memcmp(buffer, "DATA SECTION START              ", 32) == 0) {
      have_data = SU_TRUE;
      break;
    }

    if ((space = strchr(buffer, ' ')) == NULL) {
      if ((space = strchr(buffer, '=')) == NULL) {
        SU_ERROR("Invalid metadata entry at offset %p\n", p - 32);
        goto done;
      }
    }
    *space = '\0';
    value = buffer + 17;

    if (strcmp(buffer, "EVENT_INDEX") == 0) {
      if (sscanf(value, "%d", &self->index) < 1) {
        SU_ERROR("Invalid %s (%s)\n", buffer, value);
        goto done;
      }
    } else if (strcmp(buffer, "SAMPLE_RATE") == 0) {
      if (sscanf(value, "%d", &self->samp_rate) < 1) {
        SU_ERROR("Invalid %s\n", buffer);
        goto done;
      }

      /* I am SO sorry */
      if (p < size && bytes[p] == 'u')
        ++p;
    } else if (strcmp(buffer, "TIMESTAMP_SEC") == 0) {
      if (sscanf(value, "%lu", &longval) < 1) {
        SU_ERROR("Invalid %s\n", buffer);
        goto done;
      }

      self->tv.tv_sec = longval;
    } else if (strcmp(buffer, "TIMESTAMP_USEC") == 0) {
      if (sscanf(value, "%lu", &longval) < 1) {
        SU_ERROR("Invalid %s\n", buffer);
        goto done;
      }

      self->tv.tv_usec = longval;
    } else if (strcmp(buffer, "CAPTURE_LEN") == 0) {
      if (sscanf(value, "%u", &self->length) < 1) {
        SU_ERROR("Invalid %s\n", buffer);
        goto done;
      }
    }
  }

  if (!have_data) {
    SU_ERROR("File does not have a DATA section\n");
    goto done;
  }

  if (self->length > 0) {
    SU_TRYCATCH(
        self->iq = malloc(self->length * sizeof(SUCOMPLEX)),
        goto done);

    SU_TRYCATCH(
        self->snr = malloc(self->length * sizeof(SUFLOAT)),
        goto done);

    SU_TRYCATCH(
        self->doppler = malloc(self->length * sizeof(SUFLOAT)),
        goto done);

    memcpy(self->iq, bytes + p, self->length * sizeof(SUCOMPLEX));
    memcpy(self->snr, bytes + p, self->length * sizeof(SUFLOAT));
    memcpy(self->doppler, bytes + p, self->length * sizeof(SUFLOAT));
  }

  ok = SU_TRUE;

done:
  return ok;
}

SUPRIVATE void
stonefile_info(const stonefile_t *self, FILE *fp)
{
  char date[32];
  date[31] = '\0';

  fprintf(fp, "Event number: %d\n", self->index);
  fprintf(fp, "Sample rate:  %d\n", self->samp_rate);

  if (self->tv.tv_sec > 0) {
    strncpy(date, ctime(&self->tv.tv_sec), 31);
    date[strlen(date) - 1] = '\0';

    fprintf(fp, "Timestamp:    %s (+%d usec)\n", date, self->tv.tv_usec);
  }

  if (self->samp_rate > 0)
    fprintf(
        fp,
        "Duration:     %g s\n",
        self->length / (SUFLOAT) self->samp_rate);
  else
    fprintf(fp, "Duration:     %g samples\n", self->length);
}

SUPRIVATE SUBOOL
stonefile_helper_dump_float_array(
    const char *path,
    const SUFLOAT *data,
    size_t size)
{
  FILE *fp = NULL;
  SUBOOL ok = SU_FALSE;

  if ((fp = fopen(path, "wb")) == NULL) {
    SU_ERROR("Cannot open `%s' for writing: %s\n", path, strerror(errno));
    goto done;
  }

  if (fwrite(data, sizeof(SUFLOAT), size, fp) < size) {
    SU_ERROR("Write samples to `%s' failed: %s\n", path, strerror(errno));
    goto done;
  }

  ok = SU_TRUE;

done:
  if (fp != NULL)
    fclose(fp);

  return ok;
}

SUPRIVATE stonefile_t *
stonefile_new(const char *path)
{
  stonefile_t *new = NULL;
  void *map = (void *) -1;
  struct stat sbuf;
  int fd = -1;
  SUBOOL ok = SU_FALSE;

  if (stat(path, &sbuf) == -1) {
    SU_ERROR("Cannot stat `%s': %s\n", path, strerror(errno));
    goto done;
  }

  SU_TRYCATCH(new = calloc(1, sizeof(stonefile_t)), goto done);

  if ((fd = open(path, O_RDONLY)) == -1) {
    SU_ERROR("Cannot open `%s' for reading: %s\n", path, strerror(errno));
    goto done;
  }

  if ((map = mmap(
      NULL,
      sbuf.st_size,
      PROT_READ,
      MAP_PRIVATE,
      fd,
      0)) == (void *) -1) {
    SU_ERROR("Cannot mmap `%s': %s\n", path, strerror(errno));
    goto done;
  }

  SU_TRYCATCH(stonefile_parse_keys(new, map, sbuf.st_size), goto done);

  ok = SU_TRUE;

done:
  if (fd != -1)
    close(fd);

  if (map != (void *) -1)
    munmap(map, sbuf.st_size);

  if (!ok && new != NULL) {
    stonefile_destroy(new);
    new = NULL;
  }

  return new;
}

void
help(const char *a0)
{
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "  %s [OPTIONS] FILE [OUTPUT]\n\n", a0);
  fprintf(stderr, "OPTIONS:\n");
  fprintf(stderr, "  -d, --dump=SECT   Dumps section SECT to a file\n");
  fprintf(stderr, "  -h, --help        This help\n");
}

static struct option long_options[] =
{
  {"dump",   required_argument, 0, 'd'},
  {"help",   no_argument,       0, 'h'},
  {0, 0, 0, 0}
};

int
main(int argc, char **argv, char **envp)
{
  int ret = EXIT_FAILURE;
  int c;
  int option_index = 0;
  stonefile_t *file = NULL;
  const char *dump = NULL;

  if (!su_lib_init()) {
    fprintf(stderr, "%s: failed to initialize library\n", argv[0]);
    goto done;
  }

  for (;;) {
    c = getopt_long(argc, argv, "d:h", long_options, &option_index);

    if (c == -1)
      break;

    switch (c) {
      case 'd':
        dump = optarg;
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

  if (dump == NULL) {
    if (argc - optind != 1) {
      fprintf(stderr, "%s: expected one file argument\n", argv[0]);
      goto done;
    }

    SU_TRYCATCH(file = stonefile_new(argv[optind]), goto done);

    stonefile_info(file, stdout);
  } else {
    if (argc - optind != 2) {
      fprintf(stderr, "%s: expected two file arguments\n", argv[0]);
      goto done;
    }

    SU_TRYCATCH(file = stonefile_new(argv[optind]), goto done);

    stonefile_info(file, stderr);
    if (strcasecmp(dump, "doppler") == 0) {
      SU_TRYCATCH(
          stonefile_helper_dump_float_array(
              argv[optind + 1],
              file->doppler,
              file->length),
          goto done);

    } else if (strcasecmp(dump, "snr") == 0) {
      SU_TRYCATCH(
          stonefile_helper_dump_float_array(
              argv[optind + 1],
              file->snr,
              file->length),
          goto done);

    } else if (strcasecmp(dump, "iq") == 0) {
      SU_TRYCATCH(
          stonefile_helper_dump_float_array(
              argv[optind + 1],
              (SUFLOAT *) file->iq,
              2 * file->length),
          goto done);

    }
  }

  ret = EXIT_SUCCESS;

done:
  if (file != NULL)
    stonefile_destroy(file);

  return ret;
}
