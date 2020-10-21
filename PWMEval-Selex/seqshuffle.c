/*

  Perform regional shuffling on a series of FASTA-formatted sequences.
  The shuffled sequences are written to std output in FASTA-format.
  
  Giovanna Ambrosini, EPFL/SV, giovanna.ambrosini@epfl.ch

  Copyright (c) 2015
  School of Life Sciences
  Ecole Polytechnique Federale de Lausanne
  and Swiss Institute of Bioinformatics
  EPFL SV ISREC UPNAE
  Station 15
  CH-1015 Lausanne, Switzerland.

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/
#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <errno.h>
#include <sys/types.h>
#include <unistd.h>
#include <ctype.h>
#include <limits.h>
#ifdef DEBUG
#include <mcheck.h>
#endif

#define BUF_SIZE 3072
#define NUCL  5
#define LMAX  100
#define HDR_MAX 132

typedef struct _options_t {
  int help;
  int debug;
  int nohdr;
  int seed_flag;
  unsigned int seed;
} options_t;

static options_t options;

static char nucleotide[] = {'A','C','G','T', 'N'};

typedef struct _seq_t {
  char *hdr;
  int *seq;
  int len;
} seq_t, *seq_p_t;

FILE *fasta_in;

int regLen = 0;

//Arrange the n elements of ARRAY in random order.
void 
shuffle(int *array, int n)
{   
  if (n > 1) {
    for (int i = n-1; i > 0; i--) {
      // Pick a random index from 0 to i
      int j = rand() % (i+1);
      // Swap array[i] with the element at random index
      int t = array[i];
      array[i] = array[j];
      array[j] = t;
    }
  }
}

static int
process_file(FILE *input, char *iFile)
{
  char buf[BUF_SIZE], *res;
  seq_t seq;
  int mLen;

  if (input == NULL) {
    FILE *f = fopen(iFile, "r");
    if (f == NULL) {
      fprintf(stderr, "Could not open file %s: %s(%d)\n",
  	    iFile, strerror(errno), errno);
      return -1;
    }
    input = f;
  }
  if (options.debug != 0)
    fprintf(stderr, "Processing file %s\n", iFile);
  while ((res = fgets(buf, BUF_SIZE, input)) != NULL
	 && buf[0] != '>')
    ;
  if (res == NULL || buf[0] != '>') {
    fprintf(stderr, "Could not find a sequence in file %s\n", iFile);
    if (input != stdin) {
      fclose(input);
    }
    return -1;
  }
  seq.hdr = malloc(HDR_MAX * sizeof(char));
  seq.seq = malloc(BUF_SIZE * sizeof(int));
  mLen = BUF_SIZE;
  while (res != NULL) {
    /* Get the header */
    if (buf[0] != '>') {
      fprintf(stderr, "Could not find a sequence header in file %s\n", iFile);
      if (input != stdin) {
        fclose(input);
      }
      return -1;
    }
    char *s = buf;
    s += 1;
    int i = 0;
    while (*s && !isspace(*s)) {
      if (i >= HDR_MAX) {
        fprintf(stderr, "Fasta Header too long \"%s\" in file %s\n", res, iFile);
        fclose(input);
        return -1;
      }
      seq.hdr[i++] = *s++;
    }
    if (i < HDR_MAX)
      seq.hdr[i] = 0;
    /* Gobble sequence  */ 
    seq.len = 0;
    while ((res = fgets(buf, BUF_SIZE, input)) != NULL && buf[0] != '>') {
      char c;
      int n;
      s = buf;
      while ((c = *s++) != 0) {
	if (isalpha(c)) {
	  c = (char) toupper(c);
	  switch (c) {
	  case 'A':
            n = 0;
            break;
	  case 'C':
            n = 1;
            break;
	  case 'G':
            n = 2;
            break;
	  case 'T':
            n = 3;
            break;
	  case 'N':
            n = 4;
	    break;
	  default:
            n = 4;
	    ;
	  }
	  if (seq.len >= mLen) {
	    mLen += BUF_SIZE;
	    seq.seq = realloc(seq.seq, (size_t)mLen * sizeof(int));
	  }
	  seq.seq[seq.len++] = n;
	}
      }
    }
    /* We now have the (not nul terminated) sequence.
       Process it. */
    if (seq.len != 0) {
      if (regLen == 0) { // shuffle entire sequence
        shuffle(seq.seq, seq.len);
        // Print out shuffled sequence
        // Print Sequence Header
        printf(">%s_shu\n", seq.hdr);
        int i = 0;
        for (i = 0; i < seq.len; i++) {
          printf("%c", nucleotide[seq.seq[i]]);
          if ( ((i + 1) % 60) == 0 ) {
            printf("\n");
          }
        }
        printf("\n");
      } else { // regional shuffling
        int i = 0;
        int cnt = 1;
        int k = 0;
        for (i = 0; i < (seq.len-regLen); i+=regLen) {
          if (options.debug != 0) {
            fprintf(stderr, "%d shuffling piece : i=%d reg Len=%d   ", cnt, i, regLen);
            for (k = i; k < i + regLen; k++) {
              fprintf(stderr, "%c", nucleotide[seq.seq[k]]);
            }
            fprintf(stderr, "\n");
            cnt++;
          }
          shuffle(&seq.seq[i], regLen); // shuffle each region separately
        }
        if ( i < (seq.len - 1) ) {
          int res = seq.len - i - 1;
          if (options.debug != 0) {
            fprintf(stderr, "Last piece: i=%d res=%d\n", i, res);
            fprintf(stderr, "%c\n", nucleotide[seq.seq[i]]);
          }
          shuffle(&seq.seq[i + 1], res);  // shuffle residual nucleotides
        }
        // Print out shuffled sequence
        // Print Sequence Header
        printf(">%s_shu\n", seq.hdr);
        i = 0;
        for (i = 0; i < seq.len; i++) {
          printf("%c", nucleotide[seq.seq[i]]);
          if ( ((i + 1) % 60) == 0 ) {
            printf("\n");
          }
        }
        printf("\n");
      }
    }
  }
  return 0;
}

int
main(int argc, char *argv[])
{
#ifdef DEBUG
  mcheck(NULL);
  mtrace();
#endif
  options.seed_flag = 0;
  while (1) {
    int c = getopt(argc, argv, "dhr:s:");
    if (c == -1)
      break;
    switch (c) {
    case 'd':
      options.debug = 1;
      break;
    case 'h':
      options.help = 1;
      break;
    case 'r':
      regLen = atoi(optarg);
      break;
    case 's':
      options.seed = atoi(optarg);
      options.seed_flag = 1;
      break;
    case '?':
      break;
    default:
      printf ("?? getopt returned character code 0%o ??\n", c);
    }
  }
  if (optind > argc || options.help == 1) {
    fprintf(stderr,
	    "Usage: %s [options] [<] [<fasta_file>|stdin]\n"
	    "      where options are:\n"
	    "        -d        Print debug information\n"
	    "        -h        Show this help text\n"
	    "        -r <len>  Shuffle sequence(s) in regions of <len>bp (by default <len>=0).\n"
	    "        -s <seed> Set the seed (integer) for the pseudo-random number generator algorithm.\n"
	    "                  By default, time(0) is used as seed.\n"
	    "\n\tPerform regional shuffling on a set of FASTA-formatted sequences-\n"
            "\tIf regional shuffling is not defined (option -r is not set), the entire\n"
            "\tsequence(s) is(are) shuffled.\n"
            "\tThe shuffled sequence(s) is(are) written to standard output.\n\n",
	    argv[0]);
    return 1;
  }

  if (argc > optind) {
      if(!strcmp(argv[optind],"-")) {
          fasta_in = stdin;
      } else {
          fasta_in = fopen(argv[optind], "r");
          if (fasta_in == NULL) {
              fprintf(stderr, "Unable to open '%s': %s(%d)\n",
                  argv[optind], strerror(errno), errno);
             exit(EXIT_FAILURE);
          }
          if (options.debug)
             fprintf(stderr, "Processing file %s\n", argv[optind]);
      }
  } else {
      fasta_in = stdin;
  }

  // Use a different seed value so that we don't get same 
  // result each time we run this program 
  if (options.seed_flag)
    srand (options.seed); 
  else 
    srand (time(NULL)); 

  if (options.debug != 0) {
    if (fasta_in != stdin) {
      fprintf(stderr, "Fasta Sequence File : %s\n", argv[optind]);
    } else {
      fprintf(stderr, "FASTA Sequence File from STDIN\n");
    }
    fprintf(stderr, "Regional Shuffling: %d\n", regLen);
  }
  
  if (process_file(fasta_in, argv[optind++]) != 0)
    return 1;

  return 0;
}
