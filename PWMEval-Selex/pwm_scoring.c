/*

  Score a set of nucloetide sequences in FASTA format based
  on matches to a sequence motif represented by a position 
  weight matrix (PWM) or a base probability matrix (LPM)
  
  Giovanna Ambrosini, EPFL/SV, giovanna.ambrosini@epfl.ch

  Copyright (c) 2014
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
#include <math.h>
#include <errno.h>
#include <sys/types.h>
#include <unistd.h>
#include <ctype.h>
#include <getopt.h>
#include <assert.h>
#include <float.h>
#include <limits.h>
#ifdef DEBUG
#include <mcheck.h>
#endif

#define BUF_SIZE 3072
#define NUCL  5
#define LINE_SIZE 1024
#define MVAL_MAX 32

#define HDR_MAX 132
#define BEST_HIT_POS 256
/*#define MIN_SCORE -5000000 */
#define MIN_SCORE INT_MIN

typedef struct _options_t {
  int help;
  int debug;
  int norm;
  int pwm;
  int lpm;
  int seq_norm;
  int lib_norm;
  int nohdr;
  int bestscore;
  int forward;
} options_t;

static options_t options;

static char nucleotide[] = {'A','C','G','T', 'N'};
static double bg[] = {1.0,1.0,1.0,1.0,0.25};

typedef struct _seq_t {
  char *hdr;
  int *seq;
  int len;
} seq_t, *seq_p_t;

FILE *fasta_in;

int seqCnt;
double **lpm;                /* Letter Probability Matrix  */
int   **pwm;                 /* Position Weight Matrix     */
int matLen = 10;             /* Matrix Length              */

double pseudo_weight = 0.0;  /* Optional pseudo-weight for Letter Probability Matrix */ 

static int 
read_profile(char *iFile)
{
  FILE *f = fopen(iFile, "r");
  int l = 0;
  char *s, *res, *buf;
  size_t bLen = LINE_SIZE;
  int p_len = matLen;
  char mval[MVAL_MAX] = "";
  int i;

  if (f == NULL) {
    fprintf(stderr, "Could not open file %s: %s(%d)\n",
            iFile, strerror(errno), errno);
    return -1;
  }
  if (options.debug != 0)
    fprintf(stderr, "Processing file %s\n", iFile);

  if ((s = malloc(bLen * sizeof(char))) == NULL) {
    perror("process_sga: malloc");
    return(-1);
  }
  if (options.lpm) {
    /* Read Matrix file line by line */
    while ((res = fgets(s, (int) bLen, f)) != NULL) {
      size_t cLen = strlen(s);
  
      while (cLen + 1 == bLen && s[cLen - 1] != '\n') {
        bLen *= 2;
        if ((s = realloc(s, bLen)) == NULL) {
          perror("process_file: realloc");
          exit(1);
        }
        res = fgets(s + cLen, (int) (bLen - cLen), f);
        cLen = strlen(s);
      }
      if (s[cLen - 1] == '\n')
        s[cLen - 1] = 0;
  
      if (s[cLen - 2] == '\r')
        s[cLen - 2] = 0;
  
      buf = s;
      /* Get PWM fields */
      /* Get first character: if # or > skip line */
      if (*buf == '#' || *buf == '>')
        continue;
      /* Read First column value */
      while (isspace(*buf))
        buf++;
      i = 0;
      while (isdigit(*buf) || *buf == '-' || *buf == '.' || *buf == 'e' || *buf == 'E') {
        if (i >= MVAL_MAX) {
          fprintf(stderr, "Matrix value is too large \"%s\" \n", buf);
          exit(1);
       }
        mval[i++] = *buf++;
      }
      mval[i] = 0;
     if (strlen(mval) == 0) {
        fprintf(stderr, "Matrix value for colum 1 (row %d) is missing, please check the matrix format (it should be LPM)\n", l);
        return(-1);
      }
      lpm[0][l] = atof(mval);
      while (isspace(*buf))
        buf++;
      /* Read Second column value */
      i = 0;
      while (isdigit(*buf) || *buf == '-' || *buf == '.' || *buf == 'e' || *buf == 'E') {
        if (i >= MVAL_MAX) {
          fprintf(stderr, "Matrix value is too large \"%s\" \n", buf);
          exit(1);
       }
        mval[i++] = *buf++;
      }
      mval[i] = 0;
      if (strlen(mval) == 0) {
        fprintf(stderr, "Matrix value for colum 2 (row %d) is missing, please check the matrix format (it should be LPM)\n", l);
        return(-1);
      }
      lpm[1][l] = atof(mval);
      while (isspace(*buf))
        buf++;
      /* Read Third column value */
      i = 0;
      while (isdigit(*buf) || *buf == '-' || *buf == '.' || *buf == 'e' || *buf == 'E') {
        if (i >= MVAL_MAX) {
          fprintf(stderr, "Matrix value is too large \"%s\" \n", buf);
          exit(1);
       }
        mval[i++] = *buf++;
      }
      mval[i] = 0;
      if (strlen(mval) == 0) {
        fprintf(stderr, "Matrix value for colum 3 (row %d) is missing, please check the matrix format (it should be LPM)\n", l);
        return(-1);
      }
      lpm[2][l] = atof(mval);
      while (isspace(*buf))
        buf++;
      /* Read fourth column value */
      i = 0;
      while (isdigit(*buf) || *buf == '-' || *buf == '.' || *buf == 'e' || *buf == 'E') {
        if (i >= MVAL_MAX) {
          fprintf(stderr, "Matrix value is too large \"%s\" \n", buf);
          exit(1);
       }
        mval[i++] = *buf++;
      }
      mval[i] = 0;
      if (strlen(mval) == 0) {
        fprintf(stderr, "Matrix value for colum 4 (row %d) is missing, please check the matrix format (it should be LPM)\n", l);
        return(-1);
      }
      lpm[3][l] = atof(mval);
#ifdef DEBUG
      fprintf(stderr, "%3d   %f   %f   %f   %f\n", l, lpm[0][l], lpm[1][l], lpm[2][l], lpm[3][l]);
#endif
      if (l == p_len-1) {
       /* Reallocate Matrix columns */
        for ( int i = 0; i < NUCL; i++) {
          lpm[i] = realloc(lpm[i], p_len*2*sizeof(double));
          if (lpm[i] == NULL) {
            fprintf(stderr, "Out of memory\n");
            return 1;
          }
        }
        p_len *= 2;
      }
      l++;
    }
  } else {  /* Integer PWM  */
    /* Read Matrix file line by line */
    while ((res = fgets(s, (int) bLen, f)) != NULL) {
      size_t cLen = strlen(s);
  
      while (cLen + 1 == bLen && s[cLen - 1] != '\n') {
        bLen *= 2;
        if ((s = realloc(s, bLen)) == NULL) {
          perror("process_file: realloc");
          exit(1);
        }
        res = fgets(s + cLen, (int) (bLen - cLen), f);
        cLen = strlen(s);
      }
      if (s[cLen - 1] == '\n')
        s[cLen - 1] = 0;
  
      if (s[cLen - 2] == '\r')
        s[cLen - 2] = 0;
  
      buf = s;
      /* Get PWM fields */
      /* Get first character: if # or > skip line */
      if (*buf == '#' || *buf == '>')
        continue;
      /* Read First column value */
      while (isspace(*buf))
        buf++;
      i = 0;
      while (isdigit(*buf) || *buf == '-') {
        if (i >= MVAL_MAX) {
          fprintf(stderr, "Matrix value is too large \"%s\" \n", buf);
          exit(1);
       }
        mval[i++] = *buf++;
      }
      mval[i] = 0;
     if (strlen(mval) == 0) {
        fprintf(stderr, "Matrix value for colum 1 (row %d) is missing, please check the matrix format (it should be Integer)\n", l);
        return(-1);
      }
      pwm[0][l] = atoi(mval);
      while (isspace(*buf))
        buf++;
      /* Read Second column value */
      i = 0;
      while (isdigit(*buf) || *buf == '-') {
        if (i >= MVAL_MAX) {
          fprintf(stderr, "Matrix value is too large \"%s\" \n", buf);
          exit(1);
       }
        mval[i++] = *buf++;
      }
      mval[i] = 0;
      if (strlen(mval) == 0) {
        fprintf(stderr, "Matrix value for colum 2 (row %d) is missing, please check the matrix format (it should be Integer)\n", l);
        return(-1);
      }
      pwm[1][l] = atoi(mval);
      while (isspace(*buf))
        buf++;
      /* Read Third column value */
      i = 0;
      while (isdigit(*buf) || *buf == '-') {
        if (i >= MVAL_MAX) {
          fprintf(stderr, "Matrix value is too large \"%s\" \n", buf);
          exit(1);
       }
        mval[i++] = *buf++;
      }
      mval[i] = 0;
      if (strlen(mval) == 0) {
        fprintf(stderr, "Matrix value for colum 3 (row %d) is missing, please check the matrix format (it should be Integer)\n", l);
        return(-1);
      }
      pwm[2][l] = atoi(mval);
      while (isspace(*buf))
        buf++;
      /* Read fourth column value */
      i = 0;
      while (isdigit(*buf) || *buf == '-') {
        if (i >= MVAL_MAX) {
          fprintf(stderr, "Matrix value is too large \"%s\" \n", buf);
          exit(1);
       }
        mval[i++] = *buf++;
      }
      mval[i] = 0;
      if (strlen(mval) == 0) {
        fprintf(stderr, "Matrix value for colum 4 (row %d) is missing, please check the matrix format (it should be Integer)\n", l);
        return(-1);
      }
      pwm[3][l] = atoi(mval);
#ifdef DEBUG
      fprintf(stderr, "%3d   %7d   %7d   %7d   %7d\n", l, pwm[0][l], pwm[1][l], pwm[2][l], pwm[3][l]);
#endif
      if (l == p_len-1) {
       /* Reallocate Matrix columns */
        for ( int i = 0; i < NUCL; i++) {
          pwm[i] = realloc(pwm[i], p_len*2*sizeof(int));
          if (pwm[i] == NULL) {
            fprintf(stderr, "Out of memory\n");
            return 1;
          }
        }
        p_len *= 2;
      }
      l++;
    }
  }
#ifdef DEBUG
  fprintf(stderr, "PWM length: %d\n", l);
#endif
  fclose(f);
  return l;
}

static void
process_seq_lpm(seq_p_t seq, FILE *out)
{
  int i;
  int j;
  int nucl_cnt[] = {0, 0, 0, 0, 0};

  if (options.debug != 0) {
    fprintf(stderr, ">SEQ:  ");
    for (i = 0; i < seq->len; i++) {
      fprintf(stderr, "%d", seq->seq[i]);
    }
    fprintf(stderr, "\n");
  }

  if (options.seq_norm) {
    for (i = 0; i < seq->len; i++) {
      nucl_cnt[seq->seq[i]]++;
      /* fprintf(stderr, "nucl_cnt[%d] = %d\n", seq->seq[i], nucl_cnt[seq->seq[i]]); */
    }
    if (options.forward) {
      for (i = 0; i < NUCL-1; i++) {
        if (options.debug != 0)
          fprintf(stderr, "nucl_cnt[%d] = %d ; seq LEN = %d\n", i, nucl_cnt[i], seq->len); 
        bg[i] = (double)nucl_cnt[i]/(double)seq->len;
      }
    } else {
      double bcomp_at = (double) ((double)((double)((nucl_cnt[0]+nucl_cnt[3])/2)+(double)nucl_cnt[4]/4)/(double)seq->len);
      bg[0] = bcomp_at;
      bg[1] = (double) 0.5 - bcomp_at;
      bg[2] = (double) 0.5 - bcomp_at;
      bg[3] = bcomp_at;
    }
    if (options.debug != 0) {
      fprintf(stderr, "Background nucleotide frequencies: ");
      for (i = 0; i < NUCL; i++) {
        fprintf(stderr, "bg[%i] = %f ", i, bg[i]);
      }
      fprintf(stderr, "\n\n");
    }
  }
  if (options.bestscore) { // Compute the single best score
    double best_score = 0.0;
    //int best_pos = 0;
    char best_pos[BEST_HIT_POS] = "0";
    char strand = '+';
    for (i = 0; i <= seq->len-matLen; i++) {
      double prod = 1.0;
      double prod_rcomp = 1.0;
      for (j = 0; j < matLen; j++) {
        //printf ("i=%d , PWM[%d] [%d] = %.10f\n", i, seq->seq[i+j], j, lpm[seq->seq[i+j]][j]); 
        prod = prod * lpm[seq->seq[i+j]][j]/bg[seq->seq[i+j]];
        //printf ("PROD=%.10f\n", prod); 
        if (!options.forward) {
          int idx = 0;
          if (seq->seq[i+j] == 4) {
            idx = 4;
          } else {
            idx = 3-seq->seq[i+j];
          }
          //printf ("i=%d , RCPWM[%d] [%d] = %.10f\n", i, idx, matLen-j-1, lpm[idx][matLen-j-1]); 
          prod_rcomp = prod_rcomp * lpm[idx][matLen-j-1]/bg[idx];
          //printf ("RCProd=%.10f\n", prod_rcomp); 
        }
      }
      //printf ("prod: %g \t rev prod : %g\n", prod, prod_rcomp);
      double max = 0.0;
      if (options.forward)
        max = prod;
      else 
        max = prod > prod_rcomp ? prod : prod_rcomp;
      //printf ("MAX: %g (pos=%d, strand=%c, BEST: %g)\n", max, i, strand, best_score);
      if (max > best_score) {
        best_score = max;
        //best_pos = i;
        sprintf(best_pos, "%d", i);
        if (!options.forward) {
          if (max == prod) {
            strand = '+';
          } else {
            strand = '-';
            //best_pos = i + matLen;
            sprintf(best_pos, "%d", i + matLen);
          }
        }
      } else if (max == best_score && max != 0.0) {
        if (max == prod) {
            char res[16];
            sprintf(res, ",%d", i);
            strcat(best_pos, res);
        } else {
            char res[16];
            sprintf(res, ",%d", i + matLen);
            strcat(best_pos, res);
        }
      }
      //printf ("BEST SCORE (pos=%d, strand=%c): %f\n", i, strand, best_score);
    }
    if (options.debug != 0)
      fprintf(stderr, "%s\t%e\t%d\t%s\t%c\n", seq->hdr, best_score, seq->len, best_pos, strand);

    if (options.nohdr != 0)
      fprintf(out, "%g\t%d\t%s\t%c\n", best_score, seq->len, best_pos, strand);
    else
      fprintf(out, "%s\t%g\t%d\t%s\t%c\n", seq->hdr, best_score, seq->len, best_pos, strand);
  } else { // Compute sum of probabilities [both strands is the default]
    double sum = 0.0;
    for (i = 0; i <= seq->len-matLen; i++) {
      double prod = 1.0;
      double prod_rcomp = 1.0;
      for (j = 0; j < matLen; j++) {
        prod = prod * lpm[seq->seq[i+j]][j]/bg[seq->seq[i+j]];
        if (!options.forward) {
          int idx = 0;
          if (seq->seq[i+j] == 4) {
            idx = 4;
          } else {
            idx = 3-seq->seq[i+j];
          }
          prod_rcomp = prod_rcomp * lpm[idx][matLen-j-1]/bg[idx];
        }
      }
      if (options.forward) 
        sum = sum + prod;
      else
        sum = sum + prod + prod_rcomp;
    }
    if (options.debug != 0)
      fprintf(stderr, "%s\t%e\n", seq->hdr, sum);

    if (options.nohdr != 0)
      fprintf(out, "%g\n", sum);
    else
      fprintf(out, "%s\t%g\n", seq->hdr, sum);
  }
}

static void
process_seq_pwm(seq_p_t seq, FILE *out)
{
  int i;
  int j;
  char *tag_match = NULL;
  char *tag_match_pos = NULL;
  char *tag_match_rcomp = NULL;
  int score = 0;
  int rev_score = 0;
  int best_score = INT_MIN;
  int match_pos = 0;
  int strand = 0;

  if (options.debug != 0) {
    fprintf(stderr, "> ");
    for (i = 0; i < seq->len; i++) {
      fprintf(stderr, "%d", seq->seq[i]);
    }
    fprintf(stderr, "\n");
  }
  if (seq->len < matLen) {
    if (options.nohdr != 0)
      fprintf(out, "%d\t%d\t%s\t%d\t%c\n", 0, 0, "NOTAG", MIN_SCORE, '0');
    else
      fprintf(out, "%s\t%d\t%d\t%s\t%d\t%c\n", seq->hdr, 0, 0, "NOTAG", MIN_SCORE, '0');
    return;
  }
  tag_match = (char *)malloc((matLen + 1) * sizeof(char));
  tag_match_pos = (char *)malloc((matLen + 1) * sizeof(char));
  tag_match_rcomp = (char *)malloc((matLen + 1) * sizeof(char));
  int idx = 0;
  for (i = 0; i <= seq->len-matLen; i++) {
    score = 0;
    rev_score = 0;
    for (j = 0; j < matLen; j++) {
      /* Positive strand */
      score += pwm[seq->seq[i+j]][j];
      tag_match_pos[j] = nucleotide[seq->seq[i+j]];
      if (!options.forward) {
        /* Reverse Strand */
        idx = 0;
        if (seq->seq[i+j] == 4) {
           idx = 4;
        } else {
           idx = 3-seq->seq[i+j];
        }
        rev_score += pwm[idx][matLen-j-1];
        tag_match_rcomp[matLen-j-1] = nucleotide[idx];
      }
    }
    /* printf("i: %d j: %d\n", i, j); */
    tag_match_pos[j] = '\0';
    tag_match_rcomp[j] = '\0';
    /* printf("tag_matches: (pos) %s  (rev) %s\n", tag_match_pos, tag_match_rcomp); */
    /* printf("scores: pos %d  rev %d best %d\n", score, rev_score, best_score);  */
   /* Get max score between the positive and negative strands */
    int max = 0;
    int k = 0;
    if (options.forward) {
      max = score;
    } else { 
      int c = score - rev_score;
      k = (c >> 31) & 0x1;   /* check highest bit of integer c : 0 = positive, 1 = negative */
      max = score - k * c;
    }
    /* printf("max: %d\n", max);  */
    if (max > best_score) {
      best_score = max;
      match_pos = i;
      strand = k;
      if (k) {
        strcpy(tag_match, tag_match_rcomp);
      } else {
        strcpy(tag_match, tag_match_pos);
      }
      /* printf("New BEST SCORE: %d, pos %d, strand %d, match %s\n", best_score, match_pos, strand, tag_match); */
    }
  }
  char str;
  if (strand)
    str = '-';
  else
    str = '+';
  int match_end = match_pos + matLen;
  if (options.debug != 0)
    fprintf(stderr, "%s\t%d\t%d\t%s\t%d\t%c\n", seq->hdr, match_pos, match_end, tag_match, best_score, str);

  if (options.nohdr != 0)
    fprintf(out, "%d\t%d\t%s\t%d\t%c\n", match_pos, match_end, tag_match, best_score, str);
  else
    fprintf(out, "%s\t%d\t%d\t%s\t%d\t%c\n", seq->hdr, match_pos, match_end, tag_match, best_score, str);

  free(tag_match);
  free(tag_match_pos);
  free(tag_match_rcomp);
}

static int
process_file(FILE *input, char *iFile, FILE *out)
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
       Process it: once forward, and once in reverse. */
    if (seq.len != 0) {
      if (options.lpm)
        process_seq_lpm(&seq, out);
      else
        process_seq_pwm(&seq, out);
    }
  }
  free(seq.hdr);
  free(seq.seq);
  if (input != stdin) {
    fclose(input);
  }
  return 0;
}

char** str_split(char* a_str, const char a_delim)
{
    char** result = 0;
    size_t count = 0;
    char* tmp = a_str;
    char* last_comma = 0;
    char delim[2];
    delim[0] = a_delim;
    delim[1] = 0;

    /* Count how many elements will be extracted. */
    while (*tmp) {
        if (a_delim == *tmp) {
            count++;
            last_comma = tmp;
        }
        tmp++;
    }
    /* Add space for trailing token. */
    count += last_comma < (a_str + strlen(a_str) - 1);
    /* Add space for terminating null string so caller
       knows where the list of returned strings ends. */
    count++;
    result = malloc(sizeof(char*) *count);

    if (result) {
        size_t idx  = 0;
        char* token = strtok(a_str, delim);

        while (token) {
            assert(idx < count);
            *(result + idx++) = strdup(token);
            token = strtok(0, delim);
        }
        assert(idx == count - 1);
        *(result + idx) = 0;
    }
    return result;
}

int
main(int argc, char *argv[])
{
  char *matFile = NULL;
  char *bgProb = NULL;
  char** tokens;
  int i = 0;
  double bprob = 0.25; 
  options.lpm = 1;
  options.pwm = 0;

  static struct option long_options[] =
      {
          /* These options may or may not set a flag.
             We distinguish them by their indices. */
          {"best",    no_argument,       0, 'b'},
          {"debug",   no_argument,       0, 'd'},
          {"help",    no_argument,       0, 'h'},
          {"forward", no_argument,       0, 'f'},
          {"matrix",  required_argument, 0, 'm'},
          {"prob",    required_argument, 0, 'p'},
          {"seqnorm", no_argument,       0, 'q'},
          {"unorm",   no_argument,       0, 'u'},
          {"nohdr",   no_argument,       0, 'r'},
          {"pweight", required_argument, 0, 'w'},
          /* These options only set a flag. */
          {"lpm",     no_argument,       &options.lpm, 1},
          {"pwm",     no_argument,       &options.pwm, 1},
          {0, 0, 0, 0}
      };

  int option_index = 0;
  
#ifdef DEBUG
  mcheck(NULL);
  mtrace();
#endif
  while (1) {
    //int c = getopt(argc, argv, "dhl:m:p:qurw:");
    int c = getopt_long(argc, argv, "bdhfm:p:uqrw:", long_options, &option_index);
    if (c == -1)
      break;
    switch (c) {
    case 'b':
      options.bestscore = 1;
      break;
    case 'd':
      options.debug = 1;
      break;
    case 'h':
      options.help = 1;
      break;
    case 'f':
      options.forward = 1;
      break;
    case 'm':
      matFile = optarg;
      break;
    case 'p':
      bgProb = optarg;
      options.lib_norm = 1;
      break;
    case 'q':
      options.seq_norm = 1;
      break;
    case 'u':
      options.norm = 1;
      break;
    case 'r':
      options.nohdr = 1;
      break;
    case 'w':
      pseudo_weight = atof(optarg);
      break;
    case 0:
      /* If this option set a flag, do nothing else now. */
      if (long_options[option_index].flag != 0)
        break;
    case '?':
      break;
    default:
      printf ("?? getopt returned character code 0%o ??\n", c);
    }
  }
  if (optind > argc || matFile == NULL) {
    fprintf(stderr,
	    "Usage: %s [options] -m <matrix_file> [<] <fasta_file>\n"
	    "   where options are:\n"
	    "     -b[--best]             Compute best single match scores\n"
	    "     -d[--debug]            Produce debugging output\n"
	    "     -h[--help]             Show this stuff\n"
	    "     -f[--forward]          Scan sequences in forward direction [def=bidirectional]\n"
	    "     -u[--unorm]            Normalize pwm scores by a uniform background base composition (Default=0.25)\n"
	    "     -p[--prob] <bg freq>   Normalize pwm scores by library-dependent nucleotide frequencies <bg freq>: 0.29,0.21,0.21,0.29\n"
	    "                            Note that nucleotide frequencies (<bg freq>) MUST BE comma-separated\n"
	    "     -q[--seqnorm]          Normalize pwm scores by sequence-based nucleotide composition\n"
	    "     -r[--nohdr]            Output raw scores (with no FASTA header)\n"
	    "     --lpm                  Input matrix is a letter probability matrix (LPM) [Default]\n"
	    "     --pwm                  Input matrix is a position weight matrix (PWM)\n"
	    "     -w[--pweight]          Set a pseudo-weight to re-normalize the frequencies of the letter-probability matrix (LPM)\n"
	    "                            Recommended value is 0.0001 [Default=0.0]\n"
	    "\n   Score a set of nucleotide sequences in FASTA format (<fasta_file>), based on matches to a sequence motif\n"
            "   represented by an INTEGER position weight matrix [--pwm] or a base probability matrix [--lpm] (<matrix_file>).\n"
            "   Note that the background normalization options (-u, -p, -q) are only valid for base probability matrices.\n"
            "   For integer PWMs, only the best single match scores are reported, along with the position, strand, and sequence match.\n\n",
	    argv[0]);
    return 1;
  }
  if (options.pwm)
    options.lpm = 0;
  if (options.lpm) {
    /* Allocate space for profile (LPM) */
    lpm = (double **)calloc(NUCL, sizeof(double *)); /* Allocate rows */
    if (lpm == NULL) {
      fprintf(stderr, "Could not allocate matrix array: %s(%d)\n",
        strerror(errno), errno);
      return 1;
    }
    for (i = 0; i < NUCL; i++) {
      lpm[i] = calloc((size_t)matLen, sizeof(double));  /* Allocate columns for ACGT + N */
      if (lpm[i] == NULL) {
        fprintf(stderr, "Out of memory\n");
        return 1;
      }
    }
  } else {
    options.seq_norm = 0;
    options.norm = 0;
    options.lib_norm = 0;
    pwm = (int **)calloc(NUCL, sizeof(int *)); /* Allocate rows */
    if (pwm == NULL) {
      fprintf(stderr, "Could not allocate matrix array: %s(%d)\n",
          strerror(errno), errno);
      return 1;
    }
    for (i = 0; i < NUCL; i++) {
      pwm[i] = calloc((size_t)matLen, sizeof(int));  /* Allocate columns for ACGT + N */
      if (pwm[i] == NULL) {
        fprintf(stderr, "Out of memory\n");
        return 1;
      }
    }
  }
  /* Read Matrix from file */
  if ((matLen = read_profile(matFile)) <= 0)
    return 1;
  /* Fill 5th pwm column for the N nucleotide (0.25) */
  if (options.lpm) {
    for (i = 0; i < matLen; i++) {
        lpm[4][i] = bprob;
    }
    if (pseudo_weight != 0.0) {
      /* Re-normalize the matrix by adding a pseudo-weight to the bease frequencies */
      for (int j = 0; j < matLen; j++) {
        double sum = 0.0;
        for (int i = 0; i < NUCL-1; i++)
          sum += lpm[i][j] + pseudo_weight;
        for (int i = 0; i < NUCL-1; i++)
          lpm[i][j] = (lpm[i][j] + pseudo_weight)/sum;
      }
    }
  } else {
    for (i = 0; i < matLen; i++) {
        pwm[4][i] = INT_MIN;
    } 
  }
  /* Treat background nucleotide frequencies */
  if (options.norm) {
    for (int j = 0; j < matLen; j++) {
      for (i = 0; i < NUCL-1; i++) {
        bg[i] = bprob;
      }
    }
  } else if (options.lib_norm) {
    tokens = str_split(bgProb, ',');
    if (tokens) {
      int i;
      for (i = 0; *(tokens + i); i++) {
        bg[i] =  atof(*(tokens + i));
        free(*(tokens + i));
      }
      if (i != 4) {
        fprintf(stderr, "Number of TOKENS: %d\n", i);
        fprintf(stderr, "Please, specify correct library-dependent nucleotide frequencies <bg freq>: they MUST BE comma-separated!\n");
        exit(1);
      }
      free(tokens);
    }
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

  if (options.debug != 0) {
    if (fasta_in != stdin) {
      fprintf(stderr, "Fasta File : %s\n", argv[optind]);
    } else {
      fprintf(stderr, "Sequence File from STDIN\n");
    }
    fprintf(stderr, "Motif length: %d\n", matLen);
    fprintf(stderr, "Weight Matrix: \n\n");
    if (options.lpm) {
      double *p;
      for (int i = 0; i < NUCL; i++) {
        fprintf(stderr, "%c ", nucleotide[i]);
        fprintf(stderr, "[");
        /* processing the rows of the PWM */
        for (p = lpm[i]; p < lpm[i] + matLen; p++)
          fprintf(stderr, " %f ", (*p));
        fprintf(stderr, "]\n");
      }
    } else {
      int *p;
      for (int i = 0; i < NUCL; i++) {
        fprintf(stderr, "%c ", nucleotide[i]);
        fprintf(stderr, "[");
        /* processing the rows of the PWM */
        for (p = pwm[i]; p < pwm[i] + matLen; p++)
          fprintf(stderr, " %d ", (*p));
        fprintf(stderr, "]\n");
      }
    }
    fprintf(stderr, "\n");
    fprintf(stderr, "Weight Matrix: vertical representation (columns represent the four nucleotides ACGT)\n\n");
    int j = 0;
    for (j = 0; j < matLen; j++) {
      for ( int i = 0; i < NUCL-1; i++) {
        if (options.lpm) {
          double pval = lpm[i][j];
          fprintf(stderr, " %f ", pval);
        } else {
          int pval = pwm[i][j];
          fprintf(stderr, " %d ", pval);
        }
      }
      fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");

    if (options.lib_norm) {
      fprintf(stderr, "Background nucleotide frequencies:[%s]\n", bgProb);
      for (i = 0; i < NUCL; i++) {
        fprintf(stderr, "bg[%i] = %.2f ", i, bg[i]);
      }
      fprintf(stderr, "\n\n");
    }
    if (options.seq_norm) {
      fprintf(stderr, "Sequence-based nucleotide composition\n");
    }
    if (options.bestscore) {
      fprintf(stderr, "Compute best match scores instead of sum of probabilities\n");
    }
    if (options.forward) {
      fprintf(stderr, "Scanning sequences in forward direction only\n");
    }
    fprintf(stderr, "\n");
  }
  
  if (process_file(fasta_in, argv[optind++], stdout) != 0)
    return 1;
  
  if (options.lpm) {
    for (i = 0; i < NUCL; i++)
      free(lpm[i]);
    free(lpm);
  } else {
    for (i = 0; i < NUCL; i++)
      free(pwm[i]);
    free(pwm);
  }

  return 0;
}
