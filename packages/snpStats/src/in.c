/* Read genotype data in long format

filename    Name of input text file (4 fields per row) 
tmpdir      Temporary directory for sorting (no trailing /)
threshold   Quality threshold for inclusion of data
nchip       Number of chips
chip_id     Array of chip identifier strings
nsnp        Number of snps
snp_id      Array of snp identifier strings
gtypes      Array of char's, length nchip*nsnp, to hold genotype data

*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <R.h>
#include <Rinternals.h>


#define MAX_ID 128
#define MAX_GT 16

void insnp(char *filename, char *tmpdir,  
	   int *nchip, char **chip_id, int *nsnps, char **snp_id,
	   char *codes[3], double *threshold, char *gtypes, 
	   int counts[2], int *iferror) {
  /* Sort with chips varying fastest */
  char sort_command[160];
  int i = 0, j = 0;
  sprintf(sort_command, 
	  "sort  -k 2,2 -k 1,1 -T \"%s\" -o \"%s\" \"%s\"",
          tmpdir, filename, filename);
  int error = system(sort_command);
  if (error) goto sort_error;
  FILE *infile = fopen(filename, "r");
  if (!infile) 
    goto open_error;

  int not_called = 0;
  int called = 0;
  char *code_aa = codes[0];
  char *code_ab = codes[1];
  char *code_bb = codes[2];
  char chip_in[MAX_ID], snp_in[MAX_ID], gt_in[MAX_GT];
  double thr_in;
  if (fscanf(infile, " %s %s %s %lf", chip_in, snp_in, gt_in, &thr_in)!=4)
    goto read_error;
  R_xlen_t ij=0;
  for (j=0; j<(*nsnps); j++) {
    char *snp_target = snp_id[j];
    int jcmp;
    while ((jcmp = strcmp(snp_in, snp_target))<0) {
      int scanned = fscanf(infile, " %s %s %s %lf", 
		       chip_in, snp_in, gt_in, &thr_in);
      if (scanned==EOF) goto normal;
      else if (scanned!=4)
	goto read_error;  
    }   
    for (i=0; i<(*nchip); i++, ij++) {
      char * chip_target = chip_id[i];
      int icmp; 
      if (!jcmp) {
	while ((icmp = strcmp(chip_in, chip_target))<0) { 
	  int scanned = fscanf(infile, " %s %s %s %lf", 
			       chip_in, snp_in, gt_in, &thr_in);
	  if (scanned == EOF) 
	    goto normal;
	  else if (scanned!=4)
	    goto read_error;
	}
	if (!icmp) {
	  /* Assign genotype */
	  if (!strcmp(code_aa, gt_in)) {
	    gtypes[ij] = (char) 1;
	    called ++;
	  }
	  else if  (!strcmp(code_ab, gt_in)) {
	    gtypes[ij] = (char) 2;
	    called ++;
	  }
	  else if  (!strcmp(code_bb, gt_in)) {
	    gtypes[ij] = (char) 3;
	    called ++;
	  }
	  else {
	    gtypes[ij] = (char) 0;
	    not_called ++;
	  }
	}
      }
      else {
	gtypes[ij] = (char) 0;
      }
    }
  }

  /* Normal return */

 normal: {
  int full = (*nsnps)*(*nchip);
  while (ij<full)
    gtypes[ij++] = (char)0;
  counts[0] = called;
  counts[1] = not_called;
  *iferror = 0;
  return;
  }
  

    /* Error conditions */ 
  
 sort_error: *iferror = 1; return;
 open_error: *iferror = 2; return;
 read_error: *iferror = 3; return;
 
}
