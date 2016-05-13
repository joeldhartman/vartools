/*
   This file is part of the VARTOOLS implementation of David Palmer's fastchi2 
   period search algorithm. The reference for this algorithm 
       Palmer 2009, ApJ, 695, 496.
   This file provides header's for the fastchi2.c file.
*/

#define FASTCHI2_DEFAULT_FREQMIN 0.0
#define FASTCHI2_DEFAULT_DETRENDORDER 0
#define FASTCHI2_DEFAULT_OVERSAMPLE 1.0
#define FASTCHI2_DEFAULT_CHIMARGIN 0.0
#define FASTCHI2_DEFAULT_NPEAK 1


typedef struct {
  int *Nharm;
  double *freqmax;
  double *freqmin;
  int *detrendorder;
  double *t0;
  double *timespan;
  double *oversample;
  double *chimargin;
  int Npeak;
  int dorefitpeak;
  char *peroutdir;
  char *modeloutdir;
  int outputper;
  int outputmodel;
  int outputper_useformat;
  char *outputper_format;
  int outputmodel_useformat;
  char *outputmodel_format;
  int fixfreqmin;
  int fixt0;
  int fixtimespan;
  double **bestfreq;
  float **chireduction;
  int saveoutputmodel;
  double **outputmodeldata;
} _Fastchi2;
