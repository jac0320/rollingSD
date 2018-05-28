//
//  arma.h
//  rollingSP
//

#ifndef arma_h
#define arma_h

#include "rollingSP.h"

/* TODO: Description about this structure... */
typedef struct {
    int     cnt;
    intvec  seq;
    intvec	row;        /* row number array in the original problem; -1 indicates objective function */
    intvec	col; 		/* column number array in the original problem; -1 indicates right-hand side */
} varBlkType;

/* TODO: Description about this structure... */
typedef struct {
    int     total;      /* The total length of total data stream */
    int     diff;       /* Indicate if this is a differenced data stream, if 0 not differenced */
    int     start;      /* Zero time starting point at obs vector */
    int     past;       /* Length of the pre-current data */
    int     length;     /* Scope of the total data stream */
    int     open;       /* Revealed data until this index, stream[open] is known data */
    vector  *obs;       /* Indexed by [time][variate], fetch using open+i */
}streamType;

/* TODO: Description about this sturcture... */
typedef struct {
    int     cnt;
    vector  trend;
    vector  seasonal;
}decomposeType;

/* This structrue store ths arma model information. The information is obtained externally from a .arma file. More description will be here... */
typedef struct {
    int             transType;      /* 1->Univriate; 2->Multivariate */
    int             subh;           /* 0->Hourly Approximation; 1->Subhourly Approximation */
    int             window;         /* Frequency of a day's wind */
    int             p;              /* Order of AR */
    int             d;              /* Itegrated Difference: Universal Diff */
    int             q;              /* Order of MA */
    int             var;            /* Number of variates */
    double          lb;
    vector          mean;           /* Noise Process Mean = 0, indexed by Variate */
    vector          stdev;          /* Noise Process Standard Deviation = 1, indexed by Variate*/
    sparseMatrix    **Phi;          /* AR coefficient matrix */
    sparseMatrix    **Theta;        /* MA coefficient matrix */
    varBlkType      **variate;      /* Coupling Section between time series and optimization model */
    streamType      *history;        /* Historical data stream */
    streamType      *noise;         /* Not sure if I wanted this */
    vector          *histARObs;      /* [t][v] Stores the time series deterministic parts for AR section */
    vector          *histMAObs;      /* [t][v] Stores the time series deterministic parts for MA section */
    decomposeType   **decomp;       /* If decomposition method is used for pre-process this is the structure that stores all data */
}tsType;

/* This is a newer, cleaner structure of the arma model. */
typedef struct {
    int             p;
    int             q;
    int             var;
    int             subh;
    varBlkType      **variate;
    vector          mean;
    vector          stdev;
    sparseMatrix    **Phi;
    sparseMatrix    **Theta;
    vector          *history;
    vector          *noise;
    vector          *trend;
}armaType;

/* Subroutines in arma.c */
tsType *readARIMA(string inputDir, string probName, oneProblem *orig, BOOL *arimaFlag, int *stat1);
int readARIMACore(FILE *fptr, tsType *arima, int *stat1);
int readARIMACoef(FILE *fptr, tsType *arima, int *stat1);
int readARIMAVarBlock(FILE *fptr, tsType *arima, oneProblem *orig, int *stat1);
int readARIMAStream(FILE *fptr, tsType *arima, int *stat1);
int readARIMAAux(FILE *fptr, tsType *arima, int *stat1);
int simulateNoise (tsType *arima, vector *noise, long long *seed);
void reverseVector(vector observ, tsType *arima, int t);
double reverse(double tsResult, int varID, int minPos, int hourPos, tsType *arima);
void freeTsType(tsType *arima);
void freeVarBlkType(varBlkType *variate);
void freeStreamType(streamType *stream, int var);
void freeDecomposeType(decomposeType *decompose);

/* Subroutine in rvgen.c */
void generateBlocks(stocType *stoc, vector observ, int groupID, long long *seed);
void generateIndep(stocType *stoc, vector observ, int groupID, long long *seed);
int normal(vector mu, vector stdev, int numOmega, vector observ, long long *seed);
float scalit(float lower, float upper, long long *RUN_SEED);
float randUniform(long long *SEED);
int randInteger(long long *SEED, int iMax);

#endif /* arma_h */
