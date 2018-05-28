//
//  gamma.h
//  rollingSP
//
//  Created by Site Wang on 4/10/16.
//  Copyright © 2016 Site Wang. All rights reserved.
//

#ifndef gamma_h
#define gamma_h

#include "rollingSP.h"

typedef struct {
    int 	cnt;
    vector 	*vals;
}dynLambdaType;

/*
 TODO: [SITE] Write Description
 */
typedef struct{
    vector          **aData;        /* holds all information about a vector including first stage and second stage */
    int             aCnt;           /* holds prob[1]->aBar->cnt */
    int             BCnt;           /* holds prob[1]->Bbar[1]->cnt */
    vector          aDelta;         /* holds the difference of aBar from last time period to current time period */
    vector          BDelta;         /* holds the difference of dynamic decisions in Bbar from last time period to now */
    dynLambdaType 	*aDynLambda;	/* holds dual solutions corresponding to rows effected by rhs dynamics */
    dynLambdaType   *BDynLambda;    /* holds dual solutions corresponding to rows effected by dynamic decisions */
}exoType;


typedef struct{
    sparseVector    *bARMA;
    sparseVector    *changebARMA;
    int             pibCnt;
    vector          pibARMA;
    intvec          pibARMA_ck;
    vector          optX;
    vector          optXChange;
    vector          optY;
    vector          optYChange;
} endoType;

/*
 In the rolling horizon algorithm, some random variables may be produced using a time series model (ARIMA). With the ARIMA model,
 it is necessary to store a few deterministic terms separately for the ease computing. Focus on the AR part of ARIMA(p,d,q) model,
 a sliding window with length p is moving from hour t-p to t in step size n. The overlapping section between the sliding window and
 historical data (may be differenced data) is deterministic and can be stored before hand. (Note that it does not depend on 
 observations of omega).
 Phi_1 * w_n-1 + Phi_2 * w_n-2 + ... + Phi_p * w_n-p        <-This line can also be utilized for rolling
 Phi_2 * w_n-1 + Phi_3 * w_n-3 + ... + Phi_p * w_n-p+1
 ...
 where w is used to denote historical observation (deterministic)
 This will save some time (depends on amount of sub-hourly intervals and order p/q) when a new omega is generated. The stored value is valid
 for time period t only. As time evolves from t->t+, gammaType is updated with the latest observations.
 This structure can also be utilized by ARIMA model for new genrate.
 */
typedef struct {
    endoType    *endo;
    exoType     *exo;
} gammaType;

/* Subroutines in gamma.c */
gammaType *newGamma(int numPeriod, int masterCols, int subprobCols, int bLength, int lambdaCnt);
void freeExdoType(exoType *exo);
void freeEndoType(endoType *endo);
void freeGammaType(gammaType *gamma);

#endif /* gamma_h */
