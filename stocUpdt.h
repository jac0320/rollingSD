//
//  stocUpdt.h
//  rollingSP
//
//  Created by Site Wang on 4/24/16.
//  Copyright © 2016 Site Wang. All rights reserved.
//


#ifndef stocUpdt_h
#define stocUpdt_h

#include "rollingSP.h"
#include "gamma.h"

/* subroutines in stocUpdt.c */
BOOL stochasticUpdates(probType *p, numType *num, vector x, coordType *coord, sparseVector *bBar, sparseMatrix *Cbar, sparseVector *aBar, sparseMatrix **Bbar,
                       lambdaType *lambda,gammaType *gamma, sigmaType *sigma, deltaType *delta, deltaType *auxDelta, deltaType *trimDelta, omegaType *omega,
                       omegaType *auxOmega, omegaType *trimOmega, BOOL newOmegaFlag, int omegaIdx, int maxIter, int iter, vector pi, double mubBar);
double calcObjGap(probType **p,numType *num, coordType *coord, sigmaType *sigma, deltaType *delta,
                  vector Xvect, int obs);
int calcOmega(vector observ, vector auxObserv, vector trimObs, int begin, int end, omegaType *omega, omegaType *auxOmega, omegaType *trimOmega, BOOL *newOmegaFlag);
void calcDeltaCol(numType *num, coordType *coord, lambdaType *lambda, vector observ, vector auxObs, vector trimObs, int omegaIdx, deltaType *delta, deltaType *auxDelta, deltaType *trimDelta);
int calcDeltaRow(int maxIter, numType *num, coordType *coord, omegaType *omega, lambdaType *lambda, omegaType *auxOmega, omegaType *trimOmega, int lambdaIdx, deltaType *delta, deltaType *auxDelta, deltaType *trimDelta);
int calcLambda(numType *num, coordType *coord, vector Pi, lambdaType *lambda, sigmaType *sigma, BOOL *newLambdaFlag);
int calcGamma(int lambdaIdx, int iter, gammaType *gamma, sparseVector *aBar, sparseMatrix **Bbar, vector pi, BOOL newSigmaFlag, BOOL newLambdaFlag);
int calcSigma(numType *num, coordType *coord, sparseVector *bBar, sparseMatrix *CBar, vector pi, double mubBar,
              int idxLambda, BOOL newLambdaFlag, int iter, sigmaType *sigma, lambdaType *lambda, BOOL *newSigmaFlag);
int computeMu(LPptr lp, int numCols, double *mubBar);
lambdaType *newLambda(int num_iter, int numLambda, int numRVrows, coordType *coord);
sigmaType *newSigma(int numIter, int numNzCols, int numPi, coordType *coord);
deltaType *newDelta(int numIter, coordType *coord);
omegaType *newOmega(int numIter);
void printLambda(lambdaType *lambda, numType *num, int idx);
void freeLambdaType(lambdaType *lambda);
void freeSigmaType(sigmaType *sigma);
void freeGammaType(gammaType *gamma);
void freeOmegaType(omegaType *omega);
void freeDeltaType(deltaType *delta, int numLambda, int numOmega);


#endif /* stocUpdt_h */
