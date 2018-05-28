//
//  cuts.h
//  rollingSP
//
//  Created by Site Wang on 4/23/16.
//  Copyright © 2016 Site Wang. All rights reserved.
//

#ifndef cuts_h
#define cuts_h

#include "rollingSP.h"
#include "soln.h"

/* subroutines in cuts.c */
int formOptCut(probType **p, cellType *c, solnType *s, BOOL newOmegaFlag, int omegaIdx);
cutType *newCuts(int maxCuts);
oneCut *newCut(int numX, int numIstar, int numSamples, int numPeriod);
int SDCut(numType *num, coordType *coord, sigmaType *sigma, deltaType *delta, omegaType *omega, oneCut *cut, vector Xvect, vector piSRatio, BOOL *dualStableFlag, int numSamples);
int tsSDCut(numType *num, coordType *coord, lambdaType *lambda, gammaType *gamma, sigmaType *sigma, deltaType *delta, deltaType *auxDelta, deltaType *trimDelta, omegaType *omega,tsType *arima, oneCut *cut, vector Xvect, vector piSRatio, BOOL *dualStableFlag, int numSamples, int periodStart);
int formFeasCut(probType **prob, cellType *cell, solnType *soln, BOOL *newOmegaFlag, int omegaIdx);
int add2CutPool(cellType *cell, double alpha, vector beta, int betaLen, int numOmega);
int updtFeasCutPool(probType *prob, cellType *cell);
int checkFeasCutPool(cutType *cutPool, cutType *cutsAdded, int betaLen, vector incumbX, vector candidX, BOOL *infeasIncumb);
iType computeIstar(numType *num, coordType *coord, sigmaType *sigma, deltaType *delta, vector Xvect, vector PiCbarX, double *argmaxPre, double *argmaxPost, int obs, int numSamples);
iType r_computeIstar(numType *num, coordType *coord, tsType *arima, omegaType *omega, lambdaType *lambda, gammaType *gamma, sigmaType *sigma, deltaType *delta, deltaType *auxDelta, deltaType *trimDelta, vector Xvect, vector PiCbarX, double *argmaxPre, double *argmaxPost, int obs, int numSamples, int periodStart);
int addCut(numType *num, oneCut *cut, solnType *soln, cellType *c);
int addfCut(int type, LPptr lp, numType *num, int optCuts, vector incumbX, oneCut *cut, int idx);
int formIncumbCut(probType **p, cellType *c, solnType *s, int omegIdx, BOOL newOmegaFlag);
int rFormIncumbCut(probType **p, cellType *c, solnType *s, int omegIdx, vector fullObs, vector auxObs, BOOL newOmegaFlag);
int reduceCuts(cellType *cell, solnType *soln);
int dropCut(LPptr lp, cutType *cuts, solnType *soln, int cutIdx);
void freeCut(oneCut *cut);
void freeCuts(cutType *cuts);


#endif /* cuts_h */
