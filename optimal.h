//
//  optimal.h
//  rollingSP
//
//  Created by Site Wang on 4/10/16.
//  Copyright © 2016 Site Wang. All rights reserved.
//

#ifndef optimal_h
#define optimal_h

#include <time.h>
#include "rollingSP.h"
#include "cell.h"
#include "subprob.h"
#include "soln.h"
#include "cuts.h"

int lowerBoundCheck(probType **prob, cellType *cell, solnType *soln, BOOL *check);
int sampleCut(int cutCnt);

BOOL optimal(probType **prob, cellType *cell, solnType *soln);
BOOL preCheck (solnType *soln);
BOOL preTestII(probType **prob, cellType *cell, solnType *soln);
BOOL fullCheck (probType **prob, cellType *cell, solnType *soln);

cutType *chooseCuts(probType *prob, cellType *cell, solnType *soln);
void sampleOmega(int *cdf, int *observ, int k);
int randFun(int greatest);
void empiricalDistrib(omegaType *omega, int *cdf);
BOOL validOmegaIdx(omegaType *omega, int idx);
void reformCuts(cellType *cell, numType *mNum, coordType *coord, cutType *cut, intvec observ, int k);

double solveTempMaster(probType *prob, cutType *cuts, cellType *cell);
double calTempLB(probType *prob, cellType *cell, solnType *soln, cutType *T);
double calculateULm(probType *prob, cellType *cell, cutType *T, solnType *soln);

#endif /* optimal_h */
