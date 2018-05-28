//
//  evaluate.h
//  rollingSP
//
//  Created by Site Wang on 4/24/16.
//  Copyright © 2016 Site Wang. All rights reserved.
//

#ifndef evaluate_h
#define evaluate_h

#include "cell.h"
#include "rollingSP.h"
#include "subprob.h"


/* Subroutines in evaluate.c */
void evalStoOpt(stocType *stoc, probType **prob, cellType *cell, solnType *soln, FILE *SubOptVec, FILE *ePtr);
void evalDetOpt(oneProblem *orig, timeType *tim, stocType *stoc, probType **prob, cellType *cell, solnType *soln, FILE *SubOptVec, FILE *ePtr);
double calcVariance(vector x, vector mean_value, vector stdev_value, int batch_size);
void evalDualStability(vector piSRatio, double argmaxPreSum, double argmaxAllSum, BOOL *dualStableFlag, int numSamples);

#endif /* evaluate_h */
