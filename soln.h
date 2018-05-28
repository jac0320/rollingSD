//
//  soln.h
//  rollingSP
//
//  Created by Site Wang on 4/23/16.
//  Copyright © 2016 Site Wang. All rights reserved.
//

#ifndef soln_h
#define soln_h

#include "rollingSP.h"
#include "cell.h"
#include "cuts.h"

/* solnType is defined the main header file since it is a critical structure */

/* subroutines in soln.c */
BOOL checkImprovement(probType **prob,cellType *cell, solnType *soln);
int replaceIncumbent(probType **prob, cellType *cell, solnType *soln, double newCandidEst);
int newIncumbent(probType *p, cellType *c, solnType *s, double est);
int resolveInfeasibility(probType **prob, cellType *cell, solnType *soln, BOOL *newOmegaFlag, int omegaIdx);
double maxCutHeight(int lbType, double lb, cutType *cuts, int currIter, vector xk, int betaLen, int *maxCutIdx);
double cutHeight(int lbType, double lb, oneCut *cut, int currIter, vector xk, int betaLen);
solnType * newSoln(probType **p, int maxCuts, vector xk, double lb, int usefulProbIdx);
void freeSolnType(solnType *soln);
void freeRunTimeType(runTimeType *runTime);


#endif /* soln_h */
