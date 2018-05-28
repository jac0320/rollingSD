//
//  master.h
//  rollingSP
//
//  Created by Site Wang on 4/27/16.
//  Copyright © 2016 Site Wang. All rights reserved.
//

#ifndef master_h
#define master_h

#include "cell.h"
#include "cuts.h"
#include "soln.h"

/* subroutines in master.c */
int constructQP(LPptr lp, int numCols, double sigma);
int changeQPrhs(probType *prob, cellType *cell, vector xk);
int changeQPbds(LPptr lp, int numCols, vector bdl, vector bdu, vector xk);
oneProblem *newMaster(oneProblem *master, vector initSol, sparseMatrix *initCbar, cutType *cuts, int extra_cuts);
oneProblem *newStoMaster(oneProblem *master, sparseMatrix *initCbar, cutType *cuts, int extra_cuts);
oneProblem *newDETMaster(oneProblem *detProb);
void freeMaster(oneProblem *copy);
int changeEtaCol(LPptr lp, int numRows, int numCols, int k, cutType *cuts, double lb);
int solveQPMaster(numType *num, sparseVector *cBar, cellType *cell, solnType *soln);
int updateRHS(cellType *cell);

#endif /* master_h */
