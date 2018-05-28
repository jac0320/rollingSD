//
//  roll.h
//  rollingSP
//
//

#ifndef roll_h
#define roll_h

#include "rollingSP.h"
#include "stocUpdt.h"
#include "master.h"
#include "cuts.h"

int setupDetPeriod(oneProblem *detProblem, cellType *cell, probType **prob, solnType *soln, BOOL *initialPeriod);

int setupStoPeriod(cellType *cell, stocType *stoc, probType **prob, solnType *soln, BOOL *initialPeriod);
int refreshPeriod(cellType *cell, stocType *stoc, probType **prob, solnType *soln, BOOL *initialPeriod);
int warmUpTimeSeries(tsType *arima, int current);
int warmUpGamma(cellType *cell, gammaType *gamma, tsType *arima, int t, lambdaType *lambda, stocType *stoc);
int warmUpDynamic(probType **prob, solnType *soln, cellType *cell, oneProblem *master, oneProblem *subprob);
int warmStartCell(cellType *cell, probType **prob, solnType *soln);
int transformCuts(oneProblem *master, cellType *cell, int mCols, cutType *cuts, solnType *soln, int k);
int refreshSoln(solnType *soln, cellType *cell, probType **prob);

#endif /* roll_h */
