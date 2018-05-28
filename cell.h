//
//  cell.h
//  rollingSP
//
//  Created by Site Wang on 4/27/16.
//  Copyright © 2016 Site Wang. All rights reserved.
//

#ifndef cell_h
#define cell_h

#include "rollingSP.h"
#include "gamma.h"
#include "arma.h"

typedef struct {
    int			t;				// time period of current cell
    int 		k;              // Iteration Number
    int			minIter;		// Minimum number of iterations
    int			maxIter;		// Maximum number of iterations
    int			tau;			// frequency of incumbent updates
    int			lbType;			// lower bound type
    double		lb;				// lower bound
    int			feasCnt;		// counter to keep track of infeasible subproblems */
    int			LPcnt;			// Number of linear programs solved
    int         Kt;             // Final Iteration that Last Period on
    int         Kc;             // Cumulative iteration counts
    oneProblem 	*master;     	// Master Problem
    oneProblem 	*subprob;    	// Subproblem
    vector 		initPrimal;
    lambdaType 	*lambda;		/* holds dual solutions corresponding to rows effected by randomness */
    sigmaType 	*sigma;			/* holds $\pi \times \bar{b}$ and $\pi \times \bar{C} $ values */
    deltaType   *delta;			/* calculations based on realization and dual solutions observed */
    deltaType   *auxDelta;      /* Delta type used with time series model */
    deltaType   *trimDelta;
    omegaType 	*omega;			/* all realizations observed during the algorithm */
    omegaType   *auxOmega;
    omegaType   *trimOmega;
    gammaType   *gamma;         /* All about dynamic/stochasticity  */
    tsType      *arima;
    double		quadScalar;
    cutType 	*cuts;          // A set of cuts generated
    int			maxCuts;
    cutType 	*fCutPool;		// Pool of feasibility cuts
    cutType		*fCutAdded;		// set of added feasibility cuts
    int			fCol;			// column index for delta structure for which feasibility cut has been computed
    int			fRow;			// row index for delta structure for which feasibility cut has been computed
}cellType;

// cell.c
int solveSDCell(stocType *stoc, probType **prob, cellType *cell, solnType *soln);
int solveDETCell(cellType *cell, probType *probM, solnType *soln);
int injectSolution(cellType *cell, probType *probM, solnType *soln);
cellType *newCell(probType **prob, stocType *stoc);
void freeCellType(cellType *cell);

// Algorithm specific random generator
int generateARIMA(int t, cellType *cell, vector obs, vector noise, vector auxObs, vector trimObs, vector inputNoise, long long *seed);
void generateOmega(cellType *cell, stocType *stoc, vector residual, vector fullObs, vector auxObs, vector trimObs, vector inputNoise, long long *seed);


#endif /* cell_h */
