/*
 * stocUpdt.c
 *
 *  Created on: Oct 21, 2015
 *      Author: Harsha Gangammanavar
 */

#include "stocUpdt.h"

extern runType  run;
extern BOOL     newPeriod;

/* This function obtains a new vector of realizations of the random variables. It compares the new vector with all previous vectors, looking for
 * a duplication.  If it finds a duplicate, it returns the index of that duplicate; otherwise, it adds the vector to the list of distinct realizations
 * and returns the index of that realization. Note that the simulated observation does not have contain one-norm, while the values stored in
 * omegaType do */
int calcOmega(vector observ, vector auxObserv, vector trimObs, int begin, int end, omegaType *omega, omegaType *auxOmega, omegaType *trimOmega, BOOL *newOmegaFlag) {
    
#ifdef TRACE
    printf("\t\t~CalcOmega() \n");
#endif
    
	int cnt;

	/* Compare vector with all the previous observations */
	for (cnt = 0; cnt < omega->cnt; cnt++)
		if (equalVector(observ, omega->vals[cnt], end-begin, run.TOLERANCE)) {
			(*newOmegaFlag) = FALSE;
			omega->weight[cnt]++;
            auxOmega->weight[cnt]++;
            trimOmega->weight[cnt]++;
#ifdef TRACE
            printf("\t :: Updating omega %d weight to %d\n",cnt, omega->weight[cnt]);
#endif
			return cnt;
		}

	/* Add the realization vector to the list */
	omega->vals[omega->cnt] = duplicVector(observ, end-begin);
    auxOmega->vals[auxOmega->cnt] = duplicVector(auxObserv, end-begin);
    trimOmega->vals[trimOmega->cnt] = duplicVector(trimObs, end-begin);
    
	omega->weight[omega->cnt] = 1;
    auxOmega->weight[auxOmega->cnt] = 1;
    trimOmega->weight[trimOmega->cnt] = 1;
	(*newOmegaFlag) = TRUE;

    
    auxOmega->cnt++;
    trimOmega->cnt++;
    
	return omega->cnt++;
}//calcOmega()

/* This function updates all the structures necessary for forming a stochastic cut. The latest observation of omega and the latest dual solution
 * to the subproblem are added to their appropriate structures. Then Pi x b and Pi x C are computed and for the latest omega and dual vector,
 * and are added to the appropriate structures.
 * Note that the new column of delta is computed before a new row in lambda is calculated and before the new row in delta is completed,
 * so that the intersection of the new row and new column in delta is only computed once (they overlap at the bottom, right-hand corner). */
BOOL stochasticUpdates(probType *p, numType *num, vector x, coordType *coord, sparseVector *bBar, sparseMatrix *Cbar, sparseVector *aBar,
                       sparseMatrix **Bbar, lambdaType *lambda, gammaType *gamma, sigmaType *sigma, deltaType *delta, deltaType *auxDelta, deltaType *trimDelta,
                       omegaType *omega, omegaType *auxOmega, omegaType *trimOmega, BOOL newOmegaFlag, int omegaIdx, int maxIter, int iter, vector pi, double mubBar) {
#ifdef TRACE
    printf("\t\t~stochasticUpdates()\n");
#endif

	int 	lambdaIdx, sigmaIdx, gammaIdx;
	BOOL 	newLambdaFlag = FALSE, newSigmaFlag = FALSE;
	

	/* Only need to calculate column if new observation of omega found */
	if (newOmegaFlag)
		calcDeltaCol(num, coord, lambda, omega->vals[omegaIdx], auxOmega->vals[omegaIdx], trimOmega->vals[omegaIdx],
                        omegaIdx, delta, auxDelta, trimDelta);
    
	/* extract the dual solutions corresponding to rows with random elements in them */
	lambdaIdx = calcLambda(num, coord, pi, lambda, sigma, &newLambdaFlag);
    
	/* compute Pi x bBar and Pi x Cbar */
	sigmaIdx = calcSigma(num, coord, bBar, Cbar, pi, mubBar, lambdaIdx, newLambdaFlag, iter, sigma, lambda, &newSigmaFlag);
    
    /* used extracted dual solution to calculate some cut contributions with deterministic time series terms */
    /* If no time series is incorporated, this will turn out to be zero */
    gammaIdx = calcGamma(lambdaIdx, iter, gamma, aBar, Bbar, pi, newSigmaFlag, newLambdaFlag);
    
    /* Capture the associated sigmaIdx for debugging */
    omega->sigmaIdx[omegaIdx] = sigmaIdx;
    
	/* Only need to calculate row if a distinct lambda was found. We could use Pi, instead of lambda(Pi), for this calculation, */
	/* and save the time for expanding/reducing vector. */
	/* even though the lambda is the same, the current Pi might be a distinct one due to the variations in sigma*/
    // Input is specially modified for rolling Horizon problem!!!
	if (newLambdaFlag)
		calcDeltaRow((maxIter*run.HORIZON+run.HORIZON), num, coord, omega, lambda, auxOmega, trimOmega, lambdaIdx, delta, auxDelta, trimDelta);
    
#if defined(CUTS)
    int     idx;
    for(idx=0; idx<lambda->cnt; idx++)
        printLambda(lambda, num, idx);
#endif
    
#if defined(VALID)
    double  objEst;
    objEst = sigma->vals[sigmaIdx].b
                + gamma->endo->pibARMA[lambdaIdx]
                - vXv(sigma->vals[sigmaIdx].C,x,coord->colsC,num->cntCcols)
                + delta->val[lambdaIdx][omegaIdx].b
                + auxDelta->val[lambdaIdx][omegaIdx].b
                + trimDelta->val[lambdaIdx][omegaIdx].b
                - vXv(delta->val[lambdaIdx][omegaIdx].C,x,coord->rvCols, num->rvColCnt);
    printf("\t :: sigmaIdx(%d), lambdaIdx(%d), omegaIdx(%d)\n", sigmaIdx, lambdaIdx, omegaIdx);
	printf("\t :: objective estimate is : %lf\n", objEst);
#endif

	return newSigmaFlag;
}//END stochasticUpdates


/* This function calculates a new column in the delta structure, based on a new observation of omega. Thus, lambda_pi X C and lambda_pi X b
 * are calculated for all values of lambda_pi, for the new C(omega) and b(omega).  Room in the array has already been allocated, so the function
 * only fills it, in the column specified by _obs_. It is assumed that this observation is distinct from all previous ones, and thus a new column
 * must be calculated. */
void calcDeltaCol(numType *num, coordType *coord, lambdaType *lambda, vector observ, vector auxObs, vector trimObs, int omegaIdx,
                  deltaType *delta, deltaType *auxDelta, deltaType *trimDelta) {
#ifdef TRACE
    printf("\t\t\t~calcDeltaCol()\n");
#endif
    
	int piIdx;
    
	sparseVector bomega, auxbomega, trimbomega;
	sparseMatrix Comega, auxComega, trimComega;
	vector lambPi;
	vector piCrossC, piCrossCaux, piCrossCtrim;

	bomega.cnt = num->rvbOmCnt;
    bomega.col = coord->omegaRow;
    bomega.val = observ;

	Comega.cnt = num->rvCOmCnt;
    Comega.col = coord->omegaCol + num->rvbOmCnt;
	Comega.row = coord->omegaRow + num->rvbOmCnt;
    Comega.val = observ + num->rvbOmCnt;
    
    //Newly added section :: Auxilary Observation need to considered as well
    auxbomega.cnt = num->rvbOmCnt;
    auxbomega.col = coord->omegaRow;
    auxbomega.val = auxObs;
    
    auxComega.cnt = num->rvCOmCnt;
    auxComega.col = coord->omegaCol + num->rvbOmCnt;
    auxComega.row = coord->omegaRow + num->rvbOmCnt;
    auxComega.val = auxObs + num->rvbOmCnt;
    
    trimbomega.cnt = num->rvbOmCnt;
    trimbomega.col = coord->omegaRow;
    trimbomega.val = trimObs;
    
    trimComega.cnt = num->rvCOmCnt;
    trimComega.col = coord->omegaCol + num->rvbOmCnt;
    trimComega.row = coord->omegaRow + num->rvbOmCnt;
    trimComega.val = trimObs + num->rvbOmCnt;
    
    
	/* For all dual vectors, lambda(pi), calculate pi X bomega and pi X Comega */
	for (piIdx = 0; piIdx < lambda->cnt; piIdx++) {
		/* Retrieve a new (sparse) dual vector, and expand it into a full vector */
		lambPi = expandVector(lambda->vals[piIdx], coord->rvRows, num->rvRowCnt, num->rows);

		/* Multiply the dual vector by the observation of bomega and Comega */
		/* Reduce PIxb from its full vector form into a sparse vector */
		delta->val[piIdx][omegaIdx].b = vXvSparse(lambPi, &bomega);
        
		piCrossC = vxMSparse(lambPi, &Comega, num->prevCols);
		delta->val[piIdx][omegaIdx].C = reduceVector(piCrossC, coord->rvCols, num->rvColCnt);
        
        // Newly added section
        auxDelta->val[piIdx][omegaIdx].b = vXvSparse(lambPi, &auxbomega);
        
        piCrossCaux = vxMSparse(lambPi, &auxComega, num->prevCols);
        auxDelta->val[piIdx][omegaIdx].C = reduceVector(piCrossCaux, coord->rvCols, num->rvColCnt);
        
        // Newly added section
        trimDelta->val[piIdx][omegaIdx].b = vXvSparse(lambPi, &trimbomega);
        
        piCrossCtrim = vxMSparse(lambPi, &trimComega, num->prevCols);
        trimDelta->val[piIdx][omegaIdx].C = reduceVector(piCrossCtrim, coord->rvCols, num->rvColCnt);
 
        mem_free(piCrossC);
        mem_free(piCrossCaux);
        mem_free(piCrossCtrim);
		mem_free(lambPi);
	}

}//END calcDeltaCol

/* This function stores a new lambda_pi vector in the lambda structure.  Each lambda_pi represents only those dual variables whose rows in the
 * constraint matrix have random elements.  Thus  the (full) dual vector, Pi,  passed to the function is converted into the sparse vector lambda_pi.
 * This vector is then compared with all previous lambda_pi vectors, searching for a duplication. If a duplicate is found, the vector is not added
 * to the structure, and the function returns the index of the duplicate vector. Otherwise, it adds the vector to the end of the structure,
 *and returns an index to the last element in lambda. */
int calcLambda(numType *num, coordType *coord, vector Pi, lambdaType *lambda, sigmaType *sigma, BOOL *newLambdaFlag) {
#ifdef TRACE
    printf("\t\t\t~CalcLambda()\n");
#endif
	int 	pi_idx, i;
	vector	lambda_pi;

	/* Pull out only those elements in dual vector which have rv's */
	lambda_pi = reduceVector(Pi, coord->rvRows, num->rvRowCnt);
    
	/* Compare resulting lambda_pi with all previous vectors */
    for (pi_idx = 0; pi_idx < lambda->cnt; pi_idx++) {
		if (equalVector(lambda_pi, lambda->vals[pi_idx], num->rvRowCnt, run.TOLERANCE)) {
            if (lambda->alive[pi_idx] == DROPPED) {
                // Once a sleeping dual is re-activated, also reactivate's its associated sigmas
                lambda->alive[pi_idx] = SURVIVED;
//                printf("\t\t :: Reactivating lambda %d\n", pi_idx);
                for (i=0; i<sigma->cnt; i++) {
                    // This will increase the computational burden, but it only happen once
                    if (sigma->lambdaIdx[i] == pi_idx) {
                        sigma->linker[sigma->aCnt] = i; //Append the index
                        sigma->aCnt++;  //Counter increment
//                        printf("\t\t\t :: Linking sigma[%d] to end[%d]\n", i, sigma->aCnt);
                    }
                }
            }
			mem_free(lambda_pi);
			*newLambdaFlag = FALSE;
			return pi_idx;
		}
    }
    
	/* Add the vector to lambda structure */
	lambda->vals[lambda->cnt] = lambda_pi;
    lambda->alive[lambda->cnt] = SURVIVED;
	*newLambdaFlag = TRUE;

	return lambda->cnt++;
}//END calcLambda

int calcSigma(numType *num, coordType *coord, sparseVector *bBar, sparseMatrix *CBar, vector pi, double mubBar,
		int idxLambda, BOOL newLambdaFlag, int iter, sigmaType *sigma, lambdaType *lambda, BOOL *newSigmaFlag) {
    
#ifdef TRACE
    printf("\t\t\t~CalcSigma()\n");
#endif
    
	vector	piCBar, temp;
	double 	pibBar;
	int 	cnt, i;

	/* sigma = \pi_t^\top \bar{b}_t - \bar{C}_t^\top \pi_t */
	pibBar = vXvSparse(pi, bBar) + mubBar;
    temp = vxMSparse(pi, CBar, num->prevCols);
	piCBar = reduceVector(temp, coord->colsC, num->cntCcols);
	mem_free(temp);

	if (!newLambdaFlag){
		for (cnt = 0; cnt < sigma->cnt; cnt++) {
			if (DBL_ABS(pibBar - sigma->vals[cnt].b) <= run.TOLERANCE) {
				if (equalVector(piCBar, sigma->vals[cnt].C, num->cntCcols, run.TOLERANCE))
                if (sigma->lambdaIdx[cnt] == idxLambda){
                    // Retriving a sigma when its linked lambda is dropped
                    if (lambda->alive[idxLambda] == DROPPED) {
                        lambda->alive[idxLambda] = SURVIVED;    //Brought back the lambda
                        // Once detected, perform a through search on sleeping sigmas
                        // But this only happens once, this hurts the complexity
                        for (i = 0; i<sigma->cnt; cnt++) {
                            if (sigma->lambdaIdx[i] == idxLambda) {
                                sigma->linker[sigma->aCnt] = i;
                                sigma->aCnt++;
                            }
                        }
                    }
                    mem_free(piCBar);
                    (*newSigmaFlag) = FALSE;
                    return cnt;
				}
			}
		}
	}

	(*newSigmaFlag) = TRUE;
	sigma->vals[sigma->cnt].b  = pibBar;
	sigma->vals[sigma->cnt].C  = piCBar;
	sigma->lambdaIdx[sigma->cnt] = idxLambda;
	sigma->ck[sigma->cnt] = iter;
    
    sigma->linker[sigma->aCnt] = sigma->cnt; //Record the current index for link
    sigma->aCnt++;
    
	return sigma->cnt++;

}//END calcSigma()

/* This function calculates a new row in the delta structure, based on a new dual vector, lambda_pi, by calculating lambda_pi X b and
 * lambda_pi X C for all previous realizations of b(omega) and C(omega).  It is assumed that the lambda vector is distinct from all previous ones
 * and thus a new row is warranted. */
int calcDeltaRow(int maxIter, numType *num, coordType *coord, omegaType *omega, lambdaType *lambda, omegaType *auxOmega, omegaType *trimOmega,
                 int lambdaIdx, deltaType *delta, deltaType *auxDelta, deltaType *trimDelta) {
#ifdef TRACE
    printf("\t\t\t~CalcDeltaRow()\n");
#endif
	sparseVector bomega, auxbomega, trimbomega;
	sparseMatrix Comega, auxComega, trimComega;
	vector 	lamb_pi, pixC, pixCaux, pixCtrim;
	int		obs;

	bomega.cnt = num->rvbOmCnt;
    bomega.col = coord->omegaRow;
    
    Comega.cnt = num->rvCOmCnt;
    Comega.col = coord->omegaCol + num->rvbOmCnt;
    Comega.row = coord->omegaRow + num->rvbOmCnt;

	if ( !(delta->val[lambdaIdx] = (pixbCType *) arr_alloc(maxIter, pixbCType)))
		errMsg("allocation", "calcDeltaRow", "delta->val[cnt]", 0);
    
    //Newly added section
    auxbomega.cnt = num->rvbOmCnt;
    auxbomega.col = coord->omegaRow;
    
    auxComega.cnt = num->rvCOmCnt;
    auxComega.col = coord->omegaCol + num->rvbOmCnt;
    auxComega.row = coord->omegaRow + num->rvbOmCnt;
    
    if ( !(auxDelta->val[lambdaIdx] = (pixbCType *) arr_alloc(maxIter, pixbCType)))
        errMsg("allocation", "calcDeltaRow", "delta->val[cnt]", 0);
    
    //Newly added section for trim observations
    trimbomega.cnt = num->rvbOmCnt;
    trimbomega.col = coord->omegaRow;
    
    trimComega.cnt = num->rvCOmCnt;
    trimComega.col = coord->omegaCol + num->rvbOmCnt;
    trimComega.row = coord->omegaRow + num->rvbOmCnt;
    
    if ( !(trimDelta->val[lambdaIdx] = (pixbCType *) arr_alloc(maxIter, pixbCType)))
        errMsg("allocation", "calcDeltaRow", "delta->val[cnt]", 0);

 
	/* expand the compressed lambda vector */
	lamb_pi = expandVector(lambda->vals[lambdaIdx], coord->rvRows, num->rvRowCnt, num->rows);

	/* go through all the observations and compute pi x b and pi x C */
	for (obs = 0; obs < omega->cnt; obs++) {
        // Original deltaType structure
		bomega.val= omega->vals[obs];
		Comega.val = omega->vals[obs] + num->rvbOmCnt;
		delta->val[lambdaIdx][obs].b = vXvSparse(lamb_pi, &bomega);
		pixC = vxMSparse(lamb_pi, &Comega, num->prevCols);
		delta->val[lambdaIdx][obs].C = reduceVector(pixC, coord->rvCols, num->rvColCnt);
        
        // Aux-deltaType structure
        auxbomega.val = auxOmega->vals[obs];
        auxComega.val = auxOmega->vals[obs] + num->rvbOmCnt;
        auxDelta->val[lambdaIdx][obs].b = vXvSparse(lamb_pi, &auxbomega);
        pixCaux = vxMSparse(lamb_pi, &auxComega, num->prevCols);
        auxDelta->val[lambdaIdx][obs].C = reduceVector(pixCaux, coord->rvCols, num->rvColCnt);
        
        // Trim-deltaType structure
        trimbomega.val = trimOmega->vals[obs];
        trimComega.val = trimOmega->vals[obs] + num->rvbOmCnt;
        trimDelta->val[lambdaIdx][obs].b = vXvSparse(lamb_pi, &trimbomega);
        pixCtrim = vxMSparse(lamb_pi, &trimComega, num->prevCols);
        trimDelta->val[lambdaIdx][obs].C = reduceVector(pixCtrim, coord->rvCols, num->rvColCnt);
        
        mem_free(pixC); mem_free(pixCaux); mem_free(pixCtrim);
	}

	mem_free(lamb_pi);

	return 0;

}//END calcDeltaRow()

/* This function compute the reduced cost of every second stage variables. They will be used to calculate the \mu x b and then added to the \pi x b. */
int computeMu(LPptr lp, int numCols, double *mubBar) {
	vector	dj, u;
	intvec	cstat;
	int		n;

	(*mubBar) = 0.0;

	if ( !(dj = (vector) arr_alloc(numCols+1, double)))
		errMsg("allocation", "computeMu", "dual slacks", 0);
	if ( !(u = (vector) arr_alloc(numCols+1, double)))
		errMsg("allocation", "computeMu", "TDA solutions", 0);

	if (getPrimal(lp, u, numCols) ) {
		errMsg("solver", "forOptPass", "failed to obtain primal solution", 0);
		return 1;
	}
	if (getDualSlacks(lp, dj, numCols) ) {
		errMsg("solver", "computeMu", "failed to obtain dual slacks", 0);
		return 1;
	}

	/* extra column for eta if the stage problem is a QP */
	if ( !(cstat = (intvec) arr_alloc(numCols+2, int)) )
		errMsg("allocation", "computeMu", "column status", 0);
	if (getBasis(lp, cstat+1, NULL)) {
		errMsg("solver", "computeMu", "failed to get column status", 0);
		return 1;
	}

	for (n = 1; n <= numCols;  n++) {
		switch (cstat[n]) {
		case AT_LOWER:
			(*mubBar) += dj[n]*u[n];
			break;
		case AT_UPPER:
			(*mubBar) += dj[n]*u[n];
			break;
		default:
			break;
		}
	}

	mem_free(u); mem_free(cstat); mem_free(dj);

	return 0;
}//END of computeMu()

/* This function allocates a new lambda structure, with room for num_lambdas lambda vectors of size vect_size.  It returns a pointer to the structure.
 * Only some of the individual lambda vectors are expected to be allocated (according to the num_vect parameter) so that there is room for new
 * lambdas to be created. */
lambdaType *newLambda(int num_iter, int numLambda, int numRVrows, coordType *coord) {
    
	lambdaType *lambda;
	int cnt;

	if (!(lambda = (lambdaType *) mem_malloc(sizeof(lambdaType))))
		errMsg("allocation", "newLambda", "lambda",0);
    
    if (!(lambda->alive = arr_alloc(num_iter, int)))
        errMsg("allocation", "newLambda", "lambda", 0);

	if (!(lambda->vals = arr_alloc(num_iter, vector)))
		errMsg("allocation", "newLambda", "lambda->val",0);

	for (cnt = 0; cnt < numLambda; cnt++)
		if (!(lambda->vals[cnt] = arr_alloc(numRVrows+1, double)))
			errMsg("allocation", "newLambda", "lambda->val[cnt]",0);

	lambda->cnt = numLambda;

	return lambda;
}//END of newLambda

/* This function creates a new sigma structure, and allocates memory for the arrays associated with it.  It returns a pointer to this structure.
 * Some pi X T vectors are also allocated, according to the num_vals parameter  (num_vals is expected to be less than num_sigmas, so that there
 * is room for further work).  Note that  memory for sigma->col is not allocated, but is taken from prob.*/
sigmaType *newSigma(int numIter, int numNzCols, int numPi, coordType *coord) {
	sigmaType *sigma;
	int cnt;

	if (!(sigma = (sigmaType *) mem_malloc (sizeof(sigmaType))))
		errMsg("allocation", "newSigma", "sigma",0);

	if (!(sigma->lambdaIdx = (intvec) arr_alloc(numIter, int)))
		errMsg("allocation", "new_sigma", "sigma->lamb",0);

	if (!(sigma->ck = (intvec) arr_alloc(numIter, int)))
		errMsg("allocation", "new_sigma", "sigma->ck",0);

    if (!(sigma->linker = (intvec) arr_alloc(numIter, int)))
        errMsg("allocation", "new_sigma", "sigma->ck", 0);
    
	if (!(sigma->vals = arr_alloc(numIter, pixbCType)))
		errMsg("allocation", "new_sigma", "sigma->val",0);

	for (cnt = 0; cnt < numPi && cnt < numIter; cnt++)
		if (!(sigma->vals[cnt].C = arr_alloc(numNzCols+1, double)))
			errMsg("allocation", "new_sigma", "sigma->val[cnt]",0);

	sigma->cnt = numPi;
    sigma->aCnt = numPi;

	return sigma;
}//END newSigma

/***********************************************************************\
 ** This function creates a new delta structure with arrays of the specified
 ** size and returns a pointer to it.  Note that the pi X T vectors
 ** themselves are not allocated, since they will not all be filled with
 ** values.  (they are only filled as they are produced).
 ** Not even the arrays of pi_R_T_types are allocated, as this also
 ** occurs in calc_delta_row().  However, the column coordinates of the
 ** eventual multiplications are initialized, since they are known.
 \***********************************************************************/
deltaType *newDelta(int numIter, coordType *coord) {
	deltaType *d;

	if (!(d = (deltaType *) mem_malloc (sizeof(deltaType))))
		errMsg("Allocation", "newDelta", "d",0);
	if (!(d->val = (pixbCType **) arr_alloc(numIter, pixbCType *)))
		errMsg("Allocation", "newDelta", "d->val",0);

	return d;
}//END newDelta

/* This function allocates memory for an omega structure.  It allocates the memory to structure elements: a vector to hold an array of
 * observation and the weights associated with it. */
omegaType *newOmega(int numIter) {
	omegaType *omega;

	if ( !(omega = (omegaType *) mem_malloc(sizeof(omegaType))) )
		errMsg("allocation","newOmega", "omega", 0);
	if ( !(omega->weight = (intvec) arr_alloc(numIter, int)) )
		errMsg("allocation", "newOmega", "omega->weight", 0);
	if ( !(omega->vals = (vector *) arr_alloc(numIter, vector)) )
		errMsg("allocation", "newOmega", "omega->vals", 0);
    if ( !(omega->sigmaIdx = (intvec) arr_alloc(numIter, int)))
        errMsg("allocation", "newOmega", "omega->sigmaIdx", 0);
	omega->cnt = 0;

	return omega;
}//END of newOmega()

void printLambda(lambdaType *lambda, numType *num, int idx) {
    int cnt;
    
    printf("\nLambda %d:: ", idx);
    for (cnt = 1; cnt <= num->numRV; cnt++) {
        printf("%f  ", lambda->vals[idx][cnt]);
    }
    printf("\n");
    
}//END of printLambda

/*This function releases all the memory reserved for the Lambda data structure, including the Lambda itself. */
void freeLambdaType(lambdaType *lambda) {
	int cnt;
	if ( lambda ) {
        mem_free(lambda->alive);
		for (cnt = 0; cnt < lambda->cnt; cnt++) {
            if (lambda->vals[cnt]) mem_free(lambda->vals[cnt]);
        }
		if ( lambda->vals )
            mem_free(lambda->vals);
		mem_free(lambda);
    }
}
//END freeLambdaType()

/*  This function releases all the memory reserved for the sigma data structure, including the sigma itself.*/
void freeSigmaType(sigmaType *sigma) {
	int cnt;

	if ( sigma ) {
		for (cnt = 0; cnt < sigma->cnt; cnt++) {
			if (sigma->vals[cnt].C) mem_free(sigma->vals[cnt].C);
		}
		if ( sigma->lambdaIdx) mem_free(sigma->lambdaIdx);
        if ( sigma->linker ) mem_free(sigma->linker);
		if ( sigma->ck) mem_free(sigma->ck);
		if ( sigma->vals) mem_free(sigma->vals);
		mem_free(sigma);
	}
}
//END freeSigmaType()

void freeOmegaType(omegaType *omega) {
	int n;

	if (omega->weight) mem_free(omega->weight);
    if (omega->sigmaIdx) mem_free(omega->sigmaIdx);
	if (omega->vals) {
		for (n = 0; n < omega->cnt; n++ )
			if (omega->vals[n]) mem_free(omega->vals[n]);
		mem_free(omega->vals);
	}
	mem_free(omega);

}//END freeOmegaType()

void freeDeltaType(deltaType *delta, int numLambda, int numOmega) {
	int m, n;

	if (delta->val) {
		for ( m = 0; m < numLambda; m++ ) {
			for ( n = 0; n < numOmega; n++ )
				if ( delta->val[m][n].C ) mem_free(delta->val[m][n].C);
			mem_free(delta->val[m]);
		}
		mem_free(delta->val);
	}
	mem_free(delta);

}//END freeDeltaType()
