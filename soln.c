/*
 * soln.c
 *
 *  Created on: Oct 21, 2015
 *      Author: Harsha Gangammanavar
 */


#include "soln.h"
#include "master.h"
#include "subprob.h"

extern runType run;

/***********************************************************************\
 ** This function determines whether the "stagewise descent property" is
 ** satisified.  If the current approximation of f_k gives a lower difference
 ** between the candidate and incumbent x than the previous approximation
 ** gave, then the incumbent x is updated to the candidate x, and the
 ** reference to the incumbent cut is updated as well.  The function returns
 ** TRUE if the incumbent was updated; FALSE otherwise.
 \***********************************************************************/
BOOL checkImprovement(probType **prob, cellType *cell, solnType *soln) {

#ifdef TRACE
    printf("\t~checkImprovement()\n");
#endif
	double candidEst;

	/* Calculate height at new candidate x with newest cut included */
	candidEst = maxCutHeight(cell->lbType, cell->lb, cell->cuts, cell->k, soln->candidX, prob[0]->num->cols, NULL);
	candidEst += vXvSparse(soln->candidX, prob[0]->dBar);

	/* Calculate height at current incumbent x with newest cut. */
	soln->incumbEst = maxCutHeight(cell->lbType, cell->lb, cell->cuts, cell->k, soln->incumbX, prob[0]->num->cols, NULL);
	soln->incumbEst += vXvSparse(soln->incumbX, prob[0]->dBar);

	/* If we saw considerable improvement, then change the incumbent */
	if ((candidEst - soln->incumbEst) < (run.R1 * soln->gamma)) {
		replaceIncumbent(prob, cell, soln, candidEst);
		return TRUE;
	} else {
		/* Update quad_scalar when no incumbent is found. */
		cell->quadScalar = min(run.MAX_QUAD_SCALAR, cell->quadScalar / run.R2);
#ifdef TRACE
    printf("\t :: No incumbent Found... Updating quadrtic scalar = %f\n", cell->quadScalar);
#endif
		if (cell->master->type == PROB_QP)
			constructQP(cell->master->lp, prob[0]->num->cols, cell->quadScalar);

		soln->incumbStdev *= (cell->k - 1) / (double) (cell->k);
		soln->normDk_1 = soln->normDk;
		return FALSE;
	}

}//END of checkImprovement()

int replaceIncumbent(probType **prob, cellType *cell, solnType *soln, double newCandidEst){

#ifdef TRACE
    printf("\t\t\t~replaceIncumbend()\n");
#endif

	int status;

	/* replace the incumbent solution with the candidate solution */
    status = newIncumbent(prob[1], cell, soln, newCandidEst);
	if ( status ) {
		errMsg("algorithm", "replaceIncumbent", "failed to replace incumbent", 0);
		return 1;
	}

	/* update the proximal parameter based on estimated improvement */
	if ( cell->k > 1 && soln->normDk > run.TOLERANCE )
		if ( soln->normDk >= run.R3 * soln->normDk_1 ) {

			cell->quadScalar *= run.R2 * run.R3 * soln->normDk_1/ soln->normDk;
			cell->quadScalar = min(run.MAX_QUAD_SCALAR, cell->quadScalar);
			cell->quadScalar = max(run.MIN_QUAD_SCALAR, cell->quadScalar);
#if defined(TRACE)
            printf("Quadratic scalar is updated to %f\n", cell->quadScalar);
#endif
			if (cell->master->type == PROB_QP) {
				status = constructQP(cell->master->lp, prob[0]->num->cols, cell->quadScalar);
				if ( status ) {
					errMsg("algorithm", "replaceIncumbent", "failed to change the proximal term after incumbent change", 0);
					return 1;
				}
			}
		}

	/* update the right-hand side and the bounds with new incumbent solution */
	if (cell->master->type == PROB_QP) {
		status = changeQPrhs(prob[0], cell, soln->incumbX);
		if ( status ) {
			errMsg("algorithm", "replaceIncumbent", "failed to change the right-hand side after incumbent change", 0);
			return 1;
		}
		status = changeQPbds(cell->master->lp, prob[0]->num->cols, prob[0]->sp->bdl, prob[0]->sp->bdu, soln->incumbX);
		if ( status ) {
			errMsg("algorithm", "replaceIncumbent", "failed to change the bounds after incumbent update", 0);
			return 1;
		}
	}

    printf("+"); fflush(stdout);

	/* keep the two norm of solution*/
	soln->normDk_1 = soln->normDk;
	/* Since incumbent solution is now replaced by a candidate, we assume it is feasible now */
	soln->infeasIncumb = FALSE;
	/* gamma needs to be reset to 0 since there's no difference between candidate and incumbent*/
	soln->gamma = 0.0;

	return 0;
}//END of replaceIncumbent()

/* If it has been decided that the incumbent solution should change, several values must be updated,
 * including the
 * 1. incumbent X,
 * 2. incumbent cut,
 * 3. estimate at the incumbent
 * 4. the standard deviation of the estimate,
 * 5. Last update counter
 * This function performs these updates.
 * It assumes that the candidate solution becomes the incumbent solution, with a height corresponding
 * to the _est_ parameter.
 */
int newIncumbent(probType *prob, cellType *c, solnType *s, double est) {

#ifdef TRACE
    printf("\t\t\t\t~newIncumbend()\n");
#endif

	double val, stdev;
	int obs, idx, count, incumbCut;
	iType i;

	/* Copy over information about the new incumbent */
	copyVector(s->candidX, s->incumbX, prob->num->prevCols, 1);
	s->incumbEst = est;
	/* In case the incumbent is changed due to infeasibility then need to be sure that cut idx is fine */
    incumbCut = s->iCutIdx = s->cCutIdx;
	s->iCutUpdt = c->k;

	/* The cut consists of an average of many pi X omega products.  We use the mean and standard deviation
	 * of this collection of products in the check for optimality, so they are updated here.  We loop
	 * through the iStar indices of new incumbent cut to find the stdev. */
	count = 0;
	stdev = 0.0;
	if (c->cuts->cnt > 0) {
		for (obs = 0; obs < c->cuts->val[incumbCut]->omegaCnt; obs++) {
			/* Find the pi from the cut's argmax for this observation */
			i.sigma = c->cuts->val[incumbCut]->iStar[obs];
			i.delta = c->sigma->lambdaIdx[i.sigma];

			/* Calculate the height for this observation and dual vector */
            val = c->sigma->vals[i.sigma].b + c->delta->val[i.delta][obs].b;
            if (c->arima != NULL){
                val += c->gamma->endo->pibARMA[c->sigma->lambdaIdx[i.sigma]];
                val += c->auxDelta->val[i.delta][obs].b;
                val += c->trimDelta->val[i.delta][obs].b;
            }
			for (idx = 0; idx < prob->num->cntCcols; idx++)
				val -= c->sigma->vals[i.sigma].C[idx] * s->incumbX[prob->coord->colsC[idx]];
			for (idx = 0; idx < prob->num->rvCOmCnt; idx++)
				val -= c->delta->val[i.delta][obs].C[idx]* s->incumbX[prob->coord->rvCols[idx]];

			stdev += c->omega->weight[obs] * (s->incumbEst - val)*(s->incumbEst - val);
			count += c->omega->weight[obs];
		}
		s->incumbStdev = sqrt(stdev / (double) count);
	}

	s->incumbChg = TRUE;
	return 0;

}//END newIncumbent

/* This function takes the SD code into Feasibility mode (solve_cell() take the SD into Optimality mode). The SD will not return to optimality mode
 until the candidate and incumbent solution are both feasible. */
int resolveInfeasibility(probType **prob, cellType *cell, solnType *soln, BOOL *newOmegaFlag, int omegaIdx) {

#ifdef TRACE
    printf("\t\t~resolveInfeasibility()\n");
#endif

	int status;

	/* QP master will be solved in feasibility mode */
	soln->optMode = FALSE;

	while ( TRUE ) {
		/* form a feasibility cut */
		formFeasCut(prob, cell, soln, newOmegaFlag, omegaIdx);

		/* relax the proximal term and change it in the solver */
		cell->quadScalar = run.MIN_QUAD_SCALAR;
		status = constructQP(cell->master->lp, prob[0]->num->cols, cell->quadScalar);
		if ( status ) {
			errMsg("algorithm", "resolveInfeasibility", "failed to change the proximal parameter", 0);
			return 1;
		}

		/* Solver the master problem with the added feasibility cut */
		status = solveQPMaster(prob[0]->num, prob[0]->dBar, cell, soln);
		if ( status ) {
			errMsg("algorithm", "resolveInfeasibility", "failed to solve the master problem", 0);
			return 1;
		}

		/* increment the count for number of infeasible master solutions encountered */
		cell->feasCnt++;

		status = solveSubprob(prob[1], cell->subprob, soln, soln->candidX, cell->omega->vals[omegaIdx], &soln->mubBar, &soln->subFeasFlag);
		if ( status ) {
			errMsg("algorithm", "resolveInfeasibility", "failed to solve the subproblem", 0);
			return 1;
		}

		/* end the feasibility mode if a feasible candidate solution is observed */
		if (soln->subFeasFlag == TRUE)
			break;
	}

	if ( soln->infeasIncumb == TRUE )
		/* if the incumbent solution is infeasible then replace the incumbent with the feasible candidate solution */
		replaceIncumbent(prob, cell, soln, 0.0);

	/* QP master will be solved in optimality mode again */
	soln->optMode = TRUE;
	return 0;

}//END resolveInfeasibility()

/* This function loops through a set of cuts and find the highest cut height at the specified position x */
double maxCutHeight(int lbType, double lb, cutType *cuts, int currIter, vector xk, int betaLen, int *maxCutIdx) {
	double Sm;
	double ht;
    int cnt = 0, mem = 0;

	Sm = cutHeight(lbType, lb, cuts->val[0], currIter, xk, betaLen);
	for (cnt = 1; cnt < cuts->cnt; cnt++) {
		ht = cutHeight( lbType, lb, cuts->val[cnt], currIter, xk, betaLen);
        if (Sm < ht) {
			Sm = ht;
            mem = cnt;
        }
	}

    if (maxCutIdx)
        *maxCutIdx = mem;

	return Sm;
}//END maxCutHeight

/* This function calculates and returns the height of a given cut at a given X.  It includes the k/(k-1) update, but does not include
 * the coefficients due to the cell. */
double cutHeight(int lbType, double lb, oneCut *cut, int currIter, vector xk, int betaLen) {
	double height;
	double t_over_k = ((double) cut->cutObs / (double) currIter);

	/* A cut is calculated as alpha - beta x X */
	height = cut->alpha - vXv(cut->beta, xk, NULL, betaLen);

	/* Weight cut based on number of observations used to form it */
	height *= t_over_k;

	/* Updated for optimality cut height*/
	if (lbType == NONTRIVIAL)
		height += (1 - t_over_k) * lb;

	return height;
}//END cutHeight()

/* This function allocates memory for a new solnType data structure. It's fields are also allocated, and initialized. */
solnType *newSoln(probType **p, int maxCuts, vector xk, double lb, int usefulProbIdx) {

    solnType *s = NULL;

	if ( !(s = (solnType *) mem_malloc(sizeof(solnType))) )
		errMsg("allocation", "newSoln", "s", 0);

	/* solutions corresponding to master */
	if ( !(s->piM = (vector) arr_alloc(p[0]->num->rows + maxCuts + 1, double)) )
		errMsg("allocation", "newSoln", "s->piM", 0);
	if ( !(s->djM = (vector) arr_alloc(p[0]->num->cols + 2, double)) )
		errMsg("allocation", "newSoln", "s->djM", 0);
    s->meanX = duplicVector(xk, p[0]->num->cols);
	s->candidX 	= duplicVector(xk, p[0]->num->cols);
	s->incumbX 	= duplicVector(xk, p[0]->num->cols);
    s->iCutIdx = 0;
    s->cCutIdx = 0;
	s->iCutUpdt = 0;

	/* solutions corresponding to subproblem */
	if ( !(s->piS = (vector) arr_alloc(p[usefulProbIdx]->num->rows+1, double)) )
		errMsg("allocation", "newSoln", "s->piS", 0);
    if ( !(s->piSRatio = (vector) arr_alloc(run.SCAN_LEN, double)) )
        errMsg("allocation", "newSoln", "s->piSRatio", 0);

	/* elements necessary to update incumbent */
	s->gamma = 0.0;
	s->candidEst = lb + vXvSparse(s->candidX, p[0]->dBar);
	s->incumbEst = s->candidEst;
	s->incumbStdev = 0.0;
    s->FTError = 0.0;
    s->repPassed = 0;
    s->dualStableIter = 0;

	/* Initialize flags */
    s->subFeasFlag      = TRUE;
    s->optMode          = TRUE;
	s->optFlag          = FALSE;
	s->incumbChg        = FALSE;
    s->preCheckEverFlag = FALSE;
    s->preCheckFlag     = FALSE;

    /* Allocation Memory to runTimeType in solnType */
    if ( !(s->runTime = (runTimeType *) mem_malloc(sizeof(runTimeType))) )
        errMsg("allocation", "newSoln", "runTime", 0);

    s->runTime->allTime = 0.0;
    s->runTime->masterCnt = 1;
    s->runTime->optCnt = 1;
    s->runTime->cutGenCnt = 1;
    s->runTime->subprobCnt = 1;
    s->runTime->dropCutCnt = 1;
    s->runTime->detSolveCnt = 1;
    s->runTime->argmaxCnt = 1;

    s->runTime->iterIndex = NULL;

    /* The following number are recorded every period */
    s->runTime->evalOptTime = (vector) arr_alloc(run.HORIZON + 1, double);
    s->runTime->solveTime = (vector) arr_alloc(run.HORIZON + 1, double);
    s->runTime->warmUpTime = (vector) arr_alloc(run.HORIZON + 1, double);

    /* These numbers are recorded every iteration */
    s->runTime->iterIndex = (intvec) arr_alloc(run.MAX_ITER * run.HORIZON + 1, int);
    s->runTime->masterTime = (vector) arr_alloc(run.MAX_ITER * run.HORIZON + 1, double);
    s->runTime->optTime = (vector) arr_alloc(run.MAX_ITER * run.HORIZON + 1, double);
    s->runTime->fullCheck_lb = (vector) arr_alloc(run.MAX_ITER * run.HORIZON + 1, double);
    s->runTime->fullCheck_reform = (vector) arr_alloc(run.MAX_ITER * run.HORIZON + 1, double);
    s->runTime->fullCheck_sample = (vector) arr_alloc(run.MAX_ITER * run.HORIZON + 1, double);

    /* The following numbers are recorded everytime a subproblem is solved */
    s->runTime->cutGenTime = (vector) arr_alloc(run.MAX_ITER * run.HORIZON * 2 + 1, double);
    s->runTime->stoStrucTime = (vector) arr_alloc(run.MAX_ITER * run.HORIZON * 2 + 1, double);
    s->runTime->subprobTime = (vector) arr_alloc(run.MAX_ITER * run.HORIZON * 2, double);
    s->runTime->argmaxTime = (vector) arr_alloc(run.MAX_ITER * run.HORIZON * 2 + 1, double);
    s->runTime->lpCntIndex = (intvec) arr_alloc(run.MAX_ITER * run.HORIZON * 2 + 1, int);


	return s;
}//END newSoln

void freeSolnType(solnType *soln) {

	if (soln) {
        if (soln->meanX) mem_free(soln->meanX);
		if (soln->candidX) mem_free(soln->candidX);
		if (soln->incumbX) mem_free(soln->incumbX);
		if (soln->piS) mem_free(soln->piS);
        if (soln->piSRatio) mem_free(soln->piSRatio);
		if (soln->piM) mem_free(soln->piM);
		if (soln->djM) mem_free(soln->djM);
        if (soln->runTime) freeRunTimeType(soln->runTime);
		mem_free(soln);
	}

}//END freeSolnType()

void freeRunTimeType(runTimeType *runTime) {

    if (runTime->iterIndex) mem_free(runTime->iterIndex);
    if (runTime->lpCntIndex) mem_free(runTime->lpCntIndex);
    if (runTime->evalOptTime) mem_free(runTime->evalOptTime);
    if (runTime->solveTime) mem_free(runTime->solveTime);
    if (runTime->warmUpTime) mem_free(runTime->warmUpTime);
    if (runTime->optTime) mem_free(runTime->optTime);
    if (runTime->fullCheck_lb) mem_free(runTime->fullCheck_lb);
    if (runTime->fullCheck_reform) mem_free(runTime->fullCheck_reform);
    if (runTime->fullCheck_sample) mem_free(runTime->fullCheck_sample);
    if (runTime->masterTime) mem_free(runTime->masterTime);
    if (runTime->cutGenTime) mem_free(runTime->cutGenTime);
    if (runTime->subprobTime) mem_free(runTime->subprobTime);
    if (runTime->argmaxTime) mem_free(runTime->argmaxTime);
    if (runTime->stoStrucTime) mem_free(runTime->stoStrucTime);

    mem_free(runTime);

}
