/*
 * subprob.c
 *
 *  Created on: Oct 16, 2015
 *      Author: gjharsha
 */
#include <time.h>
#include "subprob.h"

extern runType run;

/* This function will solve a new subproblem. This involves replacing the right-hand side of the subproblem with new values, based upon some
 * observation of omega, and some X vector of primal variables from the master problem.  Generally, the latest observation is used.  When
 * forming a normal cut, the candidate x should be used, while the incumbent x should be used for updating the incumbent cut. */
int solveSubprob(probType *prob, oneProblem *subprob, solnType *soln, vector Xvect, vector observ, double *mubBar, BOOL *subFeasFlag) {
    
#ifdef TRACE
    printf("\t~solveSubprob();\n");
#endif
	intvec 	ind = NULL;
	vector 	rhs = NULL;
	int 	k, stat1, stat2;
	int     status;
    clock_t tic, toc;

	if ( !(ind = (intvec) arr_alloc(prob->num->rows, int)) )
		errMsg("allocation", "solveSubporb", "ind", 0);

	for(k = 0; k < prob->sp->mar; k++)
		ind[k] = k;

	/* compute the right-hand side using current observation and first-stage solution */
	rhs = computeRHS(prob->num, prob->coord, prob->bBar, prob->Cbar, Xvect, observ);
	if ( rhs == NULL ) {
		errMsg("algorithm", "solveSubprob", "failed to compute subproblem right-hand side", 0);
		return 1;
	}

	/* change the right-hand side in the solver */
	stat1 = changeRHS(subprob->lp, prob->num->rows, ind, rhs+1);
	if ( stat1 ) {
		errMsg("solver", "solveSubprob", "failed to change the right-hand side in the solver",0);
		return 1;
	}

#if defined(MODIFY)
    static int ab;
    char subprobfname[30] = "/subprob/subprob    .lp";
    subprobfname[16] = '0' + ab / 1000 % 10;
    subprobfname[17] = '0' + ab / 100 % 10;
    subprobfname[18] = '0' + ab / 10 % 10;
    subprobfname[19] = '0' + ab / 1 % 10;
    ab++;
    writeProblem(subprob->lp, subprobfname);
#endif
    

    tic = clock();
	stat1 = solveProblem(subprob->lp, subprob->name, subprob->type, &stat2);
	if ( stat1 ) {
		if ( stat2 == STAT_INFEASIBLE ) {
			printf("Subproblem is infeasible: need to create feasibility cut.\n");
			(*subFeasFlag) = FALSE;
			return 0;
		}
		else {
			errMsg("algorithm", "solveSubprob", "failed to solve subproblem in solver", 0);
			return 1;
		}
	}
	else
		(*subFeasFlag) = TRUE;
    toc = clock();
    
    /* Recording Solving Time of a Subproblem */
    soln->runTime->subprobTime[soln->runTime->subprobCnt] = ((double) (toc - tic)) / CLOCKS_PER_SEC;
    soln->runTime->subprobCnt++;

#if defined(VALID)
	double objV;
	objV = getObjective(subprob->lp, PROB_LP);
	printf("\t :: Subproblem objective function value: %lf \n", objV);
#endif
    
    /* Collect dual solution of a subproblem */
	status = getDual(subprob->lp, soln->piS, prob->num->rows);
	if(status) {
		errMsg("algorithm", "solveSubprob", "failed to get the dual", 0);
		return 1;
	}

	/* mubBar used in stochastic updates */
	stat1 = computeMu(subprob->lp, prob->num->cols, mubBar);
	if ( stat1 ) {
		errMsg("algorithm", "solveSubprob", "failed to compute mubBar for subproblem", 0);
		return 1;
	}

	mem_free(rhs);
	mem_free(ind);

	return 0;
}//END solveSubprob

/**********************************************************************************************
 * This function computes the right hand side of the subproblem, based on a given X vector and a 
 *  given observation of omega.
 *
 * It is defined as:
 * 			rhs = R(omega) - T(omega) x X
 * and is calculated as:
 * 			rhs = (Rbar - Tbar x X) + (Romega - Tomega x X)
 *
 * where the "bar" denotes the fixed or mean value, and the "omega" denotes a random variation 
 * from this mean. The function allocates an array for the vector, which must be freed by the customer.  
 * Also, the zeroth position of this rhs vector is reserved, and the actual values begin at rhs[1].
 **********************************************************************************************/
vector computeRHS(numType *num, coordType *coord, sparseVector *bBar, sparseMatrix *Cbar, vector X, vector obs) {
	int cnt;
	vector rhs;
	sparseVector bomega;
	sparseMatrix Comega;

	bomega.cnt = num->rvbOmCnt;	bomega.col = coord->omegaRow; bomega.val=obs;

	Comega.cnt = num->rvCOmCnt; Comega.col = coord->omegaCol + num->rvbOmCnt;
	Comega.row = coord->omegaRow + num->rvbOmCnt; Comega.val = obs + num->rvbOmCnt;

	if (!(rhs =(vector) arr_alloc(num->rows+1, double)))
		errMsg("Allocation", "computeRhs", "rhs",0);

	/* Start with the values of b(omega) -- both fixed and varying */
    for (cnt = 1; cnt <= bBar->cnt; cnt++)
		rhs[bBar->col[cnt]] += bBar->val[cnt];
	for (cnt = 1; cnt <= bomega.cnt; cnt++)
		rhs[bomega.col[cnt]] += bomega.val[cnt];

	/* (cumulatively) subtract values of C(omega) x X -- both fixed and varying */
	rhs = MSparsexvSub(Cbar, X, rhs);
	rhs = MSparsexvSub(&Comega, X, rhs);
    
    for (cnt = 1; cnt <= bBar->cnt; cnt++)
        if ( fabs(rhs[bBar->col[cnt]]) - 0 < run.TOLERANCE )
            rhs[bBar->col[cnt]] = 0;

	return rhs;
}//END computeRHS()

/* since the basic structure of subproblem is not modified during the course of the algorithm, we just load it onto the solver */
oneProblem *newSubprob(oneProblem *subprob) {

	subprob->lp = setupProblem(subprob->name, subprob->type, subprob->mac, subprob->mar, subprob->objsen, subprob->objx, subprob->rhsx, subprob->senx,
			subprob->matbeg, subprob->matcnt, subprob->matind, subprob->matval, subprob->bdl, subprob->bdu, NULL, subprob->cname, subprob->rname, subprob->ctype);
	if ( subprob->lp == NULL ) {
		errMsg("Problem Setup", "newSubprob", "subprob",0);
		return NULL;
	}

	return subprob;
}//END newSubprob

/* Since the subproblem is not actually copied, we don't want to free it here. We just remove it from CPLEX so the next cell can load it back in. */
void freeSubprob(oneProblem *subprob){

	freeProblem(subprob->lp);

}//END freeSubprob
