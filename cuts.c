/*
 * cuts.c
 *
 */

#include <time.h>
#include "cuts.h"
#include "stocUpdt.h"
#include "evaluate.h"

extern runType run;

/*********************************************************************************************
 Once an observation of omega has been generated and the corresponding subproblem has been
 solved, a new cut may be formed for the master problem. This function will update all the data
 (concerning omega and the dual vector) and form a new cut based upon it.
 The new cut is added to the array of cuts, and the number of cuts is incremented.
 *********************************************************************************************/
int formOptCut(probType **p, cellType *c, solnType *s, BOOL newOmegaFlag, int omegaIdx) {

#ifdef TRACE
    printf("\t~formOptCut()\n");
#endif

	oneCut 	*cut;
	int		status = 0;
    clock_t tic, toc, tiic, tooc;

	// Preparing cut space
    cut = newCut(p[0]->num->cols, c->omega->cnt, c->k, c->t);

    // Updating stochastic structure
    tic = clock();
	stochasticUpdates(p[1], p[1]->num, s->candidX, p[1]->coord, p[1]->bBar, p[1]->Cbar, p[1]->aBar,
                      p[1]->Bbar, c->lambda, c->gamma, c->sigma, c->delta, c->auxDelta, c->trimDelta,
                      c->omega, c->auxOmega, c->trimOmega, newOmegaFlag, omegaIdx, c->maxIter, c->k,
                      s->piS, s->mubBar);
    toc = clock();
    s->runTime->stoStrucTime[s->runTime->argmaxCnt] = ((double) (toc - tic)) / CLOCKS_PER_SEC;

    // Generating cuts
    tic = clock();
    tiic = clock();
    status = tsSDCut(p[1]->num, p[1]->coord, c->lambda, c->gamma, c->sigma, c->delta, c->auxDelta,
                     c->trimDelta, c->omega, c->arima, cut, s->candidX, s->piSRatio, &s->dualStableFlag, c->k, c->Kt);
    if ( status ) {
        errMsg("algorithm", "formOptCut", "failed to create a new tsSD cut", 0);
        return 1;
    }
    tooc = clock();

    s->runTime->argmaxTime[s->runTime->argmaxCnt] = ((double) (tooc - tiic)) / CLOCKS_PER_SEC;
    s->runTime->lpCntIndex[s->runTime->argmaxCnt] = c->t;
    s->runTime->argmaxCnt++;
	status = addCut(p[0]->num, cut, s, c);
	if ( status < 0 ) {
		errMsg("algorithm", "formOptCut", "failed to add new cut to the problem", 0);
		return 1;
	}
    s->cCutIdx = status;
    toc = clock();

    /* Record time for generating a optimality cut */
    s->runTime->cutGenTime[s->runTime->cutGenCnt] = ((double) (toc - tic)) / CLOCKS_PER_SEC;
    s->runTime->cutGenCnt++;

#ifdef MODIFY
    char masterfname[30] = "/master/master   .lp";
    static int mnum = 0;
    masterfname[14] = '0' + mnum / 100 % 10;
    masterfname[15] = '0' + mnum / 10 % 10;
    masterfname[16] = '0' + mnum / 1 % 10;
    ++mnum;
    writeProblem(c->master->lp, masterfname);
#endif

	return 0;
}//END formNewCut()


/*********************************************************************************************
 We want to design a subroutine here that is similar to SDCut with gammaType information incor-
 porated. Also, this new rSDCut need to consider sub-hourly interval models using time series
 models for stochastic process representation. Besides the considerations in
 gammaType (all deterministic pi*ar, pi*ma) and noise in deltaType (immediate pi*noise), we also
 need to consider all intermedia terms which produced when generating newer sub-hourly observation
 using the ones just generated.
 Another important part is that we need to consider pre-process method here as well. If we are
 using stl/cla method, the seasonality and trend is taken care of in pibAR and pibMA and we
 don't have any more troubles. However, if pre-process inhereited EMP method, there will be an
 extra term in front of (pi*noise) term and this term is changing through time. We need to make
 sure that the cust we have incorporates all these parameters.
 *********************************************************************************************/
int tsSDCut(numType *num, coordType *coord, lambdaType *lambda, gammaType *gamma, sigmaType *sigma,
            deltaType *delta, deltaType *auxDelta, deltaType *trimDelta, omegaType *omega, tsType *arima,
            oneCut *cut, vector Xvect, vector piSRatio, BOOL *dualStableFlag, int numSamples, int periodStart) {

    vector 	piCbarX, beta;
    iType 	iStar;
    BOOL    piEvalFlag = FALSE;
    int 	c, cnt, obs;
    double  argmaxPre = -DBL_MAX, argmaxPost = -DBL_MAX;
    double  argmaxPreSum = 0.0, argmaxAllSum = 0.0;

#ifdef TRACE
    printf("\t\t~tsSDCut()\n");
#endif

    /* Check if we are ready to evaluate dual stability */
    if ( numSamples > run.PI_EVAL_START && !(numSamples % run.PI_CYCLE))
        piEvalFlag = TRUE;

    /* Need to store  Pi x Cbar x X independently of observation loop */
    if ( !(piCbarX= arr_alloc(sigma->cnt, double)) )
        errMsg("Allocation", "rSDCut", "pi_Tbar_x",0);
    if ( !(beta = (vector) arr_alloc(num->prevCols, double)) )
        errMsg("allocation", "rSDCut", "beta", 0);


    /* Calculate (Pi x Cbar) x X by mult. each VxT by X, one at a time */
    for (cnt = 0; cnt < sigma->cnt; cnt++) {
        piCbarX[sigma->linker[cnt]] = 0;
        for (c = 1; c <= num->cntCcols; c++)
            piCbarX[sigma->linker[cnt]] += sigma->vals[sigma->linker[cnt]].C[c] * Xvect[coord->colsC[c]];
    }

    /* Test for omega issues */
    for (obs = 0; obs < omega->cnt; obs++) {
        /* identify the maximal Pi for each observation */
        iStar = r_computeIstar(num, coord, arima, omega, lambda, gamma, sigma, delta,
                               auxDelta, trimDelta, Xvect, piCbarX,
                               &argmaxPre, &argmaxPost, obs, numSamples, periodStart);

        if ( iStar.delta < 0 || iStar.sigma < 0) {
            errMsg("algorithm", "SDcut", "failed to identify maximal Pi for an observation", 0);
            return 1;
        }

        cut->iStar[obs] = iStar.sigma;

#ifdef CUTS_SIGMA
        printf("\t :: iStar.sigma = %d || iStar.delta = %d \n", iStar.sigma, iStar.delta);
#endif

        //Newly added for rolling horizon problem
        if (arima != NULL) {
            cut->alpha += sigma->vals[iStar.sigma].b * omega->weight[obs];
            cut->alpha += delta->val[iStar.delta][obs].b * omega->weight[obs];
            cut->alpha += gamma->endo->pibARMA[iStar.delta] * omega->weight[obs];
            cut->alpha += auxDelta->val[iStar.delta][obs].b * omega->weight[obs];
            cut->alpha += trimDelta->val[iStar.delta][obs].b * omega->weight[obs];
        } else {
            cut->alpha += sigma->vals[iStar.sigma].b * omega->weight[obs];
            cut->alpha += delta->val[iStar.delta][obs].b * omega->weight[obs];
        }

        /* C SECTION: Average using these Pi's to calculate the cut itself. */
        for (c = 1; c <= num->cntCcols; c++)
            cut->beta[coord->colsC[c]] += sigma->vals[iStar.sigma].C[c] * omega->weight[obs];
        for (c = 1; c <= num->rvColCnt; c++)
            cut->beta[coord->rvCols[c]] += delta->val[iStar.delta][obs].C[c] * omega->weight[obs];

        // Collect Dual Impacts for Later Pi Ratio Calculation
        argmaxPreSum += max((argmaxPre - 0.0), 0) * omega->weight[obs];
        argmaxAllSum += max((max(argmaxPre, argmaxPost)-0.0), 0) * omega->weight[obs];

#ifdef CUTS
        printf("\t :: obs %d :: cut->alpha = %lf; \n", obs, cut->alpha);
#endif

    }

    if (piEvalFlag == TRUE)
        evalDualStability(piSRatio, argmaxPreSum, argmaxAllSum, dualStableFlag, numSamples);

    cut->alpha /= numSamples;

#ifdef CUTS
    printf("\t :: FINAL :: cut->alpha = %lf; \n", cut->alpha);
#endif

    for (c = 1; c <= num->prevCols; c++)
        cut->beta[c] /= numSamples;

    /*coefficient of eta coloumn*/
    cut->beta[0] =1.0;

    mem_free(piCbarX); mem_free(beta);

    return 0;
}//END rSDCut

/*********************************************************************************************
 This function will add a new cut to the master problem.  If This function adds a cut to the
 master problem, in the form of a constraint.
 When the objective function is: Min cx + eta, a cut is specified as:
 * 					eta >= Alpha + Beta x X.
 * In this function, it is entered in the master problem in this form:
 *                  Beta x X  +  Eta  >=  Alpha
 *
 For now, Eta does not have a coefficient.  All the eta coefficients (not just the one for
 this cut) will be updated in changeEtaCol(), before each master program is solved.
 The function returns 0 if the cut was successfully added; 1 otherwise.
 *********************************************************************************************/
int addCut(numType *num, oneCut *cut, solnType *soln, cellType *c) {

#ifdef TRACE
    printf("\t\t~addCut()\n");
#endif

	intvec indices; 	/* column indices of each beta coefficient */
	double rhs; 		/* rhs value in regularized QP method. */
	int cnt;
	int row;
	int status;


	/* Ensure that there is room to add another cut */
	if (c->cuts->cnt >= c->maxCuts) {
		status = reduceCuts(c, soln);
		if ( status ) {
			errMsg("algorithm", "addCut", "ran out of memory to add new cut, include reduceCuts",0);
			return -1;
		}
	}

	/*delete all feasibility cuts so that the new optimality cut can be added to the
                                                                end of the optimality cuts */
	for (cnt = c->fCutAdded->cnt - 1; cnt >= 0; cnt--) {
		row = c->fCutAdded->val[cnt]->rowNum;
		if ( !removeRow(c->master->lp, row, row) ) {
			errMsg("removeRow", "dropCut", "returned FALSE",0);
			return -1;
		}
	}

	/* Initialize an array to specify columns of each coefficient in beta.
       The one-norm of beta is temporarily used as the coefficient on
	   eta (it is assumed to be replaced in the next step in solveMaster()). */
	if (!(indices = arr_alloc(num->cols+1, int)))
		errMsg("Allocation", "addCut", "coefCol",0);
	for (cnt = 1; cnt <= num->cols; cnt++)
		indices[cnt] = cnt-1;
	indices[0] = num->cols;

	/* Add the cut (it's a ">=" constraint) to the master, with coefficients as specified in beta,
       and right hand side as specified by alpha.
	  In the regularized QP method, we need to shift the rhs of the cut from 'x' to 'd' each time
      we add a cut. (We do not need to worry about it when dropping a cut.) That is, in regularized
      QP method, the rhs will become
                                    alpha - beta * incumb_x
	  instead of alpha as in the LP method. */
	if (run.MASTER_TYPE == PROB_LP && c->k > 1)
		rhs = cut->alpha;
	else
		rhs = cut->alpha - vXv(cut->beta, soln->incumbX, NULL, num->cols);

	cut->alphaIncumb = rhs;

#if defined(OFF)
    printf("\t\t :: iteration %d; cut->alphaIncumb = %lf, incumbent ? %s \n", c->k, rhs, cut->isIncumb ? "TRUE" : "FALSE");
#endif

	/* add the row in the solver */
	status = addRow(c->master->lp, num->cols+1, rhs, GE, 0, indices, cut->beta);
	if (status){
		errMsg("solver", "addCut", "failed to add new row to problem in solver",0);
		return -1;
	}
	cut->rowNum = num->rows + c->cuts->cnt;

	/* Since a new optimality cut is added, the row number of every feasibility cut will increase by 1*/
	for (cnt = 0; cnt < c->fCutAdded->cnt; cnt++)
		++c->fCutAdded->val[cnt]->rowNum;

	/* add all the removed feasibility cuts back to the problem after optimality cuts*/
	for (cnt = 0; cnt < c->fCutAdded->cnt; cnt++) {
		if (run.MASTER_TYPE == PROB_LP)
			rhs = c->fCutAdded->val[cnt]->alpha;
		else
			/* Using the same tolerance as the rest, but the original SD code uses a different
                                                                    tolerance level to check feasibility*/
			rhs = c->fCutAdded->val[cnt]->alpha +vXv(c->fCutAdded->val[cnt]->beta,
                                                       soln->incumbX , NULL, num->cols);

		/* if the lower bound is non-zero, then the cuts should be shifted appropriately */
		if ( c->lbType !=TRIVIAL )
			rhs += ((double) c->k / (double) c->cuts->val[cnt]->cutObs - 1)* c->lb;

		if (!addRow(c->master->lp, num->cols+1, rhs, GE, 1, indices, c->fCutAdded->val[cnt]->beta)) {
			errMsg("solver", "addCutToMaster", "failed to add new row to problem in solver", 0);
			return -1;
		}
	}

	/* add the newly created cut to cutsTYpe structure */
	c->cuts->val[c->cuts->cnt] = cut;

	mem_free(indices);

	return c->cuts->cnt++;
}//END addCut


/*********************************************************************************************
 Description goes here...
 *********************************************************************************************/
int formFeasCut(probType **prob, cellType *cell, solnType *soln, BOOL *newOmegaFlag, int omegaIdx) {
	int start, end;
	int idx;

	/* update the stochastic elements that are effected are effected */
	stochasticUpdates(prob[1], prob[1]->num, soln->candidX, prob[1]->coord, prob[1]->bBar, prob[1]->Cbar, prob[1]->aBar,prob[1]->Bbar,
                      cell->lambda, cell->gamma, cell->sigma, cell->delta, cell->auxDelta, cell->trimDelta, cell->omega, cell->auxOmega,
                      cell->trimOmega, (*newOmegaFlag), omegaIdx, cell->maxIter, cell->k, soln->piS, soln->mubBar);
	/* since stochastic updates for new omega have been completed, set it equal to zero to avoid recalculating */
	(*newOmegaFlag) = FALSE;

	/* add new feasibility cuts to the cut pool */
	updtFeasCutPool(prob[1], cell);

	/* identify, in the feasibility cut pool, cuts that are violated by the input solution xk */
	start = cell->fCutAdded->cnt;
	checkFeasCutPool(cell->fCutPool, cell->fCutAdded, prob[0]->num->cols, soln->incumbX, soln->candidX, &soln->infeasIncumb);
	end = cell->fCutAdded->cnt;

	/* add feasibility cuts to master problem */
	if (end > start) {
		for (idx = start; idx < end; idx++) {
			addfCut(cell->master->type,cell->master->lp, prob[0]->num, cell->cuts->cnt, soln->incumbX, cell->fCutAdded->val[idx], idx);
			writeProblem(cell->master->lp, "masterCut.lp");
		}
		/* make room for dual solutions for new feasibility cuts added */
		soln->piM = (vector) mem_realloc(soln->piM, prob[0]->num->rows+cell->cuts->cnt+cell->fCutAdded->cnt);
	}

	return 0;
}//END formFeasCut()

/*********************************************************************************************
 This function adds new feasibility cuts. It first adds feasibility cuts from old pi's
 associated with the new omega generated. Cuts from a new dual extreme ray(new pi) and all omegas
 generated so far are added to the feasible_cuts_pool structure afterwards.
 *********************************************************************************************/
int updtFeasCutPool(probType *prob, cellType *cell) {
	vector 	beta;
	double	alpha;
	int		m, n, c, deltaIdx, cutCnt;

	if ( !(beta = (vector) arr_alloc(prob->num->prevCols+1, double)) )
		errMsg("allocation", "updtFeasCutPool", "beta", 0);

	/* keep track of added cuts to feasibility-cut pool */
	cutCnt = cell->fCutPool->cnt;

	for ( n = cell->fCol; n < cell->omega->cnt; n++ ) {
		for ( m = cell->fRow; m < cell->sigma->cnt; m++ ) {
			deltaIdx = cell->sigma->lambdaIdx[m];
			alpha = 0.0;
			for (c = 0; c <= prob->num->prevCols; c++)
				beta[c] = 0.0;
			alpha = cell->sigma->vals[m].b + cell->delta->val[deltaIdx][n].b;
			for (c = 1; c <= prob->num->cntCcols; c++)
				beta[prob->coord->colsC[c]] += cell->sigma->vals[m].C[c];
			for (c = 1; c <= prob->num->rvColCnt; c++)
				beta[prob->coord->rvCols[c]] += cell->delta->val[deltaIdx][n].C[c];
			add2CutPool(cell, alpha, beta, prob->num->prevCols, cell->omega->cnt);
		}
	}
	cell->fCol = cell->omega->cnt;
	cell->fRow = cell->sigma->cnt;

	mem_free(beta);

	return cutCnt;
}//END updtFeasCutPool()

/*********************************************************************************************
 This function add a new feasibility cut to the cut pool using alpha and beta provided.
 *********************************************************************************************/
int add2CutPool(cellType *cell, double alpha, vector beta, int betaLen, int numOmega) {
	oneCut 	*cut;
	int 	cnt;

	for (cnt = 0; cnt < cell->fCutPool->cnt; cnt++) {
		if (DBL_ABS(alpha - cell->fCutPool->val[cnt]->alpha) < run.TOLERANCE) {
			if (equalVector(beta, cell->fCutPool->val[cnt]->beta, betaLen, run.TOLERANCE)) {
				/* return 0 to indicate that no cut was added to the pool */
				return 1;
			}
		}
	}

	if ( !(cut = (oneCut *) mem_malloc (sizeof(oneCut))))
		errMsg("allocation", "add2CutPool", "cut", 0);
	cut->cutObs = cell->k;
	cut->omegaCnt = numOmega;
	cut->slackCnt = 0;
	cut->isIncumb = FALSE;

	if ( !(cut->iStar = (intvec) arr_alloc(numOmega, int)) )
		errMsg("allocation", "add2CutPool", "istar", 0);
	if ( !(cut->beta = arr_alloc(betaLen+1, double)))
		errMsg("allocation", "add2CutPool", "beta", 0);

	cut->alpha = alpha;
	for (cnt = 0; cnt <= betaLen; cnt++)
		cut->beta[cnt] = beta[cnt];

	cell->fCutPool->val[cell->fCutPool->cnt++] = cut;

	return 0;
}//END add2CutPool()

/*********************************************************************************************
 The function identifies cuts from the feasibility cut pool which are voilated by the candidate
 solution, and mark them to be added to master problem.
 *********************************************************************************************/
int checkFeasCutPool(cutType *cutPool, cutType *cutsAdded, int betaLen, vector incumbX, vector candidX, BOOL *infeasIncumb) {
	double 	betaX, alpha;
	int 	idx, c;
	BOOL 	duplicCut;

	for (idx = 0; idx < cutPool->cnt; idx++) {
		duplicCut = FALSE;
		alpha = cutPool->val[idx]->alpha;
		for (c = 0; c < cutsAdded->cnt; c++) {
			if (DBL_ABS(alpha - cutsAdded->val[c]->alpha) < run.TOLERANCE) {
				if (equalVector(cutPool->val[idx]->beta, cutsAdded->val[c]->beta, betaLen, run.TOLERANCE)) {
					duplicCut = TRUE;
					break;
				}
			}
		}

		/* Add those cuts in cut pool that will be violated by incumbent solution */
		betaX = vXv(cutPool->val[idx]->beta, incumbX, NULL, betaLen);
		if (betaX < alpha) {
			(*infeasIncumb) = TRUE;
			if (duplicCut == TRUE) {
				printf("Incumbent violates one old cut from feasible cut pool (this cut also exists in feasCutsAdded)\n");
				continue;
			}
			else
				printf( "Incumbent violates one new cut from feasible cut pool (this cut is not in feasCutsAdded but will be added)\n");
			cutsAdded->val[cutsAdded->cnt++] = cutPool->val[idx];

			printf("Cut added to master due to Incumbent violation\n");
		}
		else {
			/* Check if the cut will be violated by the candidate solution*/
			if (duplicCut == TRUE)
				continue;
			betaX = vXv(cutPool->val[idx]->beta, candidX, NULL, betaLen);

			if (betaX < alpha) {
				printf("Candidate violates one cut from feasible cut pool (this cut is not in feasCutsAdded but will be added)\n");
				cutsAdded->val[cutsAdded->cnt++] = cutPool->val[idx];

				printf("Cut added to master due to candidate violation\n");
			}
		}

	}

	return 0;
}//END checkFeasCutPool()

/*********************************************************************************************
 This function will add new feasibility cut to the master problem. Unlike addCut(), we do not
 rearrange the cuts while adding.
 *********************************************************************************************/
int addfCut(int type,LPptr lp, numType *num, int optCuts, vector incumbX, oneCut *cut, int idx) {
	intvec 	indices;
	double	rhs;
	int		cnt, status;

	/*
     Initialize an array to specify columns of each coefficient in beta. The one-norm of beta
     is temporarily used as the coefficient on eta (it is assumed to be replaced in the next
     step in solveMaster()).
    */
	if (!(indices = (intvec) arr_alloc(num->cols+1, int)))
		errMsg("Allocation", "addCut", "coefCol",0);
	for (cnt = 0; cnt < num->cols; cnt++)
		indices[cnt + 1] = cnt;
	indices[0] = num->cols;

	/*
     Add the cut (it's a ">=" constraint) to the master, with coefficients as specified in
     beta, and right hand side as specified by alpha.
     In the regularized QP method, we need to shift the rhs of the cut from 'x' to 'd' each
     time we add a cut. (We do not need to worry about it when dropping a cut.) That is, in
     regularized QP method, the rhs will become
                                        alpha - beta * incumb_x
	 instead of alpha as in the LP method.
     */
	if (type == PROB_LP)
		rhs = cut->alpha;
	else
		rhs = cut->alpha - vXv(cut->beta, incumbX, indices, num->cols);

	/* add the row in the solver */
	status = addRow(lp, num->cols+1, rhs, GE, 0, indices, cut->beta);
	if (status){
		errMsg("solver", "addCut", "failed to add new row to problem in solver",0);
		return 1;
	}
	cut->rowNum = num->cols + optCuts + idx;

	mem_free(indices);

	return 0;
}//END addfCut()


/*********************************************************************************************
 As iterations proceed, the cut which was formed based on the incumbent X vector is periodically
 re-evaluated. If tau iterations have passed since the last update, then the cut is reformed.
 Or, if the value of the incumbent cut is less than the value of f_k at the incumbent X,
 the incumbent is reformed. This function will solve an additional subproblem (based upon the
 incumbent X and the most recent observation of omega), add the dual solution to the data
 structures, and re-form the incumbent cut (thus replacing it in the array of cuts).
 Note that the updates do NOT treat the omega as a new observation!

 This is a newer version of formIncumbCut() that has not been tested.
 Different from formIncumbCut(), we need the fullObs vector that stores the real right hand side.
 Different from fromIncumbCut(), the SDCut part is diversed depend on if we used time series model
 or not.
 *********************************************************************************************/
int rFormIncumbCut(probType **p, cellType *c, solnType *s, int omegIdx, vector fullObs, vector auxObs, BOOL newOmegaFlag) {

#ifdef TRACE
    printf("\t~r_formIncumbCut();\n"); printf("\t");
#endif
    oneCut 	*cut;
    int		n;
    clock_t tic, toc, tiic, tooc;
    int     status;

    /* solve the subproblem with current observation and incumbent solution */
    solveSubprob(p[1], c->subprob, s, s->incumbX, fullObs, &s->mubBar, &s->subFeasFlag);

    /* Record time spent on argmax procedures without counting the time for solving the subproblem LP. zl, 06/30/04. */
    if (s->subFeasFlag == TRUE) {
        /* Drop the previous incumbent cut */
        if (c->cuts->cnt >= c->maxCuts) {
            status = dropCut(c->master->lp, c->cuts, s, s->iCutIdx);
            if ( status ) {
                errMsg("algorithm", "formIncumbCut", "failed to drop old incumbent cut", 0);
                return 1;
            }

            /* decrease row number for each added feasibility cut */
            for ( n = 0; n < c->fCutAdded->cnt; n++ )
                --c->fCutAdded->val[n]->rowNum;
        }

    }

    if (s->subFeasFlag == FALSE) {
        /* incumbent solution is infeasible to subproblem. Replace incumbent with candidate solution */
        status = replaceIncumbent(p, c, s, s->candidEst);
        if ( status ) {
            errMsg("algorithm", "formIncumbCut", "failed to replace incumbent with candidate after incumbent infeasibility", 0);
            return 1;
        }
    } else {
        tic = clock();
        cut = newCut(p[0]->num->cols, c->omega->cnt, c->k, c->t);
        cut->isIncumb = TRUE;

        stochasticUpdates(p[1], p[1]->num, s->incumbX, p[1]->coord, p[1]->bBar, p[1]->Cbar, p[1]->aBar, p[1]->Bbar,
                          c->lambda, c->gamma, c->sigma, c->delta, c->auxDelta, c->trimDelta, c->omega, c->auxOmega,
                          c->trimOmega, newOmegaFlag, omegIdx, c->maxIter, c->k, s->piS, s->mubBar);

        tiic = clock();
        status = tsSDCut(p[1]->num, p[1]->coord, c->lambda, c->gamma, c->sigma, c->delta, c->auxDelta, c->trimDelta,
                         c->omega, c->arima, cut, s->incumbX, s->piSRatio, &s->dualStableFlag, c->k, c->Kt);
        if ( status ) {
            errMsg("algorithm", "formOptCut", "failed to create a new tsSD cut", 0);
            return 1;
        }
        tooc = clock();

        s->runTime->argmaxTime[s->runTime->argmaxCnt] = ((double) (tooc - tiic)) / CLOCKS_PER_SEC;
        s->runTime->lpCntIndex[s->runTime->argmaxCnt] = c->t;
        s->runTime->argmaxCnt++;

        status = addCut(p[0]->num, cut, s, c);
        if (status < 0){
            errMsg("Algorithm", "formIncumbCut", "Fail to add incumbent cut",0);
            return 1;
        }

        s->iCutIdx = status;
        s->iCutUpdt = c->k;
        toc = clock();
        s->runTime->cutGenTime[s->runTime->cutGenCnt] = ((double) (toc - tic)) / CLOCKS_PER_SEC;
        s->runTime->cutGenCnt++;

    }

    return 0;
}//END formIncumbentCut

/*This function loops through all the dual vectors found so far and returns the index of the one which satisfies the expression:
 * 				argmax { Pi x (R - T x X) | all Pi }
 * where X, R, and T are given.  It is calculated in this form:
 * 				Pi x bBar + Pi x bomega + (Pi x Cbar) x X + (Pi x Comega) x X.
 * Since the Pi's are stored in two different structures (sigma and delta), the index to the maximizing Pi is actually a structure
 * containing two indices.  (While both indices point to pieces of the dual vectors, sigma and delta may not be in sync with one
 * another due to elimination of non-distinct or redundant vectors. */
iType computeIstar(numType *num, coordType *coord, sigmaType *sigma, deltaType *delta, vector Xvect,
                   vector PiCbarX, double *argmaxPre, double *argmaxPost, int obs, int numSamples) {

    /*
     * This subroutine is re-designed to also output the argmax in 90% | 10%,
     * which shouldn't affect how the iStar is calculated
     */

	double 	arg, argmax;
	int 	sigPi, delPi, slice;
	iType 	ans;

    /* Calculating the divide slice position of 90% | 10% */
    slice = numSamples / 10 + 1;
    slice = numSamples - slice;

    ans.delta = -1; ans.sigma = -1; ans.height = -DBL_MAX;
    argmax = -DBL_MAX; *argmaxPre =-DBL_MAX; *argmaxPost = -DBL_MAX;

	for (sigPi = 0; sigPi < sigma->cnt; sigPi++) {
		/* Find the row in delta corresponding to this row in sigma */
		delPi = sigma->lambdaIdx[sigPi];

		/* Start with (Pi x bBar) + (Pi x bomega) + (Pi x Cbar) x X */
		arg = sigma->vals[sigPi].b + delta->val[delPi][obs].b - PiCbarX[sigPi];

		/* Subtract (Pi x Comega) x X. Multiply only non-zero VxT values */
		arg -= vXv(delta->val[delPi][obs].C, Xvect, coord->rvCols, num->rvColCnt);

		if (arg > argmax) {
            argmax = arg;
			ans.sigma = sigPi;
            ans.delta = delPi;
            ans.height = argmax;
		}

        /* Fetching the best dual base on first 90% of iterations */
        if (arg > *argmaxPre && sigma->ck[sigPi] <= slice) {
            *argmaxPre = arg;
        }

        /* Fetching the best dual base on most recent 10% iterations */
        if (arg > *argmaxPost && sigma->ck[sigPi] > slice)
            *argmaxPost = arg;
	}
	return ans;
}//END computeIstar

/*********************************************************************************************
 TEMP workspace for new implementation for rolling sp problem.

 ~~~> Comments from [SITE]
 This subroutine is an modified version of computeIstar(). It is applied in the rolling horizon
 problem. In the rolling horizon problem, the cuts are determined withing some extra terms. Hence
 in here, we need to consider these extra contributions when choosing the best dual multiplier.

 ~~~> Function description
 This function loops through all the dual vectors found so far and returns the index of the one
 which satisfies the expression:
    argmax { Pi * ( b(omega) - TxX) | for all pi}
        where b(.) is an pre-process function of data in the time series model
          and omega is the time series model.

 If no pre-process is conducte before fitting time series model and with the assumption that
 omega ~ ARMA(p,q), we have a simpler form here,

 argmax {[Pi_r,Pi'] * ([gamma.AR + gamma.MA + auxObs + epsilon; bBar] - T x X)}
        => argmax { Pi_r x gamma.AR + Pi_r x gamma.MA + Pi_r x auxObs + Pi_r x epsilon(delta)
                    + Pi' x bBar - Pi x (T x X) }
            where Pi = [Pi_r, Pi']

 We are trying to mask all pre-process before in other places. The numbers brought in should be
 in the above form with no extra pre-process. The warm-up is conducted differently with extra
 coefficients masked in these numbers. However, we will write out more instances when different
 pre-process is conduct before hand.

 In the above operation,
     -Pi x (gamma.AR) is stored in gamma->pibAR, it is periodic deterministic AR contribution
     -Pi x (gamma.MA) is stored in gamma->pibMA, it is periodic determinisitc MA contribution
     -Pi x (auxObs) only shows up(not zero) in sub-hourly model, they are simulations that are
        based on previous simulations, note that it is a dynamic number that is a function
        of (p,q,delta,phi,theta).
        Calculation is shown here...

     -Pi' x bBar is stored in sigma type, note that here Pi' represents dual multipliers for
        for all deterministic rows in subproblem
     -Pi x epsilon(delta) is stored in deltaType and stdev is stored in tsType
     -We will inheret the say how pi x (T x X) is calculated here with (Pi x Comega) x X since we
        assume that C is deterministic. Will come back to this part later

 In the following comments, we listed a 2 possible pre-processes method. The numbers comes into
 this subroutine should have been finished-reverse-pro-processed and what are carried down here to
 make sure we are calculating the right number.
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 With the assumption that omega ~ ARMA(p,q) and we used EMP preprocess method, we have,

    argmax { [Pi_r,Pi'] * ( [stdev_t * gammaAR + mu_t + stdev * gammaMA + Inter +
                stdev_t * epsilon, bBar]  - TxX) }
        => argmax { Pi_r x (stdev_t*gammaAR + mu_t) + Pi_r x (stdev_t*gammaMA)
                    + (Pi_r x Inter) + Pi' x bBar + stdev * (Pi x epsilon(delta)) - Pi x (T x X) }
        where Pi = [Pi_r, Pi']

 In the above operation,
    -Pi x (stdev_t*gamma.AR+mu_t) is stored in gamma->pibAR
    -Pi x (stdev_t*gamma.MA) is stored in gamma->pibMA
    -Pi x (stdev_t*auxObs) only shows up in sub-hourly model
        Calculation is shown here...

    -Pi' x bBar is stored in sigma type
    -Pi x epsilon(delta) is stored in deltaType and stdev is stored in tsType !!! Issue here !!!
    -We will inheret the say how pi x (T x X) is calculated here with (Pi x Comega) x X since we
            assume that C is deterministic. Will come back to this part later

  In a different circumstance, assume that we use STL/CLA method for data pre-process and
    omega ~ ARMA(p,q), we have,

    argmax {[Pi_r,Pi'] * ([season + trend + gamma.AR + gamma.MA + auxObs + epsilon; bBar] - T x X)}
        => argmax { Pi_r x (season + trend + gamma.AR) + Pi_r x gamma.MA + Pi_r x auxObs
                    + Pi_r x epsilon(delta) + Pi' x bBar - Pi x (T x X) }
        where Pi = [Pi_r, Pi']

 In the above operation,
     -Pi x (season + trend + gamma.AR) is stored in gamma->pibAR
     -Pi x (gamma.MA) is stored in gamma->pibMA
     -Pi x (auxObs) only shows up in sub-hourly model
        Calculation is shown here...

     -Pi' x bBar is stored in sigma type
     -Pi x epsilon(delta) is stored in deltaType and stdev is stored in tsType
     -We will inheret the say how pi x (T x X) is calculated here with (Pi x Comega) x X since we
     assume that C is deterministic. Will come back to this part later
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 The only thing that we need to worry scalar in front of delta when using EMP method for
 sub-hourly models.
 *********************************************************************************************/
iType r_computeIstar(numType *num, coordType *coord, tsType *arima, omegaType *omega, lambdaType *lambda,
                     gammaType *gamma, sigmaType *sigma, deltaType *delta, deltaType *auxDelta,
                     deltaType *trimDelta, vector Xvect, vector PiCbarX, double *argmaxPre,
                     double *argmaxPost, int obs, int numSamples, int periodStart) {

    double 	arg, argmax;
    int 	linkerSig, sigPi, delPi, slice;
    iType 	ans;

    /* Calculating the divide slice position of 90% | 10% */
    slice = (numSamples - periodStart) / 10 + 1;
    //slice = numSamples / 10 + 1;
    slice = numSamples - slice;

    ans.delta = -1;
    ans.sigma = -1;
    ans.height = -DBL_MAX;
    argmax = -DBL_MAX;
    *argmaxPre = -DBL_MAX;
    *argmaxPost = -DBL_MAX;

    // Replacement ::
    for (linkerSig = 0; linkerSig < sigma->aCnt; linkerSig++) {

        delPi = sigma->lambdaIdx[sigma->linker[linkerSig]];
        sigPi = sigma->linker[linkerSig];

        /* Full Scheme */
        if (arima != NULL) {
            arg = sigma->vals[sigPi].b;
            arg -= PiCbarX[sigPi];
            arg += delta->val[delPi][obs].b;
            arg += gamma->endo->pibARMA[delPi];
            arg += auxDelta->val[delPi][obs].b;
            arg += trimDelta->val[delPi][obs].b;
        } else {
            arg = sigma->vals[sigPi].b;
            arg -= PiCbarX[sigPi];
            arg += delta->val[delPi][obs].b;
        }

        /* Subtract (Pi x Comega) x X. Multiply only non-zero VxT values */
        arg -= vXv(delta->val[delPi][obs].C, Xvect, coord->rvCols, num->rvColCnt);

        /* Finding the best dual among all possible information */
        if (arg > argmax) {
            argmax = arg;
            ans.sigma = sigPi;
            ans.delta = delPi;
            ans.height = argmax;
        }

        /* Fetching the best dual base on first 90% of iterations */
        if (arg > *argmaxPre && sigma->ck[sigPi] <= slice)
            *argmaxPre = arg;
        /* Fetching the best dual base on most recent 10% iterations */
        if (arg > *argmaxPost && sigma->ck[sigPi] > slice)
            *argmaxPost = arg;
    }

//    for (sigPi = 0; sigPi < sigma->cnt; sigPi++) {
//
//        /* Find the row in delta corresponding to this row in sigma */
//        delPi = sigma->lambdaIdx[sigPi];
//
//        if (lambda->alive[delPi] != DROPPED) {
//            /* Full Scheme */
//            if (arima != NULL) {
//                arg = sigma->vals[sigPi].b;
//                arg -= PiCbarX[sigPi];
//                arg += delta->val[delPi][obs].b;
//                arg += gamma->endo->pibARMA[delPi];
//                arg += auxDelta->val[delPi][obs].b;
//                arg += trimDelta->val[delPi][obs].b;
//            } else {
//                arg = sigma->vals[sigPi].b;
//                arg -= PiCbarX[sigPi];
//                arg += delta->val[delPi][obs].b;
//            }
//
//            /* Subtract (Pi x Comega) x X. Multiply only non-zero VxT values */
//            arg -= vXv(delta->val[delPi][obs].C, Xvect, coord->rvCols, num->rvColCnt);
//
//            /* Finding the best dual among all possible information */
//            if (arg > argmax) {
//                argmax = arg;
//                ans.sigma = sigPi;
//                ans.delta = delPi;
//                ans.height = argmax;
//            }
//
//            /* Fetching the best dual base on first 90% of iterations */
//            if (arg > *argmaxPre && sigma->ck[sigPi] <= slice)
//                *argmaxPre = arg;
//            /* Fetching the best dual base on most recent 10% iterations */
//            if (arg > *argmaxPost && sigma->ck[sigPi] > slice)
//                *argmaxPost = arg;
//        }
//    }

    return ans;
}//END computeIstar

/* This function allocates memory for the arrays inside a single cut, and initializes its values accordingly.  The cut structure
 * itself is assumed to be already allocated.  Note, each beta vector contains room for its one-norm, thought it just gets filled
 * with zero anyway. */
oneCut *newCut(int numX, int numIstar, int numSamples, int numPeriod) {

#ifdef TRACE
    printf("\t\t~newCut()\n");
#endif
	oneCut *cut;

	cut = (oneCut *) mem_malloc (sizeof(oneCut));
    cut->cutPeriod = numPeriod;
	cut->cutObs   = numSamples;
	cut->omegaCnt = numIstar;
	cut->slackCnt = 0;
	cut->isIncumb = FALSE; 			/*new cut is by default not an incumbent*/
	cut->alphaIncumb = 0;

	if (!(cut->iStar = arr_alloc(numIstar, int)))
		errMsg("allocation", "newCut", "iStar", 0);
	if (!(cut->beta = arr_alloc(numX+1, double)))
		errMsg("allocation", "new_cut", "beta", 0);

	cut->alpha = 0.0;

	return cut;
}//END newCut

/* This function allocates memory for a new cut structure.  This entails the structure itself, and the _val_ array of oneCut pointers
 * inside the structure.  The actual oneCut structures are allocated according to the numBeta parameter, via calls to newCut(). */
cutType *newCuts(int maxCuts) {
	cutType *cuts;

	if (!(cuts = (cutType *) mem_malloc (sizeof(cutType))))
		errMsg("allocation", "newCuts", "cuts",0);
	if (!(cuts->val = (oneCut **) arr_alloc (maxCuts, oneCut)))
		errMsg("allocation", "newCuts", "oneCuts",0);
	cuts->cnt = 0;

	return cuts;
}//END newCuts


/* This function will remove the oldest cut whose corresponding dual variable is zero (thus, a cut which was slack in last solution). */
int reduceCuts(cellType *cell, solnType *soln) {
	int status, minObs, oldestCut, idx;

	minObs 	  = cell->k;
	oldestCut = cell->cuts->cnt;

	/* identify the oldest loose cut */
    for (idx = 0; idx < cell->cuts->cnt; idx++) {
		if ( idx == soln->iCutIdx )
			/* avoid dropping incumbent cut*/
			continue;

		if (cell->cuts->val[idx]->cutObs < minObs && DBL_ABS(soln->piM[cell->cuts->val[idx]->rowNum + 1]) <= run.TOLERANCE ) {
			minObs = cell->cuts->val[idx]->cutObs;
			oldestCut = idx;
		}
	}

	/* if the oldest loose cut is the most recently added cut */
	if ( oldestCut == cell->cuts->cnt ) {
		errMsg("algorithm", "reduceCuts", "failed to identify any cuts to drop", 0);
		return 1;
	}

	/* drop the selected cut and swap the last cut into its place */
	status = dropCut(cell->master->lp, cell->cuts, soln, oldestCut);
	if ( status ) {
		errMsg("algorithm", "reduceCuts", "failed to drop a cut", 0);
		return 1;
	}

	/* decrease row number for each added feasibility cut */
	for ( idx = 0; idx < cell->fCutAdded->cnt; idx++ )
		--cell->fCutAdded->val[idx]->rowNum;

	return 0;
}//END reduceCuts()

/* This function removes a cut from both the cutType structure and the master problem constraint matrix.  In the cuts->val array, the last
 * cut is swapped into the place of the exiting cut.  In the constraint matrix, the row is deleted, and the row numbers of all constraints
 * below it are decremented. */
int dropCut(LPptr lp, cutType *cuts, solnType *soln, int cutIdx) {

#ifdef TRACE
    printf("\t\t~dropCut()\n");
#endif

	int idx, status, deletedRow;

	deletedRow = cuts->val[cutIdx]->rowNum;
	/* Get rid of the indexed cut */
	status = removeRow(lp, deletedRow, deletedRow);
	if ( status ) {
		errMsg("solver", "dropCut", "failed to remove a row from master problem", 0);
		return 1;
	}
	freeCut(cuts->val[cutIdx]);

	/* Update the surviving cuts */
	cuts->val[cutIdx] = cuts->val[--cuts->cnt];
	for (idx = 0; idx < cuts->cnt; idx++)
		if (cuts->val[idx]->rowNum > deletedRow)
			--cuts->val[idx]->rowNum;

    /* if the swapped cut happens to be the incumbent cut, then update its index in solnType */
    if (soln->cCutIdx == cuts->cnt)
        soln->cCutIdx = cutIdx;

    if ( soln->iCutIdx == cuts->cnt )
        soln->iCutIdx = cutIdx;


	return 0;
}//END dropCut()

void freeCut(oneCut *cut) {

	if (cut) {
		if (cut->iStar)
			mem_free(cut->iStar);
		if (cut->beta)
			mem_free(cut->beta);
		mem_free(cut);
	}
}

void freeCuts(cutType *cuts) {
	int cnt;

	for (cnt = 0; cnt < cuts->cnt; cnt++)
		freeCut(cuts->val[cnt]);
	mem_free(cuts->val);
	mem_free(cuts);
}//END freeCuts
