/*
 * master.c

 *
 *  Created on: Oct 16, 2015
 *      Author: Harsha Gangammanavar
 */

#include <time.h>
#include "master.h"

extern runType run;
extern ENVptr env;

/* This function is the regularized QP version of master problem. The master problem is solved after the newest cut is added to master problem,
the incumbent cut is updated if necessary. Here the coefficients on all the cuts are updated, and finally master problem is solved. */
int solveQPMaster(numType *num, sparseVector *dBar, cellType *cell, solnType *soln) {
    
#ifdef TRACE
    printf("\t~solveQPMaster()\n");
#endif
    
	double 	ht, Sm = 0.0, d2 = 0.0; /* height at the candidate solution. */
	int 	status=0, stat1=0, i;
    clock_t tic, toc;

	/* Update eta coefficient on all cuts, based on cutObs */
	status = changeEtaCol(cell->master->lp, num->rows, num->cols, cell->k, cell->cuts, cell->lb);
	if ( status ) {
		errMsg("algorithm", "solveMaster", "failed to change the eta column coefficients", 0);
		return 1;
	}

	if ( cell->lbType == NONTRIVIAL ) {
#ifdef CALC
        printf("\t :: Cell lower bound is NONTRIVIAL. Updating cuts right hand side now...\n");
#endif
		/* update the right-hand side of cuts to reflect the non-trivial lower bound */
		status = updateRHS(cell);
		if ( status ) {
			errMsg("algorithm", "solveQPMaster", "failed to update right-hand side with lower bound information", 0);
			return 1;
		}
	}
    
	if ( soln->incumbChg )
		soln->incumbChg = FALSE;

#ifdef MODIFY
    char incumbMasterfname[30] = "/incumbMaster/imaster    .lp";
    static int imnum = 0;
    incumbMasterfname[21] = '0' + imnum / 1000 % 10;
    incumbMasterfname[22] = '0' + imnum / 100 % 10;
    incumbMasterfname[23] = '0' + imnum / 10 % 10;
    incumbMasterfname[24] = '0' + imnum / 1 % 10;
    ++imnum;
	writeProblem(cell->master->lp, incumbMasterfname);
#endif

	    
    /* solve the master problem */
	tic = clock();
	status = solveProblem(cell->master->lp, cell->master->name, cell->master->type, &stat1);
    if ( stat1 == 3 ) {
        //Reare special case met, need to re-solve it using primal simplex
        stat1 = CPXprimopt(env, cell->master->lp);
        printf("Barrier algorithm indicates infeasibility. Re-solving through primal simplex method\n");
        status = CPXgetstat(env, cell->master->lp);
        if ( status != 1 ) {
            errMsg("algorithm", "solveMaster", "failed to solve the master problem", 0);
            return 1;
        }
    }
	toc = clock();
    
    
    /* Record Time to Solve the Master Problem */
    soln->runTime->masterTime[soln->runTime->masterCnt] = ((double) (toc - tic)) / CLOCKS_PER_SEC;
    soln->runTime->masterCnt++;

	/* increment the number of problems solved during algorithm */
	cell->LPcnt++;

	/* record the objective function value */
	soln->optValM = getObjective(cell->master->lp, cell->master->type);

	/* Get the most recent optimal solution to master program */
	status = getPrimal(cell->master->lp, soln->candidX, num->cols);
	if ( status ) {
		errMsg("algorithm", "solveMaster", "failed to obtain the primal solution for master", 0);
		return 1;
	}

	/* add the incumbent back to change from \Delta X to X */
	for (i = 1; i <= num->cols; i++)
		d2 += soln->candidX[i] * soln->candidX[i];
    addVectors(soln->candidX, soln->incumbX, NULL, num->cols);
    
	/* update d_norm_k in soln_type. */
	if (cell->k == 1)
		soln->normDk_1 = d2;
	soln->normDk = d2;

	/* Get the dual solution too */
	status = getDual(cell->master->lp, soln->piM, num->rows+cell->cuts->cnt);
	if ( status ) {
		errMsg("solver", "solveQPMaster", "failed to obtain dual solutions to master", 0);
		return 1;
	}
	status = getDualSlacks(cell->master->lp, soln->djM, num->cols);
	if ( status ) {
		errMsg("solver", "solveQPMaster", "failed to obtain dual slacks for master", 0);
		return 1;
	}

	/* Find the highest cut at the candidate solution. where cut_height = alpha - beta(xbar + \Delta X) */
	if (cell->cuts->cnt > 0) {
		Sm = cutHeight(cell->lbType, cell->lb, cell->cuts->val[0], cell->k, soln->candidX, num->cols);
		for (i = 1; i < cell->cuts->cnt; i++) {
			ht = cutHeight(cell->lbType, cell->lb, cell->cuts->val[i], cell->k, soln->candidX, num->cols);
			if (Sm < ht)
				Sm = ht;
		}
	}else {
		if (cell->lbType == TRIVIAL)
			Sm += 0.0;
		else
			Sm += cell->lb;
	}
    
#ifdef TRACE
    printf("\t :: Highest cut at andidate solution with height = %f\n", Sm);
#endif
    
	soln->candidEst = Sm + vXvSparse(soln->candidX, dBar);

#ifdef TRACE
    printf("\t :: Updating candidate estimate by adding first stage cost... soln->candidEst = %f\n", soln->candidEst);
#endif

	/* Calculate gamma for next improvement check on incumbent x. If it is not solved in opt mode, gamma will not change */
	if (soln->optMode == TRUE)
		soln->gamma = soln->candidEst - soln->incumbEst;

	return 0;
}//END solveQPMaster()

/* Construct the Q diagonal matrix and copy it for quadratic problem. */
int constructQP(LPptr lp, int numCols, double sigma) {
	int status = 0, idx;
	double *qsepvec;

	if (!(qsepvec = arr_alloc(numCols+1, double)))
		errMsg("Allocation", "constructQP", "qsepvec",0);

	/* Construct Q matrix, which is simply a diagonal matrix. */
	for (idx = 0; idx < numCols; idx++)
		qsepvec[idx] = 0.5 * sigma;

	/* This is for eta column */
	qsepvec[numCols] = 0.0;

	/* Now copy the Q matrix for QP problem. */
	status = copyQPseparable(lp, qsepvec);
	if (status) {
		fprintf(stderr, "Failed to copy Q matrix.\n");
		return 1;
	}

	mem_free(qsepvec);
	return 0;
}//END constructQP

/* In the regularized QP method, we need to change the rhs of x to d. The
 * 		 A * x 			= b
 * 		 eta + beta * x >= alpha
 * Since x = xbar + d, the corresponding changes will therefore be:
 * 		 A * d = b - A * xbar
 * 		 eta + beta * d >= alpha - beta * xbar
 * But as long as the incumbent sulotion does not change, b - A * xbar and alpha - beta * xbar (for the existing cuts) won't change. So we only need
 * to change it when the incumbent changes.
 *
 * On the other hand, in each iteration, a new cut will be added (and/or some cuts may be dropped) and therefore we need to shift the rhs of the
 * added cut from _alpha_ to _alpha - beta * xbar_, which has taken care of in the routine addCut() in cuts.c. We do not need to worry about the shift
 * of rhs for the dropped cuts.
 * This function performs the change of rhs when the incumbent changes, as described above. */
int changeQPrhs(probType *prob, cellType *cell, vector xk) {
	int 	status = 0, cnt;
	vector 	rhs = NULL;
	intvec 	indices = NULL;

	if (!(rhs =(vector) arr_alloc(prob->num->rows+cell->cuts->cnt+1, double)))
		errMsg("Allocation", "changeRhs", "rhs",0);
	if (!(indices =(intvec) arr_alloc(prob->num->rows+cell->cuts->cnt, int)))
		errMsg("Allocation", "changeRhs", "indices",0);

	/* Be careful with the one_norm!! In the CxX() routine, it assumes the 0th element is reserved for the 1_norm, in the returned vector, the T sparse
	 vector, and the x vector. */
	for (cnt = 0; cnt < prob->num->rows; cnt++) {
		rhs[cnt + 1] = prob->sp->rhsx[cnt];
		indices[cnt] = cnt;
	}

	/* b - A * xbar */
	rhs = MSparsexvSub(prob->Dbar, xk, rhs);

	/*** new rhs = alpha - beta * xbar ***/
	for (cnt = 0; cnt < cell->cuts->cnt; cnt++) {
		rhs[prob->num->rows+cnt+1] = cell->cuts->val[cnt]->alpha - vXv(cell->cuts->val[cnt]->beta, xk, NULL, prob->sp->mac);
		indices[prob->num->rows+cnt] = cell->cuts->val[cnt]->rowNum;

		cell->cuts->val[cnt]->alphaIncumb = rhs[prob->num->rows+cnt+1];
	}

	/* Do the same thing for feasibility cut: new rhs = alpha - beta * xbar */
	for (cnt = 0; cnt < cell->fCutAdded->cnt; cnt++) {
		rhs[prob->num->rows + cell->cuts->cnt + cnt] = cell->fCutAdded->val[cnt]->alpha - vXv(cell->fCutAdded->val[cnt]->beta, xk, NULL, prob->sp->mac );
		indices[prob->num->rows + cell->cuts->cnt + cnt] = cell->fCutAdded->val[cnt]->rowNum;
	}

	/* Now we change the right-hand of the master problem. */
	status = changeRHS(cell->master->lp, prob->num->rows + cell->cuts->cnt, indices, rhs+1);
	if (status)	{
		errMsg("solver", "changeQPrhs", "failed to change the right-hand side in the solver", 0);
		return 1;
	}

	mem_free(rhs);
	mem_free(indices);
	return 0;
}//END changeQPrh

/* This function changes the (lower) bounds of the variables, while changing from x to d. The lower bounds of d varibles are -xbar
 * (incumbent solution). */
int changeQPbds(LPptr lp, int numCols, vector bdl, vector bdu, vector xk) {
	int 	status = 0, cnt;
	vector	lbounds, ubounds;
	intvec	lindices, uindices;
	char 	*llu, *ulu;

	if (!(lbounds = arr_alloc(numCols, double)))
		errMsg("Allocation", "changeBounds", "lbounds",0);
	if (!(lindices = arr_alloc(numCols, int)))
		errMsg("Allocation", "change_bounds", "lindices",0);
	if (!(llu = arr_alloc(numCols, char)))
		errMsg("Allocation", "changeBounds", "llu",0);

	if (!(ubounds = arr_alloc(numCols, double)))
		errMsg("Allocation", "change_bounds", "ubounds",0);
	if (!(uindices = arr_alloc(numCols, int)))
		errMsg("Allocation", "changeBounds", "uindices",0);
	if (!(ulu = arr_alloc(numCols, char)))
		errMsg("Allocation", "changeBounds", "ulu",0);

	/* Change the Upper Bound */
	for (cnt = 0; cnt < numCols; cnt++) {
		ubounds[cnt] = bdu[cnt] - xk[cnt + 1];
		uindices[cnt] = cnt;
		ulu[cnt] = 'U';
	}

	status = changeBDS(lp, numCols, uindices, ulu, ubounds);
	if (status) {
		errMsg("algorithm", "changeQP", "failed to change the upper bound in the solver", 0);
		return 1;
	}

	/* Change the Lower Bound */
	for (cnt = 0; cnt < numCols; cnt++) {
		lbounds[cnt] = bdl[cnt] - xk[cnt + 1];
		lindices[cnt] = cnt;
		llu[cnt] = 'L';
	}

	status = changeBDS(lp, numCols, lindices, llu, lbounds);
	if (status) {
		errMsg("algorithm", "changeQP", "failed to change the lower bound in the solver", 0);
		return 1;
	}

	mem_free(lbounds); mem_free(lindices); mem_free(llu);
	mem_free(ubounds); mem_free(uindices); mem_free(ulu);

	return 0;
}//END changeQPbds()

/* This function performs the updates on all the coefficients of eta in the master problem constraint matrix.  During every iteration,
 * each of the coefficients on eta are increased, so that the effect of the cut on the objective function is decreased. */
int changeEtaCol(LPptr lp, int numRows, int numCols, int k, cutType *cuts, double lb) {

#ifdef TRACE
    printf("\t\t~changeEtaCol();\n");
#endif
    
	vector	coef = NULL;
	double	etaCoef[1], etaBds[1];
	int 	status, c, etaCol[1];
	char	bdsType[1];

	etaCol[0] = numCols;
	bdsType[0] = 'L';

	/* array of coefficients for the \eta column. */
	if (!(coef = (vector) arr_alloc(cuts->cnt, double)))
		errMsg("allocation", "chgEtaCol", "coef", 0);

	for (c = 0; c < cuts->cnt; c++)
		/* Currently both incumbent and candidate cuts are treated similarly, and sunk as iterations proceed */
		coef[cuts->val[c]->rowNum - numRows] = (double) (k) / (double) cuts->val[c]->cutObs;

	/* change the eta column in the stage problem, which corresponds to the eta variable, starting with the row after the D matrix
	 * and ending with the row of the last cut. */
	status = changeCol(lp, numCols, coef, numRows, numRows+cuts->cnt);
	if ( status ) {
		errMsg("solver", "chgEtaCol", "failed to change eta column in the stage problem", 0);
		return 1;
	}

	/* if feasibility cut is added without any general cut, eta's lower bound should be changed to lower bound computed from mean value solution
	 * and its objective function coefficient should be zero. Once a general cut is encountered, revert to 1.0 for objective coefficient and -\infty
	 * for lower bound on eta */
	if ( cuts->cnt <= 1) {
		if ( cuts->cnt > 0 ) {
			etaCoef[0] = 1.0;
			etaBds[0]  = -INFBOUND;
		}
		else {
			etaCoef[0] = 0.0;
			etaBds[0] = lb;
		}
		status = changeObjx(lp, 1, etaCol, etaCoef);
		if ( status ) {
			errMsg("solver", "changeEtaCol", "failed to change the objective coefficient of eta column in objective function value", 0);
			return 1;
		}
		status = changeBDS(lp, 1, etaCol, bdsType, etaBds);
		if ( status ) {
			errMsg("solver", "changeEtaCol", "failed to change the bound for eta column", 0);
			return 1;
		}
	}

	mem_free(coef);
	return 0;
}//END chgEtaCol()

int updateRHS(cellType *cell) {

#ifdef TRACE
    printf("\t\t~updateRHS();\n");
#endif
	int 	cnt;
	vector	rhs = NULL;
	intvec	indices = NULL;

	if (!(rhs = arr_alloc(cell->cuts->cnt, double)))
		errMsg("allocation", "updateRHS", "rhs", 0);
	if (!(indices = arr_alloc(cell->cuts->cnt, int)))
		errMsg("allocation", "updateRHS", "indices", 0);

	for (cnt = 0; cnt < cell->cuts->cnt; cnt++) {
		rhs[cnt] = cell->cuts->val[cnt]->alphaIncumb;
		rhs[cnt] += ((double) cell->k / (double) cell->cuts->val[cnt]->cutObs - 1) * cell->lb;
		indices[cnt] = cell->cuts->val[cnt]->rowNum;
	}

	/* Now we change the right-hand of the master problem. */
	cnt = changeRHS(cell->master->lp, cell->cuts->cnt, indices, rhs);
	if (cnt)	{
		errMsg("solver", "changeQPrhs", "failed to change the right-hand side in the solver", 0);
		return 1;
	}

	mem_free(rhs);
	mem_free(indices);

	return 0;
}//END updateRHS


/* This function will allocate and initialize memory for a copy of the master problem passed to it as _master_. In addition to the original
 * master, it will add a new column, called "eta", whose coefficients are all zero and whose cost is one. It will also include all cuts which exist
 * in the _cuts_ structure as additional rows in the problem. Finally, it will allocate enough room in all relevant problem arrays to hold up to
 * _extra_cuts_ more constraints, assuming these constraints to be non-sparse (non-zero coefficients in every column).  It returns a pointer to
 * the newly created master problem.
 *
 * Note how carefully the coefficients are placed in the matval array. CPLEX requires them to be ordered by columns, but we know we will be adding
 * rows. So, each segment of the array corresponding to a given column contains enough empty locations at its end to hold all present and future
 * cuts (rows).  This way, when we add a row, CPLEX does not have to reorder and shift all the coefficients down. */
oneProblem *newMaster(oneProblem *master, vector initSol, sparseMatrix *initCbar, cutType *cuts, int extra_cuts) {
	oneProblem *copy;
	int 	r, i, j, idx, cnt, len;
	long 	colOffset, rowOffset;
	char 	cutName[NAMESIZE] = { "Old    " }; 	/* diff. from add_cut */
	char 	*q;
    vector  tempRHS = NULL;

	if (!(copy = (oneProblem *) mem_malloc (sizeof(oneProblem))))
		errMsg("allocation", "new_master", "copy", 0);

	/* Initialize dimensions of copy based on master and new cuts. */
    copy->name = NULL;
    copy->objname = NULL;
    copy->objx = NULL;
    copy->bdl = NULL;
    copy->bdu = NULL;
    copy->rhsx = NULL;
    copy->senx = NULL;
    copy->matbeg = NULL;
    copy->matcnt = NULL;
    copy->cname = NULL;
    copy->cstore = NULL;
    copy->rname = NULL;
    copy->rstore = NULL;
    copy->matval = NULL;
    copy->matind = NULL;
    copy->ctype = NULL;
    copy->lp = NULL;
	copy->matsz = master->matsz + (cuts->cnt + extra_cuts) * (master->mac + 1);
	copy->marsz = master->mar + cuts->cnt + extra_cuts;
	copy->mar = master->mar + cuts->cnt;
	copy->macsz = master->mac + 1;
	copy->mac = master->mac + 1;
	copy->cstorsz = master->cstorsz + NAMESIZE;
	copy->rstorsz = master->rstorsz + NAMESIZE * (cuts->cnt + extra_cuts);
	copy->objsen = master->objsen;
	copy->numInt = master->numInt;
	if ( run.MASTER_TYPE == PROB_QP && master->type == PROB_MILP )
		printf("Warning :: solving the relaxed MIP using 2-SD");
	else if ( run.MASTER_TYPE == PROB_LP ) {
		errMsg("setup", "newMaster", "requested to solve master as an LP, change the configuration file", 0);
		return NULL;
	}
	copy->type = run.MASTER_TYPE;

	/* Make all allocations of known sizes, as calculated above */
	if (!(copy->name = arr_alloc(NAMESIZE, char)))
		errMsg("allocation", "new_master", "copy->name",0);
	if (!(copy->objname = arr_alloc(NAMESIZE, char)))
		errMsg("allocation", "new_master", "copy->objname",0);
	if (!(copy->objx = arr_alloc(copy->macsz, double)))
		errMsg("allocation", "new_master", "copy->objx",0);
	if (!(copy->bdl = arr_alloc(copy->macsz, double)))
		errMsg("allocation", "new_master", "copy->bdl",0);
	if (!(copy->bdu = arr_alloc(copy->macsz, double)))
		errMsg("allocation", "new_master", "copy->bdu",0);
	if (!(copy->rhsx = arr_alloc(copy->marsz, double)))
		errMsg("allocation", "new_master", "copy->rhsx",0);
	if (!(copy->senx = arr_alloc(copy->marsz, char)))
		errMsg("allocation", "new_master", "copy->senx",0);
	if (!(copy->matbeg = arr_alloc(copy->macsz, int)))
		errMsg("allocation", "new_master", "copy->matbeg",0);
	if (!(copy->matcnt = arr_alloc(copy->macsz, int)))
		errMsg("allocation", "new_master", "copy->matcnt",0);
	if (!(copy->cname = arr_alloc(copy->macsz, string)))
		errMsg("allocation", "new_master", "copy->cname",0);
	if (!(copy->cstore = arr_alloc(copy->cstorsz, char)))
		errMsg("allocation", "new_master", "copy->cstore",0);
	if (!(copy->rname = arr_alloc(copy->marsz, string)))
		errMsg("allocation", "new_master", "copy->rname",0);
	if (!(copy->rstore = arr_alloc(copy->rstorsz, char)))
		errMsg("allocation", "new_master", "copy->rstore",0);
	if (!(copy->matval = arr_alloc(copy->matsz, double)))
		errMsg("allocation", "new_master", "copy->matval",0);
	if (!(copy->matind = arr_alloc(copy->matsz, int)))
		errMsg("allocation", "new_master", "copy->matind",0);
	if (!(copy->ctype = arr_alloc(copy->macsz, char)))			/*Added by zahra*/
		errMsg("allocation", "new_master", "copy->ctype",0);

	/* First copy information directly from the original master problem. */
	/* Copy the master problem's column and row names */
	/* Assume uninitialized elements are zero, or '\0', from calloc */
	i = 0;
	for (q = master->cname[0]; q < master->cname[0] + master->cstorsz; q++)
		copy->cstore[i++] = *q;

	i = 0;
	for (q = master->rname[0]; q < master->rname[0] + master->rstorsz; q++)
		copy->rstore[i++] = *q;

	strcpy(copy->name, master->name);
	strcpy(copy->objname, master->objname);

	/* Calculate difference in pointers for master/copy row and column names */
	colOffset = copy->cstore - master->cname[0];
	rowOffset = copy->rstore - master->rname[0];

	/* Copy the all column information from the original master problem */
	cnt = 0;
	for (j = 0; j < master->mac; j++) {
		copy->objx[j] = master->objx[j];
		copy->ctype[j] = master->ctype[j];
		copy->bdu[j] = master->bdu[j];
		copy->bdl[j] = master->bdl[j];
		copy->cname[j] = master->cname[j] + colOffset;
		copy->matbeg[j] = cnt;
		copy->matcnt[j] = master->matcnt[j];
		for (idx = master->matbeg[j]; idx < master->matbeg[j] + master->matcnt[j]; idx++) {
			copy->matval[cnt] = master->matval[idx];
			copy->matind[cnt] = master->matind[idx];
			cnt++;
		}
		cnt += cuts->cnt + extra_cuts;
	}

	/* Copy all information concerning rows of master */
	for (r = 0; r < master->mar; r++) {
		copy->rhsx[r] = master->rhsx[r];
		copy->senx[r] = master->senx[r];
		copy->rname[r] = master->rname[r] + rowOffset;
	}

	/* Initialize information for the extra column in the new master. */
	strcpy(copy->cstore + master->cstorsz, "eta");
	copy->cname[master->mac] = copy->cstore + master->cstorsz;
	copy->objx[master->mac] = 1.0;
	copy->ctype[master->mac] = 'C';
	copy->bdu[master->mac] = INFBOUND;
	copy->bdl[master->mac] = -INFBOUND;
	copy->matbeg[master->mac] = cnt;
	copy->matcnt[master->mac] = cuts->cnt;

	for (idx = 0; idx < cuts->cnt; idx++) {
		copy->matval[idx + cnt] = 1.0;
		copy->matind[idx + cnt] = master->mar + idx;
	}

	/* Now copy information from the cuts into the new master problem. */

	if (cuts->cnt) {
		/* Copy the constraint coefficients */
		for (j = 0; j < master->mac; j++)
		{
			cnt = copy->matbeg[j] + copy->matcnt[j];
			copy->matcnt[j] += cuts->cnt;

			for (i = 0; i < cuts->cnt; i++)
			{
				copy->matval[cnt + i] = cuts->val[i]->beta[j + 1];
				copy->matind[cnt + i] = master->mar + i;
			}
		}

		/* Give names to the cut constraints, and add rhs values */
		len = master->rstorsz;
		for (cnt = 0; cnt < cuts->cnt; cnt++)
		{
			cutName[3] = '0' + cnt / 1000 % 10;
			cutName[4] = '0' + cnt / 100 % 10;
			cutName[5] = '0' + cnt / 10 % 10;
			cutName[6] = '0' + cnt / 1 % 10;
			copy->rname[master->mar + cnt] = copy->rstore + len;
			strcpy(copy->rname[master->mar + cnt], cutName);
			len += strlen(cutName) + 1;
			copy->rhsx[master->mar + cnt] = cuts->val[cnt]->alpha;
			copy->senx[master->mar + cnt] = 'G';
		}
	}
    
    /* Updates right hand side only when encountering rolling horizon problems */
    if (run.HORIZON > 1 && initSol) {
        if (!(tempRHS = arr_alloc(copy->marsz+1, double)))
            errMsg("allocation", "new_master", "tempRHS",0);
        for (i=0; i<copy->marsz; i++) tempRHS[i+1] = copy->rhsx[i];
         /* modify the right-hand side with initial solution */
        if (initSol != NULL)
            tempRHS = MSparsexvSub(initCbar, initSol, tempRHS);
        for (i=0; i<copy->marsz; i++) copy->rhsx[i] = tempRHS[i+1];
    }

	/* Load the copy into CPLEX now */
	copy->lp = setupProblem(copy->name, copy->type, copy->mac, copy->mar, copy->objsen, copy->objx, copy->rhsx, copy->senx, copy->matbeg, copy->matcnt,
			copy->matind, copy->matval, copy->bdl, copy->bdu, NULL, copy->cname, copy->rname, copy->ctype);
	if ( copy->lp == NULL ) {
		errMsg("Problem Setup", "newMaster", "failed to setup master problem in the solver",0);
		return NULL;
	}
    
    /* Initial Master Problem */
#if defined(MODIFY)
    writeProblem(copy->lp, "initMaster.lp");
#endif
    
    if (tempRHS) mem_free(tempRHS);
	return copy;
}//END newMaster

oneProblem *newStoMaster(oneProblem *master, sparseMatrix *initCbar, cutType *cuts, int extra_cuts) {
    oneProblem *copy;
    int 	r, i, j, idx, cnt, len;
    long 	colOffset, rowOffset;
    char 	cutName[NAMESIZE] = { "Old    " }; 	/* diff. from add_cut */
    char 	*q;
    vector  tempRHS = NULL;
    
    if (!(copy = (oneProblem *) mem_malloc (sizeof(oneProblem))))
        errMsg("allocation", "new_master", "copy", 0);
    
    /* Initialize dimensions of copy based on master and new cuts. */
    copy->name = NULL;
    copy->objname = NULL;
    copy->objx = NULL;
    copy->bdl = NULL;
    copy->bdu = NULL;
    copy->rhsx = NULL;
    copy->senx = NULL;
    copy->matbeg = NULL;
    copy->matcnt = NULL;
    copy->cname = NULL;
    copy->cstore = NULL;
    copy->rname = NULL;
    copy->rstore = NULL;
    copy->matval = NULL;
    copy->matind = NULL;
    copy->ctype = NULL;
    copy->lp = NULL;
    copy->matsz = master->matsz + (cuts->cnt + extra_cuts) * (master->mac + 1);
    copy->marsz = master->mar + cuts->cnt + extra_cuts;
    copy->mar = master->mar + cuts->cnt;
    copy->macsz = master->mac + 1;
    copy->mac = master->mac + 1;
    copy->cstorsz = master->cstorsz + NAMESIZE;
    copy->rstorsz = master->rstorsz + NAMESIZE * (cuts->cnt + extra_cuts);
    copy->objsen = master->objsen;
    copy->numInt = master->numInt;
    if ( run.MASTER_TYPE == PROB_QP && master->type == PROB_MILP )
        printf("Warning :: solving the relaxed MIP using 2-SD");
    else if ( run.MASTER_TYPE == PROB_LP ) {
        errMsg("setup", "newMaster", "requested to solve master as an LP, change the configuration file", 0);
        return NULL;
    }
    copy->type = run.MASTER_TYPE;
    
    /* Make all allocations of known sizes, as calculated above */
    if (!(copy->name = arr_alloc(NAMESIZE, char)))
        errMsg("allocation", "new_master", "copy->name",0);
    if (!(copy->objname = arr_alloc(NAMESIZE, char)))
        errMsg("allocation", "new_master", "copy->objname",0);
    if (!(copy->objx = arr_alloc(copy->macsz, double)))
        errMsg("allocation", "new_master", "copy->objx",0);
    if (!(copy->bdl = arr_alloc(copy->macsz, double)))
        errMsg("allocation", "new_master", "copy->bdl",0);
    if (!(copy->bdu = arr_alloc(copy->macsz, double)))
        errMsg("allocation", "new_master", "copy->bdu",0);
    if (!(copy->rhsx = arr_alloc(copy->marsz, double)))
        errMsg("allocation", "new_master", "copy->rhsx",0);
    if (!(copy->senx = arr_alloc(copy->marsz, char)))
        errMsg("allocation", "new_master", "copy->senx",0);
    if (!(copy->matbeg = arr_alloc(copy->macsz, int)))
        errMsg("allocation", "new_master", "copy->matbeg",0);
    if (!(copy->matcnt = arr_alloc(copy->macsz, int)))
        errMsg("allocation", "new_master", "copy->matcnt",0);
    if (!(copy->cname = arr_alloc(copy->macsz, string)))
        errMsg("allocation", "new_master", "copy->cname",0);
    if (!(copy->cstore = arr_alloc(copy->cstorsz, char)))
        errMsg("allocation", "new_master", "copy->cstore",0);
    if (!(copy->rname = arr_alloc(copy->marsz, string)))
        errMsg("allocation", "new_master", "copy->rname",0);
    if (!(copy->rstore = arr_alloc(copy->rstorsz, char)))
        errMsg("allocation", "new_master", "copy->rstore",0);
    if (!(copy->matval = arr_alloc(copy->matsz, double)))
        errMsg("allocation", "new_master", "copy->matval",0);
    if (!(copy->matind = arr_alloc(copy->matsz, int)))
        errMsg("allocation", "new_master", "copy->matind",0);
    if (!(copy->ctype = arr_alloc(copy->macsz, char)))			/*Added by zahra*/
        errMsg("allocation", "new_master", "copy->ctype",0);
    
    
    
    /* First copy information directly from the original master problem. */
    /* Copy the master problem's column and row names */
    /* Assume uninitialized elements are zero, or '\0', from calloc */
    i = 0;
    for (q = master->cname[0]; q < master->cname[0] + master->cstorsz; q++) {
        copy->cstore[i++] = *q;
    }
    
    i = 0;
    for (q = master->rname[0]; q < master->rname[0] + master->rstorsz; q++)
        copy->rstore[i++] = *q;
    
    strcpy(copy->name, master->name);
    strcpy(copy->objname, master->objname);
    
    
    /* Calculate difference in pointers for master/copy row and column names */
    colOffset = copy->cstore - master->cname[0];
    rowOffset = copy->rstore - master->rname[0];
    
    /* Copy the all column information from the original master problem */
    cnt = 0;
    for (j = 0; j < master->mac; j++) {
        copy->objx[j] = master->objx[j];
        copy->ctype[j] = master->ctype[j];
        copy->bdu[j] = master->bdu[j];
        copy->bdl[j] = master->bdl[j];
        copy->cname[j] = master->cname[j] + colOffset;
        copy->matbeg[j] = cnt;
        copy->matcnt[j] = master->matcnt[j];
        for (idx = master->matbeg[j]; idx < master->matbeg[j] + master->matcnt[j]; idx++) {
            copy->matval[cnt] = master->matval[idx];
            copy->matind[cnt] = master->matind[idx];
            cnt++;
        }
        cnt += cuts->cnt + extra_cuts;
    }
    
    /* Copy all information concerning rows of master */
    for (r = 0; r < master->mar; r++) {
        copy->rhsx[r] = master->rhsx[r];
        copy->senx[r] = master->senx[r];
        copy->rname[r] = master->rname[r] + rowOffset;
    }
    
    /* Initialize information for the extra column in the new master. */
    strcpy(copy->cstore + master->cstorsz, "eta");
    copy->cname[master->mac] = copy->cstore + master->cstorsz;
    copy->objx[master->mac] = 1.0;
    copy->ctype[master->mac] = 'C';
    copy->bdu[master->mac] = INFBOUND;
    copy->bdl[master->mac] = -INFBOUND;
    copy->matbeg[master->mac] = cnt;
    copy->matcnt[master->mac] = cuts->cnt;
    
    for (idx = 0; idx < cuts->cnt; idx++) {
        copy->matval[idx + cnt] = 1.0;
        copy->matind[idx + cnt] = master->mar + idx;
    }
    
    /* Now copy information from the cuts into the new master problem. */
    
    if (cuts->cnt) {
        /* Copy the constraint coefficients */
        for (j = 0; j < master->mac; j++)
        {
            cnt = copy->matbeg[j] + copy->matcnt[j];
            copy->matcnt[j] += cuts->cnt;
            
            for (i = 0; i < cuts->cnt; i++)
            {
                copy->matval[cnt + i] = cuts->val[i]->beta[j + 1];
                copy->matind[cnt + i] = master->mar + i;
            }
        }
        
        /* Give names to the cut constraints, and add rhs values */
        len = master->rstorsz;
        for (cnt = 0; cnt < cuts->cnt; cnt++)
        {
            cutName[3] = '0' + cnt / 1000 % 10;
            cutName[4] = '0' + cnt / 100 % 10;
            cutName[5] = '0' + cnt / 10 % 10;
            cutName[6] = '0' + cnt / 1 % 10;
            copy->rname[master->mar + cnt] = copy->rstore + len;
            strcpy(copy->rname[master->mar + cnt], cutName);
            len += strlen(cutName) + 1;
            copy->rhsx[master->mar + cnt] = cuts->val[cnt]->alpha;
            copy->senx[master->mar + cnt] = 'G';
        }
    }
    
    /* Load the copy into CPLEX now */
    copy->lp = setupProblem(copy->name, copy->type, copy->mac, copy->mar, copy->objsen, copy->objx, copy->rhsx, copy->senx, copy->matbeg, copy->matcnt,
                            copy->matind, copy->matval, copy->bdl, copy->bdu, NULL, copy->cname, copy->rname, copy->ctype);
    if ( copy->lp == NULL ) {
        errMsg("Problem Setup", "newMaster", "failed to setup master problem in the solver",0);
        return NULL;
    }
    
    /* Initial Master Problem */
#if defined(MODIFY)
    writeProblem(copy->lp, "initMaster.lp");
#endif
    
    if (tempRHS) mem_free(tempRHS);
    return copy;
}//END newMaster

oneProblem *newDETMaster(oneProblem *detProb) {
    
    oneProblem *copy;
    int r;
    
    if (!(copy = mem_malloc(sizeof(oneProblem)))) {
        errMsg("allocation","newCell","failed to allocate memory to deterministic oneProblem",0);
        return NULL;
    }
    
    /* Initialization of every member in deterministic master problem */
    copy->name = NULL;
    copy->objname = NULL;
    copy->objx = NULL;
    copy->bdl = NULL;
    copy->bdu = NULL;
    copy->rhsx = NULL;
    copy->senx = NULL;
    copy->matbeg = NULL;
    copy->matcnt = NULL;
    copy->cname = NULL;
    copy->cstore = NULL;
    copy->rname = NULL;
    copy->rstore = NULL;
    copy->matval = NULL;
    copy->matind = NULL;
    copy->ctype = NULL;
    copy->matsz = 0;
    copy->marsz = 0;
    copy->mar = detProb->mar;
    copy->macsz = 0;
    copy->mac = detProb->mac;
    copy->cstorsz = 0;
    copy->rstorsz = 0;
    copy->objsen = 0;
    copy->numInt = 0;
    copy->type = PROB_LP;
    
    /* Copy all information concerning rows of master */
    if (!(copy->rhsx = arr_alloc(detProb->marsz, double)))
        errMsg("allocation", "new_master", "copy->rhsx",0);
    if (!(copy->senx = arr_alloc(detProb->marsz, char)))
        errMsg("allocation", "new_master", "copy->senx",0);
    for (r = 0; r < detProb->mar; r++) {
        copy->rhsx[r] = detProb->rhsx[r];
        copy->senx[r] = detProb->senx[r];
    }
    
    copy->lp = setupProblem(detProb->name, detProb->type, detProb->mac, detProb->mar, detProb->objsen,
                                 detProb->objx, detProb->rhsx, detProb->senx,detProb->matbeg, detProb->matcnt,
                                 detProb->matind, detProb->matval, detProb->bdl, detProb->bdu, NULL,
                                 detProb->cname, detProb->rname, detProb->ctype);
    
    if (copy->lp == NULL) {
        errMsg("setup problem", "newDETMaster", "copy->lp", 0);
        return NULL;
    }
    
    return copy;
}


/* The master being used by each cell is just a copy of the real master which is stored in prob.  Once a cell is done, its copy
 ** (which includes all the cell's cuts) may be freed.  This function unloads the problem from the CPLEX and frees the memory. */
void freeMaster(oneProblem *copy) {

	freeProblem(copy->lp);
	freeOneProblem(copy);

}//END freeMaster
