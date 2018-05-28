                   //
//  modify.c
//  rollingSP
//
//

#include "cell.h"
#include "master.h"
#include "output.h"
#include "optimal.h"

extern runType  run;
extern string   outputDir;
extern string   evalsolDir;

int solveSDCell(stocType *stoc, probType **prob, cellType *cell, solnType *soln) {

#ifdef TRACE
	printf("~solveCell()\n");
#endif

	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     There are 3 observation vectors.
     When NOT USING time series model, observ is the same as fullObs, auxObs will be all zero.
     Because prob[1]->bBar still maintaines mean values in such circumstance, we need to keep
     fullObs the same as observ because solveSubprob() restores the values back using bBar.

     When USING time series model, observ is a noise tracker and fullObs is for subproblem right
     hand side. This is because warmUp() will modified prob[1]->bBar's random rows to be 0.

     Index  Variable   What about it?

       1     observ  - Residual Observation -> deltaType  White noise

       2     fullObs - Full Observation 
                        -> Subproblem RHS -> Historical + Recursive + Auxilary + Residual

       3     auxObs  - Dependent Auxiliary Observation
                        -> Calculated throught (fullObs) - (Residual) - (Historical + Recursive)
                                                                        |    Stored in gamma   |
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

	vector      observ = NULL;
	vector      fullObs = NULL;
	vector      auxObs = NULL;
    vector      trimObs = NULL;
	int         status = 0, omegaIdx, obsLength;
	BOOL        newOmegaFlag = TRUE;
    
#if defined(TSVALID)
    int         i;
#endif

#ifdef OMEGA
	/* When OMEGA is on, system will record all simulated fullObs into one file */
	FILE *omegaF = NULL;
	char omegaFName[BLOCKSIZE * 2];
	status = sprintf(omegaFName, "omega_%d.obs", cell->t);
    omegaF = openFile(outputDir, omegaFName, "w");
#endif

#ifdef SUBOBJ
    FILE *subObjF = NULL;
    char subObjFName[BLOCKSIZE * 2];
    status = sprintf(subObjFName, "subObj.out");
    subObjF = openFile(outputDir, subObjFName, "w");
#endif

	printf("---------------------------------------\n");
	printf("~~~~~~~~~ Time Period :: %d ~~~~~~~~~~~\n", cell->t);
	printf("---------------------------------------\n");

	/* Allocate memory to observation vectors */
	if ( !strcmp(stoc->type, "ARIMA") || cell->arima != NULL ) {
		// Time Series Model :: Obsevation Space = (# of Variates) * (# of Sub-hourly Steps) + 1
		obsLength = (cell->arima->subh)*(cell->arima->var) + 1;
	} else {
		// STOC File :: Observation Space = (# of r.v. in subproblem)
		obsLength = prob[1]->num->numRV + 1;
	}
	if ( !(fullObs = (vector) arr_alloc(obsLength, double)) ) {
		errMsg("allocation", "solveCell", "fullObs",0);
		return 1;
	}
	if ( !(observ = (vector) arr_alloc(obsLength, double)) ) {
		errMsg("allocation", "solveCell", "observ",0);
		return 1;
	}
	if ( !(auxObs = (vector) arr_alloc(obsLength, double)) ) {
		errMsg("allocation", "solveCell", "auxObs",0);
		return 1;
	}
    if ( !(trimObs = (vector) arr_alloc(obsLength, double)) ) {
        errMsg("allocation", "solveCell", "trimObs", 0);
        return 1;
    }

	/******************** Algorithm Start From Here ********************/
	while (cell->k < run.MAX_ITER * cell->t || (run.WARM_UP == 0 && cell->k < run.MAX_ITER)) {

        /******* 0. termination when running without warm start *******/
        if (run.WARM_UP == 0 && cell->k >= run.MAX_ITER)
            break;
        
		/******* 1. Optimality tests *******/
		soln->optFlag = optimal(prob, cell, soln);
		if (soln->optFlag)
			break;

		/******* 0. Entering New Iteration *******/
		cell->k++;

		/******* ITERATIVE COUNTER ******/
		if ( cell->k % run.PRINT_CYCLE == 0 || cell->k == (cell->t-1)*run.MAX_ITER + 1 )
			printf("\nIteration-%d : \n", cell->k);

		/******* 2. Generate new observation *******/
		generateOmega(cell, stoc, observ, fullObs, auxObs, trimObs, NULL, &run.RUN_SEED); //Maybe a different name for fullObs
        
		omegaIdx = calcOmega(observ-1, auxObs-1, trimObs-1, 0, prob[1]->num->numRV, cell->omega, cell->auxOmega, cell->trimOmega, &newOmegaFlag);
		//NOTE : omega->vals[idx] starts from 1... Omega keeps residual/observ start in position 1
        
#if defined(TSVALID)
        for (i=1; i<=cell->arima->subh*cell->arima->var; i++)
            if (fabs((auxObs[i-1] + cell->gamma->endo->bARMA->val[i] + observ[i-1] + trimObs[i-1]) - (fullObs[i-1])) > run.EPSILON)
                errMsg("TSERROR", "solveSDCell", "TS_VALIDATION ERROR", 0);
#endif
        
#if defined(OMEGA)
        printVector(fullObs-1,stoc->numOmega, NULL);
        printVector(fullObs-1,stoc->numOmega, omegaF);
#endif

		/******* 3. setup and solve subproblem *******/
		status = solveSubprob(prob[1], cell->subprob, soln, soln->candidX, fullObs-1, &soln->mubBar, &soln->subFeasFlag);
		if ( status ) {
			errMsg("algorithm", "solveCell", "failed to solve subproblem.\n", 0);
			return 1;
		}

#if defined(SUBOBJ)
        fprintf(subObjF,"%.3f\n", getObjective(cell->subprob->lp, PROB_LP));
#endif

		if ( !(soln->subFeasFlag) ) {
			/******* 4. resolve infeasibility *******/
			/* candidate solution is infeasible at subproblem: entering feasibility mode */
			// Keep a break point here for debugging, we should never go here.
			status = resolveInfeasibility(prob, cell, soln, &newOmegaFlag, omegaIdx);
			if ( status ) {
				errMsg("algorithm", "solveCell", "failed to resolve infeasibility", 0);
				return 1;
			}
		}

		/******* 5. form and update optimality cut *******/
		status = formOptCut(prob, cell, soln, newOmegaFlag, omegaIdx);
		if ( status ) {
			errMsg("algorithm", "solveCell", "failed to create a new optimality cut", 0);
			return 1;
		}
		/* change the new omega flag to reflect that stochastic updates have been completed */
		newOmegaFlag = FALSE;

		/******* 6. form incumbent cut *******/
		if ( cell->k > 1 && cell->cuts->cnt > 0 ) {
			if ( cell->k - soln->iCutUpdt == run.TAU ) {
				status = rFormIncumbCut(prob, cell, soln, omegaIdx, fullObs-1, auxObs-1, newOmegaFlag);
				if ( status ) {
					errMsg("algorithm", "solveCell", "failed to add incumbent cut", 0);
					return 1;
				}
			}
#if defined(SUBOBJ)
            fprintf(subObjF,"%.3f\n", getObjective(cell->subprob->lp, PROB_LP));
#endif
		}
		else {
			status = constructQP(cell->master->lp, prob[0]->num->cols, cell->quadScalar);
			if ( status ) {
				errMsg("algorithm", "solveCell", "failed to change the proximal term", 0);
				return 1;
			}
			status = changeQPrhs(prob[0], cell, soln->incumbX);
			if ( status ) {
				errMsg("algorithm", "solveCell",
						"failed to change the right-hand side to convert the problem into QP", 0);
				return 1;
			}
			status = changeQPbds(cell->master->lp, prob[0]->num->cols, prob[0]->sp->bdl,
					prob[0]->sp->bdu, soln->incumbX);
			if ( status ) {
				errMsg("algorithm", "solveCell",
						"failed to change the bounds to convert the problem into QP", 0);
				return 1;
			}
#if defined(MODIFY)
			writeProblem(cell->master->lp, "initQPMaster.lp");
#endif
		}

		/******* 7. check improvement in predicted values at candidate solution *******/
		if ( !(soln->incumbChg) && cell->k > 1)
			/* If the incumbent has not changed in the current iteration */
			checkImprovement(prob, cell, soln);

        /******* 8. Solve master problem ******/
		status = solveQPMaster(prob[0]->num, prob[0]->dBar, cell, soln);
		if ( status ) {
			errMsg("algorithm", "solveCell", "failed to solve master problem", 0);
			return 1;
		}
        
#if defined(DEBUGGING)
        printf("%d, %.2f, %s, %s, %.0f, %d, %d, %d, %.3f, %d, %d, %.2f, %.2f, %.2f, %.2f, %.2f \n",
               cell->k, cell->quadScalar,
               soln->dualStableFlag ? "TRUE" : "FALSE",
               soln->preCheckEverFlag ? "TRUE" : "FALSE",
               soln->gamma,
               cell->lambda->cnt, cell->sigma->cnt, cell->sigma->aCnt,
               soln->normDk_1,
               soln->cCutIdx, soln->iCutIdx,
               soln->optValM, vXvSparse(soln->candidX, prob[0]->dBar), soln->incumbX[0],
               soln->candidEst, soln->incumbEst);
#endif
        
	}

#if defined(OMEGA)
	fclose(omegaF);
#endif
    
#if defined(SUBOBJ)
    fclose(subObjF);
#endif

	mem_free(observ);
	mem_free(fullObs);
	mem_free(auxObs);
	mem_free(trimObs);
	return 0;
}//END algoIntSD()

/* 
 * A deterministic problem is solved here. Solving a deterministic problem is eaiser comparing to the stochastic version.
 * From one time period to another, we need to update the following information in the system
 *  - 1 - Dynamic Information of both "first stage" primal solution and "second stage" right hand side
 *  - 2 - Stochastic Information for random rows in the "second stage"
 * Instead of true randomness, we update the "second stage" right hand side with historical data as the best forecast we
 * considered.
 * The exogenous dynamic is still read from the .dyn file.
 * The main difference here is the problem structure. Instead of having a two-stage structure, we only decompose the problem
 * in a dummy stage that is no much use for us and a main stage that conataines the whole problem. Given this information,
 * the way we update the problem through differen time period need to be carefully constructed.
 */
int solveDETCell(cellType *cell, probType *probM, solnType *soln) {

#ifdef TRACE
	printf("~solveDETCell()\n");
#endif

	vector  xt;
	int     i;
	int     stat1 = 0, status = 0;

	printf("---------------------------------------\n");
	printf("  ~~~~~~~~ Time Period :: %d ~~~~~~~~  \n", cell->t);
	printf("---------------------------------------\n");

	if ( !(xt = arr_alloc(probM->sp->mac+1, double)) )
		errMsg("allocation", "solveDETCell", "xt", 0);

	static int detCnt;
	char detProbfName[30] = "/detProblem/detProblem   .lp";
	detProbfName[22] = '0' + detCnt / 100 % 10;
	detProbfName[23] = '0' + detCnt / 10 % 10;
	detProbfName[24] = '0' + detCnt / 1 % 10;
	detCnt++;
	writeProblem(cell->master->lp, detProbfName);

	/* Step 1. Solve the deterministic problem */
	status = solveProblem(cell->master->lp, probM->sp->name, probM->sp->type, &stat1);
	if (status) {
		errMsg("algorithm", "solveDETCell", "Failed to solve the deterministic cell", 0);
		return 1;
	}

	/* Step 2. Get Objective Value */
	soln->optValM = getObjective(cell->master->lp, probM->sp->type);
	printf(" :: Problem Solved with objective function value = %f\n", soln->optValM);

	/* Step 3. Get the primal soltuion fo the deterministic problem */
	status = getPrimal(cell->master->lp, xt, probM->sp->mac);
	if (status) {
		errMsg("algorithm", "solveDETCell", "Failed to obtain the primal solution for deterministic cell", 0);
		return 1;
	}
	printf(" :: Solution oneNorm = %f\n", xt[0]);

	/* Step 4. Parse the solution to fit soln->incumbX */
	for (i=0; i<=probM->num->cols; i++) {
		soln->incumbX[i] = xt[i];
	}

	mem_free(xt);
	return 0;
}//END of solveDETCell()

int injectSolution(cellType *cell, probType *probM, solnType *soln){
    
#ifdef TRACE
    printf("~injectSolution()\n");
#endif
    
    char        extsolfile[BLOCKSIZE];
    FILE        *fptr = NULL;
    float       solVal;
    int         i = 1;
    
    /* Locate the solution file */
    sprintf(extsolfile, "%s", evalsolDir);
    printf("%s\n",extsolfile);
    
    fptr = fopen(extsolfile ,"r");
    if ( fptr == NULL ) {
        errMsg("read", "injectSolution", "failed to open external solution file", 0);
        return 1;
    }
    
    /* Step 1. Read eavluator solution file */
    while( fscanf(fptr, "%f\n", &solVal) > 0 ) {
        soln->incumbX[0] += solVal;
        soln->incumbX[i] = solVal;
        i++;
        if (i >= probM->num->cols) {
            errMsg("read", "injectSolution", "unexcepted number exceeding the problem column length", 0);
            return 1;
        }
    }
    
    fclose(fptr);
    
    return 0;
}//END of injectSolution

/********************************************************************************************
 This function creates a brand spankin new cell structure. By that I mean that the structures 
 within the cell, like lambda, sigma, and cuts, do not have any data already initialized. 
 For a single cell problem, or for first generation ancestors, this is standard; but for the
 parallel implementation, the new cell must be read in from a file, and its individual 
 structures will not be empty.
 *********************************************************************************************/
cellType *newCell(probType **prob, stocType *stoc) {
	cellType *c = NULL;
	int length = 0;

	if (!(c = (cellType *) mem_malloc (sizeof(cellType))))
		errMsg("allocation", "newCell", "cell",0);

    /* scalar elements in cellType */
	c->t = 0;											/* beginning of the horizon */
	c->k = 0; 											/* Cell starts with no cuts of its own */
	c->LPcnt = 0; 										/* Cell starts with no LP solved.*/
	c->feasCnt = 0;										/* number of infeasible subproblems */
	c->quadScalar = run.MIN_QUAD_SCALAR; 				/* The quadratic scalar, 'sigma'*/

	c->maxIter = run.MAX_ITER;
	c->minIter = run.MIN_ITER;
	c->tau = run.TAU;

	/* cell lower bounds */
	c->lb = prob[0]->lb;
	if ( c->lb  == 0 )
		c->lbType = TRIVIAL;
	else
		c->lbType = NONTRIVIAL;

	/* Setup all the stochastic elements used by the algorithm */
	length      = (c->maxIter + c->maxIter / c->tau + 1) * run.HORIZON;

	if (run.DETERMINISTIC == 1) {

		c->maxCuts = run.CUT_MULT * MAX_CUTS(prob[0]->num->cols);
		c->cuts = NULL;
		c->fCutAdded = NULL;
		c->fCutPool = NULL;
		c->fRow = 0;
		c->fCol = 0;

		/* Constructing a determinsitic solving cell */
		c->lambda   = NULL;
		c->sigma    = NULL;
		c->delta    = NULL;
		c->auxDelta = NULL;
        c->trimDelta= NULL;
		c->omega    = NULL;
        c->auxOmega = NULL;
        c->trimOmega= NULL;
		c->gamma    = newGamma(run.HORIZON+1, prob[0]->num->cols, prob[0]->num->cols, prob[0]->bBar->cnt, length);
		c->arima    = NULL;
		c->master   = newDETMaster(prob[0]->sp);
		if ( c->master == NULL ) {
			errMsg("Problem Setup", "newCell", "deterministic master",0);
			return NULL;
		}
        c->subprob  = NULL;
        c->initPrimal = NULL; //Delete this later
        

	} else {

		/* initialize optimality and feasibility cuts */
		// Be careful with this :: Rolling horizon problem doesn't change the maximum cuts in master problem
		c->maxCuts = run.CUT_MULT * MAX_CUTS(prob[0]->num->cols);
		c->cuts = newCuts(c->maxCuts);
		c->fCutAdded = newCuts(c->maxCuts);
		c->fCutPool = newCuts(c->maxCuts);
		c->fRow = 0;
		c->fCol = 0;

		/* Constructing a stochastic solving cell */
		// Certain parameters need to be adjusted to fit rolling horizon problem //
		c->lambda   = newLambda(length, 0, prob[1]->num->rvRowCnt, prob[1]->coord);
		c->sigma    = newSigma(length, prob[1]->num->rvColCnt, 0, prob[1]->coord);
		c->delta    = newDelta(length, prob[1]->coord);
		c->auxDelta = newDelta(length, prob[1]->coord);
        c->trimDelta= newDelta(length, prob[1]->coord);
		c->omega    = newOmega(c->maxIter * run.HORIZON + run.HORIZON);
        c->auxOmega = newOmega(c->maxIter * run.HORIZON + run.HORIZON);
        c->trimOmega= newOmega(c->maxIter * run.HORIZON + run.HORIZON);
		c->gamma    = newGamma(run.HORIZON+1, prob[0]->num->cols, prob[1]->num->cols, prob[1]->bBar->cnt, length);

		c->arima    = NULL;
        
		/* Setup master problem */
		c->master = newStoMaster(prob[0]->sp, prob[0]->Cbar, c->cuts, c->maxCuts);
		if ( c->master == NULL ) {
			errMsg("setup", "newCell", "failed to setup the master problem", 0);
			return NULL;
		}
        
		/* Setup subproblem */
		c->subprob = newSubprob(prob[1]->sp);
		if ( c->subprob == NULL ) {
			errMsg("setup", "newCell", "failed to setup the subproblem", 0);
			return NULL;
		}
        
        c->initPrimal = NULL;
	}

#if defined(MODIFY)
        if ( writeProblem(c->master->lp, "initMaster.lp") ) {
        	errMsg("solver", "newCell", "failed to write cell master problem", 0);
        	return NULL;
        }
        if ( c->subprob != NULL)
        	if ( writeProblem(c->subprob->lp, "initSubprob.lp") ) {
        		errMsg("solver", "newCell", "failed to write cell subproblem", 0);
        		return NULL;
        	}
#endif

	return c;
}//END newCell()

/****************************************************************************************************************************************
 This subroutine is designed to simulate ARMA observation. It takes three holders in, simulate, and carry out these holders with values.
 Note that the simulation use the deterministic information stored in tsType rather than gammaType. To simulate arma observations, one need
 to simulate noise first. With the simulated noise, it is not neceaary to simulate from all the begining given p and q value. Certain
 historical based knowledge can be used directly to speed up the simualtion. But notice, these historical information don't count for the
 dependent results that brought with newly simulated noise, which is the recursive parts. Hence, with the base historical infomration, we
 can start our simulation from there. Especially when in sub-hourly model(multiple time-scale model), an observation is decomposed into
 four parts.
 ARMA(p,q) = Histocial Deterministic Value + Recursive Deterministic Value + Dependent Auxilary Value + Residual Value
 We use the first part in this subroutine and simulate the rest, then calculate the dependent auxilary value and finally send out
 obs -> FULL OBSERVATION  ||    auxObs -> Dependent Auxilary Value    ||   noise -> Residual Value
 ****************************************************************************************************************************************/
int generateARIMA (int t, cellType *cell, vector obs, vector noise, vector auxObs, vector trimObs, vector inputNoise, long long *seed) {
    
#ifdef TRACE
    printf("\t\t~generateARIMA()\n");
#endif
    
    int i, j, k;
    int status;
    
    vector *tsObs;          //Stores time series simulation results
    vector *white;          //Stores generated noise
    
    /* Prepare the observation vector. Each variate in time series model have multiple sub-houly observation points */
    /* This raw observation vector is arranged for the ease of time series calculation. */
    if ( !(tsObs = (vector *) arr_alloc(cell->arima->subh+1, vector)))
        errMsg("allocation", "algoIntSD", "truObs", 1);
    
    white = (double **)malloc((cell->arima->subh+1)*sizeof(double*));
    
    for ( i=1; i<=cell->arima->subh; i++ ) {
        if (!(tsObs[i] = arr_alloc(cell->arima->var+1, double)))
            errMsg("allocation", "solveCell", "truObs",1);
        if (!(white[i] = arr_alloc(cell->arima->var+1, double)))
            errMsg("allocation", "solveCell", "truObs",1);
    }
    
    if (inputNoise == NULL) {
        /* Simulate Noise Vector */
        status = simulateNoise(cell->arima, white, seed);
        if (status) {
            errMsg("simulate","simulate_obs","Failed to simulate noise",0);
            return 1;
        }
    } else {
        /* If noise is mannully given, generate observation using that. */
        k=1;
        for (i=1; i<=cell->arima->subh; i++) {
            for (j=1; j<=cell->arima->var; j++) {
                white[i][j] = inputNoise[k];
                k++;
            }
        }
    }
    
    
    /* Generate new observations stored in obs matrix indexing using obs[n][var] */
    for (i=1; i<=cell->arima->subh; i++) {
        /* Note here, for every time period, we need to go through \Phi_1 -> \Phi_p.
         But some parts of this walk thourgh is already calulated.
         We take the calculate part(multivaritate->vector),
         calculate the new part (multivariate->vector), and add the noise to it to obtain the obs. */
        //----->AR model part = \sum_{i=1}^{p} Phi[i] * Omega[t-i]
        if (i == 1) {
            if (cell->gamma == NULL) {
                /* Calculate throug looping all p */
                for (j=1; j<cell->arima->p; j++)
                    tsObs[i]= MSparsexvAdd(cell->arima->Phi[j], cell->arima->history->obs[cell->arima->history->open-j+1], tsObs[i]);
            } else {
                /* The first case move values in detObs to obs directly */
                for (j=1; j<=cell->arima->var; j++) tsObs[i][j] = cell->arima->histARObs[i][j]; //Loop through all variates
            }
        } else {
            /* Calculate the part this is dependented on newly generated observations */
            for (j=1; j<=cell->arima->p; j++) {
                if ( j<i ) {
                    tsObs[i] = MSparsexvAdd(cell->arima->Phi[j], tsObs[i-j], tsObs[i]);
                } else {
                    if (cell->gamma == NULL) {
                        tsObs[i] = MSparsexvAdd(cell->arima->Phi[j], cell->arima->history->obs[cell->arima->history->open+i-j+1], tsObs[i]);
                    } else {
                        /* Once touching the deterministic part, use value from detObs in arima */
                        vPlusv(tsObs[i],cell->arima->histARObs[i],1,cell->arima->var);
                        break;
                    }
                }
            }
        }
        //----->MA model part
        if (i == 1) {
            if (cell->gamma == NULL) {
                for (j=1; j<cell->arima->q; j++)
                    tsObs[i] = MSparsexvAdd(cell->arima->Theta[j], cell->arima->noise->obs[cell->arima->history->open-j+1], tsObs[i]);
            } else {
                /* The first case move values in detObs to obs directly */
                for (j=1; j<=cell->arima->var; j++) tsObs[i][j] += cell->arima->histMAObs[i][j]; //Loop through all variates
            }
        } else {
            /* Calculate the part this is dependented on newly generated observations */
            for (j=1; j<=cell->arima->q; j++) {
                if ( j<i ) {
                    tsObs[i] = MSparsexvAdd(cell->arima->Theta[j], white[i-j], tsObs[i]); // MA part, this need modification
                } else {
                    if (cell->gamma == NULL) {
                        tsObs[i] = MSparsexvAdd(cell->arima->Theta[j], cell->arima->noise->obs[cell->arima->history->open+i-j+1], tsObs[i]);
                    } else {
                        /* Once touching the deterministic part, use value from detObs in arima */
                        vPlusv(tsObs[i],cell->arima->histMAObs[i],1,cell->arima->var);
                        break;
                    }
                }
            }
        }
        
        // Add Noise with consideration of preprocess  -> x_hat_t = ALL_OTHER + epsilon_t <-
        for (j=1; j<=cell->arima->var; j++) {
            tsObs[i][j] += white[i][j];
        }
    }
    
    
    k=0;  //Output position should start from 0
    /* Reverse truObs Results to real observation */
    for (i=1; i<=cell->arima->subh; i++) {
        for (j=1; j<=cell->arima->var; j++) {
            
            tsObs[i][j] = reverse(tsObs[i][j], j, i, t, cell->arima);
            
            /* Position simulated results in the way that optimization model can read */
            /* !!IMPORTANT!!: To do this, it is important to keep varieties .arma file in the right order. */
            // Output Full Observation : the right hand side value
            obs[k] = tsObs[i][j] * run.SCALAR; //Scale exit
            
            // Output Only Noise :: No reverse pre-process is conducted
            if (noise != NULL) {
                noise[k] = white[i][j] * run.SCALAR; //Scale exit
            }
            
            // Output intermedia contribution as an aux vector for later cut generation
            if (auxObs != NULL) {
                auxObs[k] = tsObs[i][j];
                auxObs[k] -= white[i][j];
                auxObs[k] -= (cell->gamma->endo->bARMA->val[k+1]);
            }
            
            // Record the trimmed value if pointer exist
            if (trimObs) {
                if (tsObs[i][j] < cell->arima->lb) {
                    trimObs[k] = cell->arima->lb - tsObs[i][j];
                }else {
                    trimObs[k] = 0.0;
                }
            }
            
            // Safety scheme
            if (obs[k] < cell->arima->lb)
                obs[k] = cell->arima->lb;
            
            k++;
        }
    }
    
    for (i=1; i<=cell->arima->subh; i++) {
        mem_free(tsObs[i]);
        mem_free(white[i]);
    }
    mem_free(tsObs);
    mem_free(white);
    
    
    return 0;
}


/* This function releases all the memory reserved for the cell data structure, including the cell itself.*/
void freeCellType(cellType *cell) {

    if ( cell ) {
		if (cell->delta) freeDeltaType(cell->delta, cell->lambda->cnt, cell->omega->cnt);
		if (cell->auxDelta) freeDeltaType(cell->auxDelta, cell->lambda->cnt, cell->omega->cnt);
        if (cell->trimDelta) freeDeltaType(cell->trimDelta, cell->lambda->cnt, cell->omega->cnt);
		if (cell->lambda ) freeLambdaType(cell->lambda);
		if (cell->sigma ) freeSigmaType(cell->sigma);
		if (cell->omega) freeOmegaType(cell->omega);
        if (cell->auxOmega) freeOmegaType(cell->auxOmega);
        if (cell->trimOmega) freeOmegaType(cell->trimOmega);
		if (cell->master) freeMaster(cell->master);
		if (cell->cuts)	freeCuts(cell->cuts);
		if (cell->fCutAdded) freeCuts(cell->fCutAdded);
		if (cell->fCutPool) freeCuts(cell->fCutPool);
		if (cell->gamma) freeGammaType(cell->gamma);
		if (cell->arima) freeTsType(cell->arima);
		if (cell->initPrimal) mem_free(cell->initPrimal);

		mem_free(cell);
	}

}//END freeCellType()

