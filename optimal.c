//
//  optimality.c
//  rollingSP
//
//
#include <time.h>
#include "optimal.h"
#include "output.h"
#include "arma.h"
#include "master.h"

extern runType run;
extern string  outputDir;

/***************************************************************************************
 This Subroutine is desgined to validate lower bound on warmed cuts between time periods.
 And how do we do it...
 ***************************************************************************************/
int lowerBoundCheck(probType **prob, cellType *cell, solnType *soln, BOOL *check) {
    
#ifdef OPTIMAL
    printf("\t~bootStrapCheck()\n");
#endif
    
    int     status, status2;
    int     i, obs;
    int     cutID = 0;
    int     masterRowCnt;
    double  Est=0.0, smEst=0.0, obj=0.0;
    vector  fullObs;
    vector  rhs = NULL;
    intvec  indices;

    
    masterRowCnt = getNumRows(cell->master->lp);

    printf("\t :: Current master has %d ROWs. \n", masterRowCnt);

    /* Memory allocation for calculation vectors */
    if ( !(fullObs = (vector) arr_alloc((cell->arima->subh)*(cell->arima->var) + 1, double)) )
        errMsg("validating", "bootStrapCheck", "fullObs", 0);
    if ( !(indices = (intvec) arr_alloc(prob[1]->num->rows, int)) )
        errMsg("validating", "solveSubporb", "ind", 0);
    for(i = 0; i < prob[1]->sp->mar; i++) indices[i] = i;
    
    printf("CutID cutObs OmegaCnt LBound UBound Gap \n");
    
    // Loop through all cuts
    for (i=0; i<cell->cuts->cnt; i++) {
        
        cutID = i;
        Est = 0.0;
        smEst = 0.0;
        
        // Compute the cutHeight at current incumbent solution
        Est = cutHeight(cell->lbType, cell->lb, cell->cuts->val[cutID], cell->k, soln->incumbX, prob[0]->num->cols);
        
        // Loop through all observations that is associated with this cut
        for (obs=0; obs<cell->cuts->val[cutID]->omegaCnt; obs++) {
            
            // Regenerate Observation (Hourly Model/Sub-Hourly Model)
            generateARIMA(cell->t, cell, fullObs, NULL, NULL, NULL, cell->omega->vals[obs], NULL);

            // Compute RHS Value of the right-hand side value
            rhs = computeRHS(prob[1]->num, prob[1]->coord, prob[1]->bBar, prob[1]->Cbar, soln->incumbX, fullObs-1);
            if (rhs == NULL) {
                errMsg("validating","bootStrapCheck","failed to compute subproblem right-hand side",0);
                return 1;
            }
            
            // Change the right-hand side in the solver
            status = changeRHS(cell->subprob->lp, prob[1]->num->rows, indices, rhs+1);
            if (status) {
                errMsg("solver", "bootStrapCheck", "Failed to change RHS of subproblem", 0);
                return 1;
            }
            
            // Solve the newly constructed Subproblem
            status = solveProblem(cell->subprob->lp, cell->subprob->name, cell->subprob->type, &status2);
            if (status) {
                errMsg("validating", "bootStrapCheck", "failed to solve checking subproblem",0);
                return 1;
            }
            
            // Fetch the objective value from solver
            obj = getObjective(cell->subprob->lp, PROB_LP);
            
            // Accumulative Construct smEst
            smEst += ((double) cell->omega->weight[obs] / (double) cell->k) * obj;
            
            mem_free(rhs);
        }

        printf("%d \t %d \t %d \t %.2f \t %.2f \t %.5f \n",
               cutID, cell->cuts->val[cutID]->cutObs, cell->cuts->val[cutID]->omegaCnt, Est, smEst, (smEst-Est)/smEst);


        if ( Est - smEst >= run.TOLERANCE ) {
            *check = FALSE;
            //return 0;
        } else {
            *check = TRUE;
        }
        
    }
    mem_free(fullObs);
    mem_free(indices);
    return 0;
}

/**********************************************************
 This is a small utility function used in bootStrapCheck.
 It samples a cut from all existing cuts.
 **********************************************************/
int sampleCut(int cutCnt) {
    int sampleID;
    sampleID = (int) ( cutCnt * randUniform(&run.EVAL_SEED) );
    if (sampleID < 1) {
        sampleID = 1;
    } else if (sampleID > cutCnt) {
        sampleID = cutCnt-1;
    }
    return sampleID;
}


/***********************************************************************
 * This function determines whether or not the current incumbent
 * solution is considered to be optimal.  In the early iterations
 * (less than MIN_ITER iterations), the function automatically assumes
 * that the incumbent is NOT optimal.  Provided that MIN_ITER iterations
 * have passed, it performs an optimality pre-test, to determine whether
 * a full optimality check is a worthwhile pursuit.  If so, it does
 * the full test, which involves reforming cuts and solving new master
 * problems.  The function returns TRUE if the incumbent is optimal;
 * FALSE otherwise.
 \***********************************************************************/
BOOL optimal(probType **prob, cellType *cell, solnType *soln) {
    
#ifdef OPTIMAL
    printf("\t~optimal()\n");
#endif
    
    BOOL optCheck = FALSE;
    clock_t tic, toc;
    
    tic = clock();
    
    if ( run.HORIZON == 1 ) { //Original SD code condition
        if (cell->k>run.MIN_ITER && soln->dualStableFlag == TRUE) {
            optCheck = TRUE;
        } else {
            optCheck = FALSE;
        }
    } else {
        if ( run.WARM_UP == 0 || (run.WARM_UP == 1 && run.PI_EVAL_CLEAR == 1) ){
            if (cell->k > run.MIN_ITER && soln->dualStableFlag == TRUE) {
                optCheck = TRUE;
            } else {
                optCheck = FALSE;
            }
        } else {
            if ( cell->t == 1 ) {
                if (cell->k > run.MIN_ITER && soln->dualStableFlag == TRUE) {
                    optCheck = TRUE;
                } else {
                    optCheck = FALSE;
                }
            } else {
                if ( cell->k > cell->Kc + run.PI_EVAL_DELAY && soln->dualStableFlag == TRUE ) {
                    optCheck = TRUE;
                } else {
                    optCheck = FALSE;
                }
            }
        }
    }
    
    
    if (optCheck == TRUE) {
        /* Record when dual is stablized */
        if (soln->dualStableIter != 0)
            soln->dualStableIter = cell->k;
        if ( preCheck(soln) ) {
            if (run.TEST_TYPE == 0) {
                return TRUE;
            }
            else if ( fullCheck(prob, cell, soln) )
            {
                printf(">\n");
                fflush(stdout);
                // Record time for each optimality test if conducted
                toc = clock();
                soln->runTime->optTime[soln->runTime->optCnt] = ((double) (toc-tic)) / CLOCKS_PER_SEC;
                soln->runTime->iterIndex[soln->runTime->optCnt] = cell->t;
                soln->runTime->optCnt++;
                return TRUE;
            }
            else
            {
                printf("<");
                fflush(stdout);
            }
        }
    }
    
    toc = clock();
    
    //Record time for each optimality test
    soln->runTime->optTime[soln->runTime->optCnt] = ((double) (toc-tic)) / CLOCKS_PER_SEC;
    soln->runTime->iterIndex[soln->runTime->optCnt] = cell->t;
    soln->runTime->optCnt++;

#ifdef DETAIL_CHECK
    printf("Finishing optimality checking completed...\n");
#endif
    
    return FALSE;
}//END of optimal()


/************************************************************************************************
 * A shorter/easier version of the optimality check is conducted before entering the full check.
 * This check asks whether the height at the candidate is close enough to the height at the 
 * incumbend to warrant an optimality test. 
 ************************************************************************************************/
BOOL preCheck(solnType *soln) {
    
#ifdef OPTIMAL
    printf("\t\t~pre-check()\n");
#endif
    
    // Check if candidate must be within some small percentage of incumbent cut
    if ( soln->candidEst > 0 ) {
        soln->preCheckFlag = (soln->candidEst > (1-run.PRE_EPSILON) * soln->incumbEst);
    } else {
        soln->preCheckFlag = (soln->candidEst > (1+run.PRE_EPSILON) * soln->incumbEst);
    }
    
    // If preCheck is ever passed, flip the flag of everFlag
    if (soln->preCheckFlag) {
    	soln->preCheckEverFlag = soln->preCheckFlag;
        
#ifdef OPTIMAL
        printf("\t\t :: candidEst(%f) > (1 - %f)IncumbEst(%f) \n",
                soln->candidEst, run.PRE_EPSILON, soln->incumbEst);
#endif
        
    }
    
    return soln->preCheckFlag;
}//END of preCheck()


/***********************************************************************\
 
     ----------- THIS SUBROUTINE IS NOT USED CURRENTLY -------------
 
 ** If we pass the first pre_test, we undertake a second pre_test to
 ** check to determine if the full test is worthwhile.  This function
 ** does bootstrapped based calculations to check whether or not the
 ** full test is likely to be passed.  If this second pretest fails, the
 ** full test would fail as well, so it need not be performed.
 **
 ** JH NOTE: Save bootstrap seed and re-initialize it if pre_test_2 passes
 
     ----------- THIS SUBROUTINE IS NOT USED CURRENTLY -------------
 
 \***********************************************************************/
BOOL preTestII(probType **prob, cellType *cell, solnType *soln) {
    
    cutType *T;
    double ULm;
    int *cdf, *observ;
    int m, sum;
    long long seed;
    
    int j;
    double ht, Sm;
    double eps_factor = 1.0;
    double p_pass_factor = 1.0;
    
#ifdef TRACE
    printf("Inside pre_test_2\n");
#endif
    
#ifdef RUN
    printf("\n Performing the second pre-test \n\n");
#endif
    
#ifdef OPT
    printf("omega->cnt=%d. sigma->cnt=%d. lambda->cnt=%d.\n",
           s->omega->cnt, c->sigma->cnt, c->lambda->cnt);
    print_num(p->num);
#endif
    
    seed = run.EVAL_SEED;
    sum = 0;
    T = chooseCuts(prob[0], cell, soln);
    if (!(observ = (intvec) arr_alloc(cell->k+1, int)) )
        errMsg("allocation", "preTestII", "Failed to allocate memory to observ", 0);
    if (!(cdf = arr_alloc(cell->omega->cnt+1, int)))
        errMsg("allocation", "preTestII", "Failed to allocate memroy to cdf", 0);
    
    empiricalDistrib(cell->omega, cdf);

    /* Find out how many of the resamplings satisfy optimality requirement */
    for (m = 0; m < run.M; m++) {
        sampleOmega(cdf, observ, cell->k);

        reformCuts(cell,prob[1]->num, prob[1]->coord, T, observ, cell->k);
#ifdef OPT
        printf("SS-PRINT: passed reform_cuts \n");
#endif
        
        /* Find the highest reformed cut at the incumbent solution */
        Sm = T->val[0]->alpha - vXv(T->val[0]->beta, soln->incumbX, NULL, prob[0]->num->cols);
        for (j = 1; j < T->cnt; j++)
        {
            ht = T->val[j]->alpha - vXv(T->val[j]->beta, soln->incumbX, NULL, prob[0]->num->cols);
            if (Sm < ht)
                Sm = ht;
        }
        Sm += vXv(prob[0]->dBar->val, soln->incumbX, NULL, prob[0]->num->cols);
        
        ULm = calculateULm(prob[0], cell, T, soln);
        
        /*
         printf("ULm=%lf, Sm=%lf inc = %lf Sm-ULm/inc = %lf \n",
         ULm, Sm, s->incumb_est, (Sm-ULm)/s->incumb_est);
         */
        
        /* JH: What value for epsilon here?  */
        if ((Sm - ULm) / soln->incumbEst <= eps_factor * run.PRE_EPSILON)
            sum++;
        
        /*  Skip out of the loop if there's no hope of meeting condition. - SS
         *** JH deleting for now 6/10/02 ***
         if (sum + sd_global->config.M - m - 1 < sd_global->config.PERCENT_PASS*sd_global->config.M) return FALSE;
         if (sum + sd_global->config.M - m - 1 < 0.75*sd_global->config.M)
         {      sd_global->config.MIN_ITER = c->k + 20;
         return (FALSE);
         }
         */
        
    }
    
#ifdef OPT
    printf("\nsum=%d\n",sum);
#endif
#ifdef RUN
    printf("\nc->k = %d, sum_pre2 = %d\n", c->k, sum);
#endif
    
    mem_free(cdf);
    mem_free(observ);
    freeCuts(T);
    
    /* JH: What value for percent_pass here?  */
    /*  return (sum >= p_pass_factor*sd_global->config.PERCENT_PASS * sd_global->config.M); */
    if (sum >= p_pass_factor * run.PERCENT_PASS * run.M)
    {
        run.EVAL_SEED = seed;
        return (TRUE);
    }
    /* this else shouldn't be executed, given SS's pass above */
    /*  jh whoopsie 6/10/02 
     else { if (sum <= 0.75*p_pass_factor*sd_global->config.PERCENT_PASS*sd_global->config.M)
     sd_global->config.MIN_ITER = c->k + 20;
     return (FALSE); 
     } 
     */
    else {
        return (FALSE); /* part of whoopsie */
    }
}

/***********************************************************************\
 ** This function performs a complete statistical test of optimality.
 ** First, it selects cuts whose height at the incumbent is "close" to
 ** incumbent cut's height.  Then, it performs M resamplings of the
 ** observations in omega, and reforms the selected cuts with respect
 ** to these observations (as if *they* were observed instead of the actual
 ** omega).  For each of the M resamplings, a master program containing the
 ** reformed cuts is solved, and if almost all of the solutions to these
 ** master programs
 \***********************************************************************/
BOOL fullCheck(probType **prob, cellType *cell, solnType *soln) {

#ifdef OPTIMAL
	printf("\t\t~fullCheck()\n");
#endif
    
    cutType *T;
    double  Sm, Lm, ht;
    intvec  cdf, observ;
    int     m, sum, i;
    int     numFailed = 0;  //Number of failed replications
    double  errSum;         //Sum of errors of failed replications
    clock_t tic, toc;

#ifdef OPTIMAL
    printf("\t\t :: Conducting full test at iteration %d...\n", cell->k);
#endif
        
    /* Step 0. Full Check Initialization */
    if ( !(observ = (intvec) arr_alloc(cell->k+1, int)) )
        errMsg("allocation","fullCheck","failed to allocate memory to observ",0);
    if (!(cdf = (intvec) arr_alloc(cell->omega->cnt+1, int)) )
        errMsg("allocation","fullCheck","failed to allocate memory to cdf",0);
    
    /* Step 1. Choose Cuts For Rebuilding */
    sum = 0;
    T = chooseCuts(prob[0], cell, soln);
    
    /* Step 2. Obtain empirical distribution of all observed residuals */
    empiricalDistrib(cell->omega, cdf);
    
    /* Step 3. Check on how many of the resamplings(Replication) satisfy optimality requirement */
    for ( m=0; m<run.M; m++) {
        
        tic = clock();
        /* Step 3.1 Resample from omegas for k observations */
        sampleOmega(cdf, observ, cell->k);
        toc = clock();
        soln->runTime->fullCheck_sample[soln->runTime->optCnt] += ((double) (toc-tic))/CLOCKS_PER_SEC;
        
        /* Step 3.2 Reconstruct cuts using the resampled observations */
        tic = clock();
        reformCuts(cell, prob[1]->num, prob[1]->coord, T, observ, cell->k);
        toc = clock();
        soln->runTime->fullCheck_reform[soln->runTime->optCnt] += ((double) (toc-tic))/CLOCKS_PER_SEC;
        
        /* Step 3.3 Find the highest reformed cut at the incumbent solution */
        Sm = T->val[0]->alpha - vXv(T->val[0]->beta,soln->incumbX, NULL, prob[0]->sp->mac);
        for (i=1; i<T->cnt; i++) {
            ht = T->val[i]->alpha - vXv(T->val[i]->beta,soln->incumbX, NULL, prob[0]->sp->mac);
            if (Sm < ht) {
                Sm = ht;
            }
        }

#ifdef OPTIMAL
        printf("\t\t :: Highest Reform Cut's Height = %f\n", Sm);
#endif
        
        /* Step 3.4 Calculate Lower Bound on Resampled Envrionment */
        /* If master is QP, do not include incumbdX x c in Sm */
        tic = clock();
        if (run.MASTER_TYPE == PROB_LP) {
            Sm += vXvSparse(soln->incumbX,prob[0]->dBar);
            Lm = solveTempMaster(prob[0], T, cell);
        } else {
            Lm = calTempLB(prob[0], cell, soln, T);
        }
        toc = clock();
        soln->runTime->fullCheck_lb[soln->runTime->optCnt] += ((double) (toc-tic))/CLOCKS_PER_SEC;

        
        /* Step 3.5 Checking of the Tests */
        if (DBL_ABS((Sm-Lm) / soln->incumbEst) <= run.EPSILON) {
            sum++;
        } else {
            /* Sum up errors of failed replication in full test. Only at Max_iter. */
            if ( cell->k >= cell->maxIter * run.HORIZON ) {
                numFailed++;
                errSum += (Sm-Lm) / soln->incumbEst;
            }
        }
    
#if defined(OPTIMAL)
        printf("\t\t :: k=%d m=%d, LB=%f, UB=%f, sum = %d, s->incumbEst = %f, (Sm-Lm)/incumbEst = %f\n",
                                cell->k, m+1, Lm, Sm, sum, soln->incumbEst, (Sm-Lm)/soln->incumbEst);
#endif
        
        /* If the flag for bootstrap test is disabled, then we simply claimed all tests passed the test */
        if ( m+1-sum >= (1-run.PERCENT_PASS) * run.M) {
            if ( cell->k >= run.MAX_ITER * run.HORIZON ) {
                soln->FTError = errSum / numFailed;
            }
            mem_free(observ);
            mem_free(cdf);
            freeCuts(T);
            return FALSE;
        }
    }

    /* Record the # of replication passed on the optimal solution */
    soln->repPassed = sum;
    
#ifdef OPTIMAL
    printf("\t\t :: Optimal Solution - %d replication passed....\n", soln->repPassed);
#endif
    
    mem_free(observ);
    mem_free(cdf);
    freeCuts(T);
    
    return (sum >= run.PERCENT_PASS * run.M);
    
}//END of fullCheck()


/***********************************************************************\
 ** This function selects all cuts whose height at the incumbent x
 ** is close to the height of the incumbent cut.  These cuts together
 ** are likely to provide good approximations of f at incumb_x, when
 ** they are reformed with new observations.  The function returns a
 ** new cut structure which contains room for cuts to be reformed.
 ** Only the _istar_ and _cut_obs_ fields of each cut have been initialized.
 \***********************************************************************/
cutType *chooseCuts(probType *prob, cellType *cell, solnType *soln) {

#ifdef OPTIMAL
    printf("\t\t\t~chooseCuts()\n");
#endif
    
    cutType *T;
    int cnt,i;
    
    T = newCuts(cell->maxCuts);
    
    /* Define what it means to be "close" to the incumbent */
    // double unknown; //Need to figure out...
    // unknown = vXvSparse(soln->incumbX, prob->dBar);
    
    /* Loop through the cuts available for reforming; pick close ones */
    for (cnt = 0; cnt < cell->cuts->cnt; cnt++)
    {
        if (soln->piM[cell->cuts->val[cnt]->rowNum + 1] > 0.00001)
        {
            T->val[T->cnt] = newCut(prob->num->cols, cell->cuts->val[cnt]->omegaCnt, cell->cuts->val[cnt]->cutObs, cell->t);
            //if (!(T->val[T->cnt]->iStar = (intvec) arr_alloc(cell->cuts->val[cnt]->omegaCnt+1, int)))
            //      errMsg("allocation", "chooseCuts", "one of the iStar", 0);
            for (i=0;i<cell->cuts->val[cnt]->omegaCnt;i++) //Just made change <= -> <
                T->val[T->cnt]->iStar[i] = cell->cuts->val[cnt]->iStar[i];
            T->val[T->cnt]->rowNum = cell->cuts->val[cnt]->rowNum;
#if defined(OPTIMAL)
            printf("\t\t\t :: Taking the cut: %d\n", cnt);
            printCut(T, prob->num, T->cnt);
#endif
            T->cnt++;
        }
    }
    
#if defined(OPTIMAL)
    printf("\t\t\t :: Selected %d cuts IN TOTAL...\n", T->cnt);
#endif
    
    return T;
}//End of chooseCuts()


/***********************************************************************\
 ** This function forms an empirical distribution on the observations
 ** stored in omega, and calculates an integer cdf to represent the
 ** distribution.  An observation which has been seen n times will have n
 ** times the probability of being chosen as an observation seen only once.
 \***********************************************************************/
void empiricalDistrib(omegaType *omega, int *cdf) {
    
#ifdef OPTIMAL
    printf("\t\t\t~empiricalDistrib()\n");
#endif
    
    int cnt;
    
    /* Calculate an integer cdf distribution for observations */
    /* If the cnt is not a valid omega idx, we know that weight[cnt] is 0 */
    cdf[0] = omega->weight[0];
    for (cnt = 1; cnt < omega->cnt; cnt++)
        cdf[cnt] = cdf[cnt - 1] + omega->weight[cnt];
    
}//END of empiricalDistrb()


/***********************************************************************\
 ** This function randomly selects a new set of observations from the
 ** old set of observations stored in omega.  Entries in omega which
 ** have been observed multiple times have a proportionally higher
 ** chance of being selected for the new set.  The function fills an
 ** array, assumed to be of a size equal to the number of iterations
 ** which have passed, with the new set of observations.
 ***********************************************************************/
void sampleOmega(int *cdf, int *observ, int k) {
    
#ifdef OPTIMAL
    printf("\t\t\t~sampleOmega()\n");
#endif
    
    int cnt, obs;
    int sample;
    
    /* Choose k observations according to cdf (k = number of iterations) */
    for (obs = 0; obs < k; obs++) {
        sample = randFun(k);
        for (cnt = 0; sample > cdf[cnt]; cnt++)
            /* Loop until sample falls below cdf */;
        observ[obs] = cnt;
    }
    
}//END of sampleOmega()

/**********************************************************************
 ** This function returns a uniform random number between [0, greatest-1]
 ** using our own random number generator.  This is not so good, since
 ** we are all running off the same seed... we ought to have different
 ** streams of random numbers.
 **********************************************************************/
int randFun(int greatest) {
    return (int) (randUniform(&run.EVAL_SEED) * greatest);
}//END of randFun()

/***********************************************************************\
 ** This function will calculate a new set of cuts based on the
 ** observations of omega passed in as _observ_, and the istar's
 ** which have already been stored in the _istar_ field of each cut.
 ** If an istar field does not exist for a given observation, then
 ** a value of zero is averaged into the calculation of alpha & beta.
 \***********************************************************************/
void reformCuts(cellType *cell, numType *num, coordType *coord, cutType *cut, intvec observ, int k) {

#ifdef OPTIMAL
    printf("\t\t\t~reformCuts()\n");
#endif
    
    int cnt, obs, idx, count;
    iType   istar;
    
    /* Loop through all the cuts and reform them */
    for (cnt = 0; cnt < cut->cnt; cnt++) {
#ifdef OPTIMAL
        printf("\t\t\t :: Reforming cut #%d\n", cnt);
#endif
        /* Begin with re-initializing cut coefficients of zero */
        for (idx = 0; idx <= num->prevCols; idx++)
            cut->val[cnt]->beta[idx] = 0.0;
        cut->val[cnt]->alpha = 0.0;
        
        count = 0;
        /* Reform this cut based on resampled observations */
        for (obs = 0; obs < k; obs++) {
            
            /* Only compute if this observation actually exists in omega */
            if (validOmegaIdx(cell->omega, observ[obs])) {
                
                /* Only sum values if the cut has an istar for this observation */
                if ( (observ[obs] < cut->val[cnt]->omegaCnt) ){
                    
                    istar.sigma = cut->val[cnt]->iStar[observ[obs]];
                    istar.delta = cell->sigma->lambdaIdx[istar.sigma];
                    
                    cut->val[cnt]->alpha += cell->sigma->vals[istar.sigma].b + cell->delta->val[istar.delta][observ[obs]].b;
                    
                    // This part is modified with consideration when time series is used for generating omega
                    if ( run.HORIZON > 1 ) {
                        cut->val[cnt]->alpha += cell->gamma->endo->pibARMA[cell->sigma->lambdaIdx[istar.sigma]] * 1.0;
                        cut->val[cnt]->alpha += cell->auxDelta->val[istar.delta][observ[obs]].b * 1.0;
                        cut->val[cnt]->alpha += cell->trimDelta->val[istar.delta][observ[obs]].b * 1.0;
                    }
                    
                    for (idx = 1; idx <= num->cntCcols; idx++)
                        cut->val[cnt]->beta[coord->colsC[idx]] += cell->sigma->vals[istar.sigma].C[idx] * 1.0;
                    for (idx = 1; idx <= num->rvColCnt; idx++)
                        cut->val[cnt]->beta[coord->rvCols[idx]] += cell->delta->val[istar.delta][observ[obs]].C[idx] * 1.0;
                
                    count++;
                }
            }
            
        }
        
        /* Take the average of the alpha and beta values */
        for (idx = 0; idx <= num->prevCols; idx++) {
            cut->val[cnt]->beta[idx] /= (double) k;
        }
        cut->val[cnt]->alpha /= (double) k;
        if ( cell->lbType == NONTRIVIAL ) {
            cut->val[cnt]->alpha += (1 - (double) count / (double) k) * cell->lb; //run.Eta0
        }

#ifdef OPTIMAL
        printCut(cut, num, cnt);
#endif
        
    }
}//END of reformCuts()

/***********************************************************************\
 ** This function returns TRUE if the omega->idx array is storing an
 ** realization vector at the index passed to it.  It return FALSE
 ** if the index is out-of-bounds or does not contain a realization.
 \***********************************************************************/
BOOL validOmegaIdx(omegaType *omega, int idx) {
    if (idx < 0 || idx >= omega->cnt)
        return FALSE;
    
    /* Check the range of valid indices */
    if (idx < 0 || idx >= omega->cnt)
        return FALSE;
    
    /* Check the contents of the idx array */
    if (omega->vals[idx] == NULL)
        return FALSE;
    
    return TRUE;
}


/***********************************************************************\
 ** This function loads, solves, and frees a master program with cuts
 ** specified by _T_.  It returns the value of the objective function
 ** at the optimal solution.
 \***********************************************************************/
double solveTempMaster(probType *prob, cutType *cuts, cellType *cell) {
    
    oneProblem *tempMaster;
    double ans;
    int status, stat1;
    
    tempMaster = newStoMaster(prob->sp, NULL, cuts, 0);
    
    status = solveProblem(tempMaster->lp, tempMaster->name, tempMaster->type, &stat1);
    
    ans = getObjective(tempMaster->lp, tempMaster->type);
    
    /* Free lp pointer */
    freeProblem(tempMaster->lp);
    
    freeOneProblem(tempMaster);
    
    return ans;
}



/****************************************************************************\
 This function is to calculate the lower bound on the optimal value which
 is used in stopping rule in full_test() in optimal.c in the case of
 regularized approach.
\****************************************************************************/
double calTempLB(probType *prob, cellType *cell, solnType *soln, cutType *T) {

#ifdef OPTIMAL
    printf("\t\t\t~calTempLB();\n");
#endif
    
    vector bk;         /* vector: b - A*incumb_x. */
    vector lambda;     /* vector: the dual of the primal constraints. */
    double bk_lambda;   /* scalar: bk*lambda. */
    
    sparseMatrix *ATrans;   /* sparse_matrix: the transpose of A. */
    vector ATransLambda;   /* vector: - A_Trans * lambda. */
    
    vector theta;      /* the dual of the reformed cut constraints. */
    vector Vk;         /* the vector of scalars of cut constraints. Vk = alpha - beta * incumb_x  */
    double Vk_theta;   /* scalar: Vk*theta. */
    vector Bk_theta;   /* vector: Bk_Transpose * theta, where Bk_Transpose is the matrix of cut coefficients. */
    vector Bk_col;     /* vector: Store one column of Bk while calculating Bk_theta. */
    vector q_vec;      /* vector: c + Bk_theta - A_Trans_lambda. */
    double q_term;     /* scalar: q_vec * q_vec. */
    vector lb;         /* vector: store the lower bounds.  */
    vector ub;         /* vector: store the upper bounds.  */
    
    double Lm;         /* The calculated lower bound of the optimal value. */
    
    /* double eta_zero; */
    int cnt, i;
    
    /* Memory Allocation*/
    if (!(bk = (vector) arr_alloc(prob->num->rows + cell->fCutAdded->cnt + 1, double)))
        errMsg("Allocation", "cal_temp_lb", "bk",0);
    if (!(lambda = (vector) arr_alloc(prob->num->rows + cell->fCutAdded->cnt + 1, double)))
        errMsg("Allocation", "cal_temp_lb", "lambda",0);
    if (!(ATransLambda = (vector) arr_alloc(prob->num->cols+1, double)))
        errMsg("Allocation", "cal_temp_lb", "A_lambda",0);
    if (!(ATrans = (sparseMatrix *) mem_malloc(sizeof(sparseMatrix))))
        errMsg("Allocation", "cal_temp_lb", "A_Trans",0);
        if (!(ATrans->val = (vector) arr_alloc(prob->Dbar->cnt+1, double)))
            errMsg("Allocation", "cal_temp_lb", "A_Trans->val",0);
        if (!(ATrans->row = (intvec) arr_alloc(prob->Dbar->cnt+1, int)))
            errMsg("Allocation", "cal_temp_lb", "A_Trans->row",0);
        if (!(ATrans->col = (intvec) arr_alloc(prob->Dbar->cnt+1, int)))
            errMsg("Allocation", "cal_temp_lb", "A_Trans->col",0);
    if (!(theta = (vector) arr_alloc(T->cnt+1, double)))
        errMsg("Allocation", "cal_temp_lb", "theta",0);
    if (!(Vk = (vector) arr_alloc(T->cnt+1, double)))
        errMsg("Allocation", "cal_temp_lb", "Vk",0);
    if (!(Bk_theta = (vector) arr_alloc(prob->num->cols+1, double)))
        errMsg("Allocation", "cal_temp_lb", "Bk_theta",0);
    if (!(Bk_col = (vector) arr_alloc(T->cnt+1, double)))
        errMsg("Allocation", "cal_temp_lb", "Bk_col",0);
    if (!(q_vec = (vector) arr_alloc(prob->num->cols+1, double)))
        errMsg("Allocation", "cal_temp_lb", "q_vec",0);
    if (!(lb = (vector) arr_alloc(prob->num->cols+1, double)))
        errMsg("Allocation", "cal_temp_lb", "lb",0);
    if (!(ub = (vector) arr_alloc(prob->num->cols+1, double)))
        errMsg("Allocation", "cal_temp_lb", "ub",0);
    
    
    /*************************** Calculate upper bound, Sm. ***************************/
    //Don't consider feasbility cuts in rolling horizon problem

    /**********************************************************************************/
    
    
    
    /*************************** Calculate lower bound, Lm. ***************************/
    
    // Before we actually get the lower bound, we need to prepare a bunch of variables ready
    /* --> Original Part means the probelm inherite rows in master problem
     *
     * 1. bk = b - A x incumbX, composed using original part and cuts part
     * 2. lambda = (-) soln->piM, composed using original part and cuts part
     # 3. bk_lambda = bk x lambda
     * 4. A^\top
     * 4. A^\top x lambda
     * 5.
     */
    
    /*------- 1. Compute b_k = b - A x incumbX -------*/
    for (cnt = 0; cnt < prob->num->rows; cnt++) {
        bk[cnt + 1] = cell->master->rhsx[cnt];
    }
    MSparsexvSub(prob->Dbar, soln->incumbX, bk);
    for (cnt = 0; cnt < cell->fCutAdded->cnt; cnt++) {
        bk[prob->num->rows + cnt + 1] = cell->fCutAdded->val[cnt]->alpha
                - vXv(cell->fCutAdded->val[cnt]->beta, soln->incumbX, NULL, prob->num->cols);
    }
    /*------- 2. Compute lambda = (-soln->piM) -------*/
    /* Dual values' sign need to be flipped here before assigning to lambda*/
    /* Obtain lambda from s->Master_pi of original master problem constraints. */
    for (cnt = 0; cnt < prob->num->rows; cnt++) {
        lambda[cnt + 1] = - soln->piM[cnt + 1];
    }
    
    /* Obtain lambda from s->Master_pi of master problem feasibility cuts. */
    for (cnt = 0; cnt < cell->fCutAdded->cnt; cnt++) {
        lambda[prob->num->rows + cnt + 1] = - soln->piM[cell->fCutAdded->val[cnt]->rowNum + 1];
    }
    
    /*------- 3. Compute bk * lambda -------*/
    bk_lambda = vXv(bk, lambda, NULL, prob->num->rows + cell->fCutAdded->cnt);
    
    /*------- 4. Compute A^\top -------*/
    ATrans->cnt = prob->Dbar->cnt;
    for (cnt = 1; cnt <= ATrans->cnt; cnt++) {
        ATrans->val[cnt] = prob->Dbar->val[cnt];
        ATrans->row[cnt] = prob->Dbar->col[cnt];
        ATrans->col[cnt] = prob->Dbar->row[cnt];
    }
    
    /*-------- 5. Compute ATransLambda = A^\top * lambda -------*/
    MSparsexvSub(ATrans, lambda, ATransLambda);
    for (i = 0; i < prob->num->cols; i++) {
        for (cnt = 0; cnt < cell->fCutAdded->cnt; cnt++) {
            ATransLambda[i + 1] -= cell->fCutAdded->val[cnt]->beta[i+1] * lambda[prob->num->rows + cnt + 1];
        }
    }
    // Adding dual slacks back...
    for (i = 0; i < prob->num->cols; i++) {
        ATransLambda[i + 1] += soln->djM[i + 1];
    }
    
    /* ----- */
    for (cnt = 0; cnt < T->cnt; cnt++) {
        /*
         ** Obtain theta from s->Master_pi. Obtain Vk from
         **       T->val[i]->alpha - T->val[i]->beta * incumb_x.
         ** Note: be careful of the corresponding row number of each cut. And also
         ** be careful that the first element of s->Master_pi was reserved for one
         ** norm.
         */
        /* 6. Compute theta = the dual of the reformed cut constraints */
        theta[cnt + 1] = ((double) cell->k / (double) T->val[cnt]->cutObs) * soln->piM[T->val[cnt]->rowNum + 1];
        /* 7. Compute Vk = alpha - beta x incumbX */
        Vk[cnt + 1] = T->val[cnt]->alpha - vXv(T->val[cnt]->beta, soln->incumbX, NULL, prob->num->cols);
    }

    /* 8. Compute Vk_theta = Vk x theta */
    Vk_theta = vXv(Vk, theta, NULL, T->cnt);
    
    for (i = 1; i <= prob->num->cols; i++) {
        /* 9. Compute each column Bk_col = beta (reformed cuts) */
        for (cnt = 0; cnt < T->cnt; cnt++)
            Bk_col[cnt + 1] = T->val[cnt]->beta[i];
        /* 10. Compute each column Bk_theta = beta x theta */
        Bk_theta[i] = vXv(Bk_col, theta, NULL, T->cnt);
    }
    
    /* 11. Compute q_vec = Bk_theta - ATransLambda */
    for (i = 1; i <= prob->num->cols; i++)
        q_vec[i] = prob->dBar->val[i] - Bk_theta[i] - ATransLambda[i];
    
    /* 12. Compute q_vec x q_vec */
    q_term = vXv(q_vec, q_vec, NULL, prob->num->cols);

    /**********************************************************
     ** Now it is the time to calculate the lower bound, Lm. **
     **********************************************************/
    
    /* 13. Compute the lower bound */
    Lm = Vk_theta + bk_lambda - q_term / cell->quadScalar / 2.0;
    
#if defined(OPT)
    printf("\n\t\t\t :: LB[%.2f] = Vk_theta[%.2f] + bk_lambda[%.5f] - q_term[%.5f] / quadScalar[%f] / 2\n"
           ,Lm, Vk_theta, bk_lambda, q_term, cell->quadScalar);
#endif
    
    /**********************************************************************************/

    mem_free(bk);
    mem_free(lambda);
    mem_free(ATransLambda);
    mem_free(ATrans->col);
    mem_free(ATrans->row);
    mem_free(ATrans->val);
    mem_free(ATrans);
    mem_free(theta);
    mem_free(Vk);
    mem_free(Bk_theta);
    mem_free(Bk_col);
    mem_free(q_vec);
    mem_free(lb);
    mem_free(ub);
    return Lm;
    
}

/*********************************************************\
 ** This function calculates an upper bound on a lower
 ** bound on the value of the resampled master problem,
 ** using the dual solution ** to the actual master problem
 ** and an exact penalty function.
 ** JH 5/98  WORK
 \*********************************************************/
double calculateULm(probType *prob, cellType *cell, cutType *T, solnType *soln) {
    int j, t;
    double zbar, dbar, value;

    writeProblem(cell->subprob->lp, "OPTsubprob.lp"); /* 2011.10.30 */
    
    zbar = soln->candidEst;
    for (t = 0; t < cell->cuts->cnt; ++t) {
        zbar -= cell->cuts->val[t]->alpha
        * soln->piM[cell->cuts->val[t]->rowNum + 1];
    }
    
    value = zbar;
    for (t = 0; t < T->cnt; ++t){
        value += T->val[t]->alpha * soln->piM[T->val[t]->rowNum + 1];
    }
    
    for (j = 1; j <= prob->num->cols; ++j){
        dbar = soln->djM[j];
        for (t = 0; t < cell->cuts->cnt; ++t){
            dbar += cell->cuts->val[t]->beta[j] * soln->piM[cell->cuts->val[t]->rowNum + 1];
        }
        
        for (t = 0; t < T->cnt; ++t){
            dbar -= T->val[t]->beta[j] * soln->piM[T->val[t]->rowNum + 1];
        }

        value += dbar * soln->candidX[j];
        
    }
    return (value);
}


