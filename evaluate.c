//
//  evaluate.c
//  rollingSP
//
//  Created by Site Wang on 4/5/16.
//  Copyright © 2016 Site Wang. All rights reserved.
//
#include <time.h>
#include "evaluate.h"
#include "arma.h"

extern runType run;

/*****************************************************************************************************************
 * This subroutines is used to evaluate the optimal solution with the expectation value. The optimal incumbent
 * solution is used as the optimality have been declared. With that solution, we re-sample the stochastic process
 * until the expectation of that solution is stablized. One should expect to see that the estimate objective
 * on the incumbent solution falls into the confidence interval of the expectation objective value.
 *****************************************************************************************************************/
void evalStoOpt(stocType *stoc, probType **prob, cellType *cell, solnType *soln, FILE *SubOptVec, FILE *ePtr) {
    
    
    vector 	fullObs = NULL, Subsol = NULL, nullObs = NULL, rhs = NULL;
    intvec  indices;
    double 	cx, obj, mean, variance, stdev, temp, CI[2], totalTime;
    int		i, status = 0, cnt, stat2;
    clock_t	tic, toc;
    
#ifdef EVAL
    printf("~evalOpt()\n");
#endif
    
    if ( !(fullObs = (vector) arr_alloc(prob[1]->num->numRV, double)) )
        errMsg("allocation", "evaluateOpt", "observ", 0);
    if ( !(nullObs = (vector) arr_alloc(prob[1]->num->numRV, double)) )
        errMsg("allocation", "evaluateOpt", "nullObs", 0);
    if ( !(Subsol = (vector) arr_alloc(prob[1]->num->cols+1, double)) )
        errMsg("allocation", "evaluateOpt", "Subsol", 0);
    if ( !(indices = (intvec) arr_alloc(prob[1]->num->rows, int)) )
        errMsg("validating", "solveSubporb", "ind", 0);
    for(i = 0; i < prob[1]->sp->mar; i++) indices[i] = i;
    
    /* first stage cost */
    cx = vXvSparse(soln->incumbX, prob[0]->dBar);
    
#ifdef EVAL
    printf(" :: Evaluating optimal solution...\n");
#endif

    tic = clock();
    cnt = 0; mean = 0; variance = 0.0; stdev = INFBOUND; //Initialization
    
    while (3.92 * stdev > run.EVAL_ERROR * DBL_ABS(mean) || cnt < run.EVAL_MIN_ITER ) {
        
        if (run.HORIZON > 1) {
            generateARIMA(cell->t, cell, fullObs, NULL, NULL, NULL, NULL, &run.RESAMPLE_SEED);
        } else {
            generateOmega(cell, stoc, nullObs, fullObs, nullObs, NULL, NULL, &run.RESAMPLE_SEED);
        }
        
#ifdef EVAL
        printVector(fullObs-1, stoc->numOmega, NULL);
#endif
        
        // Compute RHS Value of the right-hand side value
        rhs = computeRHS(prob[1]->num, prob[1]->coord, prob[1]->bBar, prob[1]->Cbar, soln->incumbX, fullObs-1);
        if (rhs == NULL) {
            errMsg("validating","bootStrapCheck","failed to compute subproblem right-hand side",0);
        }
        
        // Change the right-hand side in the solver
        status = changeRHS(cell->subprob->lp, prob[1]->num->rows, indices, rhs+1);
        if (status) {
            errMsg("solver", "bootStrapCheck", "Failed to change RHS of subproblem", 0);
        }
        mem_free(rhs);
        
        status = solveProblem(cell->subprob->lp, cell->subprob->name, cell->subprob->type, &stat2);
        if ( status ) {
            if ( stat2 == STAT_INFEASIBLE ) {
                /* subproblem is infeasible add a feasibility cut */
                printf("Warning:: Subproblem is infeasible: need to create feasibility cut.\n");
                exit(1);
            }
            else {
                errMsg("algorithm", "evaluateOpt", "failed to solve subproblem in solver", 0);
                exit(1);
            }
        }
        
        status= getPrimal(cell->subprob->lp, Subsol, prob[1]->num->cols);
        if ( status ) {
            errMsg("algorithm", "evaluateOpt", "failed to get the primal vector",0);
            exit(1);
        }
        
        if ( status ) {
            errMsg("algorithm", "evaluateOpt", "failed to write the primal vector",0);
            exit(1);
        }
        
        /* use subproblem objective and compute evaluation statistics */
        obj = getObjective(cell->subprob->lp, PROB_LP);
        if ( cnt == 0 )
            mean = obj;
        else {
            temp = mean;
            mean = mean + (obj - mean) / (double) (cnt + 1);
            variance  = (1 - 1 / (double) cnt) * variance + (cnt + 1) * (mean - temp) * (mean - temp);
            stdev = sqrt(variance/ (double) cnt);
        }
        
        cnt++;
        
        /* Print the results every once in a while for long runs */
        if (!(cnt % 100)) {
            printf(".");
            fflush(stdout);
        }
        if (!(cnt % 10000))
            printf("\n\nobs:%d; mean:%lf; error: %lf \n 0.99 CI: [%lf , %lf]\n", cnt, cx + mean, 3.92 * stdev / mean,
                   cx + mean - 2.576 * stdev, cx + mean + 2.576 * stdev);
    }
    toc = clock();
    
    /* Record Evaluation Time */
    soln->runTime->evalOptTime[cell->t] = ((double) (toc-tic)) / CLOCKS_PER_SEC;
    
    totalTime = (toc-tic)/CLOCKS_PER_SEC;
    mean += cx;
    CI[0] = mean - 2.576 * stdev;
    CI[1] = mean + 2.576 * stdev;
    
    /* Print the value of the solution to a file and the screen */
    printf("\nFinal Estimate :: period %d :: iter %d :: obs:%d, mean:%lf, stdev:%lf, 0.99 C.I.: [%lf , %lf] \n\n", cell->t, cell->k, cnt, mean, stdev, CI[0], CI[1]);
    
    if (cell->t == 1) {
        fprintf(ePtr, "t \t k \t obsCnt \t mean \t stdev \t CI_L \t SD_EST \t CI_HI \t totalTime \n");
    }
    fprintf(ePtr, "%d \t %d \t %d \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \n",
            cell->t, cell->k, cnt, mean, stdev, CI[0], soln->incumbEst, CI[1], totalTime);
    
    mem_free(fullObs);
    mem_free(Subsol);
    mem_free(nullObs);
    mem_free(indices);
    
}//END evaluateOpt()


/*********************************************************************************************************************
 * This subroutine is designed for evaluating deterministic solution vs. expectation solution. Similar to evalOpt,
 * the solution obatined by sampling the stochastic process is evaluated on its variance. Once we hit minimum eval
 * iterations AND get the variance stablized. It is possible for us to conclude the expectation solution is obtained.
 * Then we construct a confidence interval for our statistically obtained expectation model as a reference interval.
 * Difference from the stochastic version. The procedure here needs to change a few indices since the structure of
 * the problem is decomposed differently when solving it deterministically.
 *********************************************************************************************************************/
void evalDetOpt(oneProblem *orig, timeType *tim, stocType *stoc, probType **prob, cellType *cell, solnType *soln, FILE *SubOptVec, FILE *ePtr) {
    
#ifdef EVAL
    printf("~evalDetOpt()\n");
#endif
    
    probType **stageProb; // Re-decompose the problem for
    
    vector 	fullObs, rhs = NULL;
    intvec  indices;
    double 	cx, obj, mean, variance, stdev, temp, CI[2], totalTime;
    int		i, status = 0, cnt, stat2;
    clock_t	tic, toc;
    
    tim->numStages++;   //Restore the stages cnt
    stageProb = newProb(orig, stoc, tim, 0, 0.1); //Assumed 0 lower bound
    if (prob == NULL) {
        errMsg("decompose","setup_rSP","Failed to decompose the problem",0);
        goto EVAL_FAIL;
    }
    
    if ( !(fullObs = (vector) arr_alloc(stageProb[1]->num->numRV, double)) )
        errMsg("allocation", "evaluateOpt", "observ", 0);
    if ( !(indices = (intvec) arr_alloc(stageProb[1]->num->rows, int)) )
        errMsg("validating", "solveSubporb", "ind", 0);
    for(i = 0; i < stageProb[1]->sp->mar; i++) indices[i] = i;
    
    /* first stage cost: there is no firsts stage when solving deterministically */
    cx = vXvSparse(soln->incumbX, stageProb[0]->dBar);
    
#ifdef EVAL
    printf(" :: Evaluating deterministic optimal solution -> cx = %lf\n", cx);
#endif
    
    //Initialization
    tic = clock();
    cnt = 0;
    mean = 0;
    variance = 0.0;
    stdev = INFBOUND;
    
    // Update stageProb subprob structure for generating subproblem
    for (i=1; i<=stageProb[1]->bBar->cnt; i++)
        stageProb[1]->bBar->val[i] = prob[0]->bBar->val[stageProb[0]->num->rows+i];
    
    // Wipe out the stochastic rows bBar to be zero
    for (i=1; i<=stageProb[1]->num->numRV; i++)
        stageProb[1]->bBar->val[stageProb[1]->coord->omegaRow[i]] = 0.0;
    
    while (3.92 * stdev > run.EVAL_ERROR * DBL_ABS(mean) || cnt < run.EVAL_MIN_ITER ) {
        
        /* Simulate a solution */
        generateARIMA(cell->t, cell, fullObs, NULL, NULL, NULL, NULL, &run.RESAMPLE_SEED);
        
        rhs = computeRHS(stageProb[1]->num, stageProb[1]->coord, stageProb[1]->bBar, stageProb[1]->Cbar, soln->incumbX, fullObs-1);
        if (rhs == NULL) {
            errMsg("algorithm", "evalDetOpt", "failed to compute subproblem right-hand side", 0);
            goto EVAL_FAIL;
        }
        
        // Change the right-hand side in the solver
        status = changeRHS(stageProb[1]->sp->lp, stageProb[1]->num->rows, indices, rhs+1);
        if (status) {
            errMsg("solver", "bootStrapCheck", "Failed to change RHS of subproblem", 0);
            goto EVAL_FAIL;
        }
        
        status = solveProblem(stageProb[1]->sp->lp, stageProb[1]->sp->name, stageProb[1]->sp->type, &stat2);
        if ( status ) {
            if ( stat2 == STAT_INFEASIBLE ) {
                /* subproblem is infeasible add a feasibility cut */
                printf("Warning:: Subproblem is infeasible: need to create feasibility cut.\n");
                goto EVAL_FAIL;
            }
            else {
                errMsg("algorithm", "evaluateOpt", "failed to solve subproblem in solver", 0);
                goto EVAL_FAIL;
            }
        }
        
        /* use subproblem objective and compute evaluation statistics */
        obj = getObjective(stageProb[1]->sp->lp, PROB_LP);
#ifdef EVAL
        printf(" :: Subproblem Objective -> %lf\n", obj);
#endif
        if ( cnt == 0 )
            mean = obj;
        else {
            temp = mean;
            mean = mean + (obj - mean) / (double) (cnt + 1);
            variance  = (1 - 1 / (double) cnt) * variance + (cnt + 1) * (mean - temp) * (mean - temp);
            stdev = sqrt(variance/ (double) cnt);
        }
        cnt++;
        /* Print the results every once in a while for long runs */
        if (!(cnt % 100)) {
            printf(".");
            fflush(stdout);
        }
        if (!(cnt % 10000))
            printf("\n\nobs:%d mean:%lf  error: %lf \n 0.95 CI: [%lf , %lf]\n", cnt, cx + mean, 3.92 * stdev / mean,
                   cx + mean - 2.576 * stdev, cx + mean + 2.576 * stdev);
        
        mem_free(rhs);
    }
    toc = clock();
    
    /* Record Evaluation Time */
    soln->runTime->evalOptTime[cell->t] = ((double) (toc-tic)) / CLOCKS_PER_SEC;
    
    totalTime = (toc-tic)/CLOCKS_PER_SEC;
    mean += cx;
    CI[0] = mean - 2.576 * stdev;
    CI[1] = mean + 2.576 * stdev;
    
    /* Print the value of the solution to a file and the screen */
    printf("\nFinal Estimate :: period %d :: iter %d :: obs:%d, mean:%lf, 0.99 C.I.: [%lf , %lf] \n\n", cell->t, cell->k, cnt, mean, CI[0], CI[1]);
    
    if (cell->t == 1) {
        fprintf(ePtr, "t \t k \t cnt(obs) \t mean \t stdev \t CI_L \t DET_EST \t CI_HI \t totalTime \t dualStable \n");
    }
    fprintf(ePtr, "%d \t %d \t %d \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %d \n",
            cell->t, 1, cnt, mean, stdev, CI[0], soln->optValM, CI[1], totalTime, soln->dualStableIter);
    
    
EVAL_FAIL:
    if (indices) mem_free(indices);
    if (fullObs) mem_free(fullObs);
    freeProbType(stageProb, tim->numStages);
    tim->numStages--;
    
}//END evaluateOpt()


/************************************************************************************************************************* 
 * This subroutine is used to evaluate the dual stability when generating a cut. It will evaluate the cut height produced
 * using the first 90% iterations and compare it with the real cut height. If the difference is pretty small and getting
 * stablized, that indicates the dual solution we generated before is being utilized to generate supporting hyperplanes.
 * Hence, we can claim the dual stability which give us the a green light for doing full check. 
 *************************************************************************************************************************/
void evalDualStability(vector piSRatio, double argmaxPreSum, double argmaxAllSum, BOOL *dualStableFlag, int numSamples) {
    
    double  variance = 1.0;
    
    /* Calculate Correponding position's subproblem dual solution ratio */
    piSRatio[numSamples % run.SCAN_LEN] = argmaxPreSum / argmaxAllSum;
    if ( numSamples - run.PI_EVAL_START > run.SCAN_LEN )
        variance = calcVariance(piSRatio, NULL, NULL, 0);
    if ( (DBL_ABS(variance) >= run.DUAL_EPSILON) || (piSRatio[numSamples % run.SCAN_LEN] < run.PERCENT_PASS) ) {
        *dualStableFlag = FALSE;
    } else {
        *dualStableFlag = TRUE;
    }
}


/*********************************************************
 ** This function calculate the variance of the vector x.
 *********************************************************/
double calcVariance(vector x, vector mean_value, vector stdev_value, int batch_size) {
    double mean, vari, temp;
    int count, length;
    double stdev;
    stdev = 10000000.0;
    temp = 0.0;
    mean = x[0];
    vari = 0.0;
    
    if (mean_value != NULL) {
        length = batch_size;
    } else {
        length = run.SCAN_LEN;
    }
    
    for (count = 1; count < length; count++) {
        temp = mean;
        mean = mean + (x[count] - mean) / (double) (count + 1);
        vari = (1 - 1 / (double) count) * vari
        + (count + 1) * (mean - temp) * (mean - temp);
    }
    
    if (mean_value != NULL) {
        *mean_value = mean;
    }
    if (stdev_value != NULL) {
        stdev = sqrt(vari / (double) count);
        *stdev_value = stdev;
    }
    
    return vari;
    
}
