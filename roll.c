//
//  roll.c
//  rollingSP
//
//

#include "roll.h"
#include "output.h"
#include "optimal.h"
#include "evaluate.h"

extern runType  run;
extern BOOL     newPeriod;
double *base;

/***********************
 * This subroutine is designed for setting up a next-period-deterministic problem.
 When solving a problem deterministically,
 * we need to update the following:
 * 1. Dynamic primal solution to first stage right hand side;
 * 2. Dynamic subproblem right hand side;
 * 3. Stochastic Right Hand Side for the suproblem;
 * To do the thrid step, we need to perform the same procedure as in stochastic version where we warm up the time
 * series model pointers to indicate where we are actually standing at among the data stream.
 * 
 * In deterministic version, the problem is not segemented into master problem and subproblem. It is considered to
 * be one big problem. Hence, when warming up the problem's dynamic parts, we need to be careful with the how to index
 * the coordination of random rows and dynamic rows. These things are taken care of within warmUpDynamic(). Both
 * deterministic version and stochastic version calls the same subroutine here.
 *************/
int setupDetPeriod(oneProblem *detProblem, cellType *cell, probType **prob, solnType *soln, BOOL *initialPeriod) {
    
#ifdef ROLL
    printf("~detWarmUp()\n");
#endif
    
    intvec indices;
    vector rhs;
    int current;
    int status = 0, rows;
    int i, v, h, rowLoc;
    double tsHistory, newRHS;
    
    /********** Step 1. Parameter/Pointer Adjustment **********/
    if (*initialPeriod == FALSE) {
        newPeriod = TRUE;
        cell->t++;
        cell->arima->history->open += cell->arima->subh;            // More data is revealed
        cell->arima->history->past += cell->arima->subh;
#ifdef  ROLL
        printf(" :: Adjusting stream pointer open = %d; (VAR1) history[open] = %lf; noise[open] = %lf\n",
               cell->arima->history->open, cell->arima->history->obs[cell->arima->history->open][1],
               cell->arima->noise->obs[cell->arima->noise->open][1]);
#endif
    } else {
        // Initial Period Explicit Operations :: Double Check -> Restric cell->t to 1
        cell->t = 1;
    }
    
    current = cell->arima->history->open;
    
    /* Update the time series model for data strem pointers to move */
    status = warmUpTimeSeries(cell->arima, current);
    if (status) {
        errMsg("warming", "warmUp", "Failed to warm up the stochastic parts", 0);
        return 1;
    }
    
    /* Total rows of the deterministic problem */
    rows = getNumRows(detProblem->lp);
    
    /* Get Problem's Right Hand Side for Updating the bBar */
    if ( !(rhs = (vector) arr_alloc(rows, double)) )
        errMsg("allocation", "detWarmUp", "rhs", 0);
    if ( getRhsx(detProblem->lp, 0, rows, rhs) ) {
        errMsg("solver","warmUpDynamic","Failed to reach current master problem right hand side",0);
        return 1;
    }
    
    /* Update problem's right hand side with data stored in time series model */
    for (v=1; v<=cell->arima->var; v++) {
        for (h=1; h<=cell->arima->variate[v]->cnt; h++) {
            tsHistory = cell->arima->history->obs[current+cell->arima->variate[v]->seq[h]][v];
            newRHS = reverse(tsHistory, v, cell->arima->variate[v]->seq[h], cell->t, cell->arima);
            rowLoc = cell->arima->variate[v]->row[h] - 1;
            //Note that tsType saves the index that works with cellTypem which starts from 1
#ifdef ROLL
            printf(" :: Updating rowLoc = %d, oldRHS = %f with newRHS = %f\n", rowLoc, rhs[rowLoc], newRHS);
#endif
            // Everything is ready, update the right hand side
            rhs[rowLoc] = max(0,newRHS * run.SCALAR);
        }
    }
    
    if ( !(indices = (intvec) arr_alloc(rows, int)) )
        errMsg("allocation","warmUpCell","failed to allocate memory to indices",0);
    for (i=0; i<rows; i++) indices[i] = i;
    
    // Send Ready rhs back to master problem
    status = changeRHS(detProblem->lp, rows, indices, rhs);
    if (status) {
        errMsg("warmup","warmUpCell","Failed to update the warm-up right hand side",0);
        return 1;
    }

#ifdef ROLL
    writeProblem(detProblem->lp, "afterTSWarmStart.lp");
#endif
    
    /* Update all the dynamic parts in deterministic problem */
    status = warmUpDynamic(prob, soln, cell, detProblem, NULL);
    if (status) {
        errMsg("warmUp", "detWarmUp", "Failed to warmup the problem", 0);
        return 1;
    }
    
    *initialPeriod = FALSE;
    
    mem_free(indices); mem_free(rhs);
    
    return 0;
}//END of detWarmUp()

/*************************************************************************************************************
 -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
 Things to do in this subroutine
 1. Parameter/Pointer Adjustments
    - open/past in historical dataset
    - cell->t incremental
 2. Calculate Deterministic Parts in Time Series Models
    - Calculate AR part (historical) for easy simulation of the observations
    - Calculate MA part (historical) for easy simulation of the observations
    * Results are stored into tsType organized using matrix[row-time steps][col-variate] *
 3. Warm Start Dynamics of the Problem
    - Store the most recent incumbent solution and track the changes
    - Track the master problem solution and its changes
    - Update the master problem right hand side with changes of the solution between time periods
    - Track the subprob right hand side value changes
    - Update the subprob right hand side value with the changes
 4. Update gamma Structure
    - Calculate ARMA time series components (both historical and recursive parts)
    - Track changes of the current time series conponents to get ready for warming starting the cell
    - Update pibARMA contents with the latest information
    - For no good reason, Update prob[1]->bBar, change its random row values to 0.0. Do it here since it is a
        nice loop for further validation.
 5. Warm Start Stochastic Components in Problem
    - Update the sigma structure with the changes in subproblem right hand side computed in step 3
    - Update the delta structure with the historical noise
 6. [MAIN] Warm start the algorithmic cell information
    - Transform the cuts by dropping a certain portion of the generated cuts
    - Take the changes of time series components calculated in step 4 and warm up each cut
    - Take the changes of subproblem right hand side calculated in step 3 and warm up each cut
    - Each warm up is conducted by using the defining dual associated with the cut for every obs that the cut
        is based on
 7. Refresh the solnType for indicators to work
    - Clean the piSRatio according to user's configuration
    - Refresh indicators/flags/parameters stored in the solnType
 8. SMILE :-D
 -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
 Note that STEPs [6], [7] shouldn't be conducted when in the initial period.
 Note that STEPs [3] is reduce to minimum functionality (few things completed) when in initial period.
 -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
 *************************************************************************************************************/
int setupStoPeriod(cellType *cell, stocType *stoc, probType **prob, solnType *soln, BOOL *initialPeriod) {
    
#ifdef  ROLL
    printf("~setupStochastic()\n");
#endif
    
    int     i, j;
    int     current = cell->arima->history->open;
    int     status;
    
    /********** Step 1. Parameter/Pointer Adjustment **********/
    if (*initialPeriod == FALSE) {
        newPeriod = TRUE;
        cell->t++;
        cell->arima->history->open += cell->arima->subh;        // More data is revealed
        cell->arima->history->past += cell->arima->subh;
        cell->arima->noise->open += cell->arima->subh;          // Keep every stream synced
        cell->arima->noise->past += cell->arima->subh;
        cell->Kt = cell->k;                                     //Record the last time period iteration count
        cell->Kc = cell->k;
#ifdef  ROLL
        printf(" :: Adjusting stream pointer open = %d; (VAR1) history[open] = %lf; noise[open] = %lf\n",
               cell->arima->history->open, cell->arima->history->obs[cell->arima->history->open][1],
               cell->arima->noise->obs[cell->arima->noise->open][1]);
#endif
        current = cell->arima->history->open;
    } else {
        // Initial Period Explicit Operations
        /*
         When a time series model is incorporated, the mean value for random rows lose its meaning.
         Since this mean value will be utilize when changing subprob right hand side, we want to clear
         the value here before even solving the initial time period cell.
         */
        cell->t = 1;            //Starting with time period counter counter to 1
        cell->Kt = cell->k;     // Initialize the previous iteration counter the cumulative
        cell->Kc = cell->k;     // Initialize the cumulative iteration counter

        for (j=1;j<=cell->arima->subh;j++){
            for (i=1;i<=cell->arima->var;i++){
                prob[1]->bBar->val[cell->arima->variate[i]->row[j]] = 0.0;
            }
        }
    }

#ifdef  ROLL
    printLine();
    if (*initialPeriod == FALSE) {
        printf("Rolling Start from time period %d -> time period %d\n", cell->t, cell->t+1);
    } else {
        printf("Rolling the initial time period time period 0 -> time period 1\n");
        cell->t = 1;
    }
    printf("Calculating new time series deterministic value... \n");
#endif
    
    /********** Collecting Historical Noise for deltaType **********/
        // Note here is the time line
        // .. (current - 1)|        current        |  (current + 1) ..
        // -- past -- .. --|     current(known)    | -- .. -- future -- .. --
        // ..............  | First Historical Data | ........................
#ifdef  ROLL
    printf(" :: Available Historical Data Points :: \n");
    for (j=1; j<=cell->arima->var;j++) {
        for (i=1; i<=cell->arima->p+1; i++) {
            printf("\t :: Idx%d=>%lf; ", current-i+1, cell->arima->history->obs[current-i+1][j]);
        }
        printf("\n");
    }
    printf(" :: Available Historical Noise Points :: \n");
    for (j=1; j<=cell->arima->var; j++) {
        for (i=1; i<=cell->arima->q+1; i++) {
            printf("\t :: Idx%d=>%lf ",current-i+1, cell->arima->noise->obs[current-i+1][j]);
        }
        printf("\n");
    }
    printf("\n");
#endif
    
#ifdef  CALC
    printf("Validating Omega/Delta Integrity.\n");
    printf("Omega->cnt = %d", cell->omega->cnt-1);
    printf("Checking corresponding delta information. \n");
    if (cell->t > 1) {
        if ( !(cell->delta->val[0][cell->omega->cnt-1].C) ) {
            errMsg("Integrity Corrupted", "warmUp", "no delta information detected.", 0);
            return 1;
        } else if ( cell->delta->val[cell->lambda->cnt-1][cell->omega->cnt-1].C == NULL ) {
            errMsg("INtegrity Corrupted", "warmUp", "no delta information detected", 0);
            return 1;
        } else {
            printf("Passed.\n");
        }
    }
#endif
    
    /********** Step 2. Caluclate deterministic parts related to histocial data in time series **********/
    status = warmUpTimeSeries(cell->arima, current);
    if (status) {
        errMsg("warming", "warmUp", "Failed to warm up the stochastic parts", 0);
        return 1;
    }
    
    /********** Step 4. Updating gammaStructure bARMA & pibARMA **********/
    status = warmUpGamma(cell, cell->gamma, cell->arima, cell->t, cell->lambda, stoc);
    if (status) {
        errMsg("warming", "warmUp", "Failed to warm up gamma", 0);
        return 1;
    }

    /********** Step 3. Update Dynamic Parts (aBar/sigma) **********/
    status = warmUpDynamic(prob, soln, cell, cell->master, cell->subprob);
    if (status) {
        errMsg("warming", "setupStochastic", "Failed to warm up the dynamic parts", 0);
        return 1;
    }
#if defined(ROLL)
    //Write down dynamically updated problem for rolling debug
    char dynMasterfname[30] = "/roll/rollDynMaster  .lp";
    char dynSubprobfname[30] = "/roll/rollDynSubprob  .lp";
    dynMasterfname[19] = '0' + cell->t / 10 % 10;
    dynMasterfname[20] = '0' + cell->t / 1 % 10;
    dynSubprobfname[20] = '0' + cell->t / 10 % 10;
    dynSubprobfname[21] = '0' + cell->t / 1 % 10;
    writeProblem(cell->master->lp, dynMasterfname);
    writeProblem(cell->subprob->lp, dynSubprobfname);
#endif

    /********** Step 5. Update the stochastic components **********/
    if (*initialPeriod == FALSE) {
        status = warmStartCell(cell, prob, soln);
        if (status) {
            errMsg("warming", "warmUp", "Failed to warm up cell", 0);
            return 1;
        }

#if defined(ROLL)
        /* Write down warmed problem for rolling debug */
        char wsMasterfname[30] = "/roll/warmStarted  .lp";
        wsMasterfname[17] = '0' + cell->t / 10 % 10;
        wsMasterfname[18] = '0' + cell->t / 1 % 10;
        writeProblem(cell->master->lp, wsMasterfname);
#endif
        
        status = refreshSoln(soln,cell,prob);
        if (status) {
            errMsg("warming", "warmUp", "Failed to clean up the solnType between time periods", 0);
        }
    }
    
#ifdef TUNE
    int     obs, lambdaIdx, sigmaIdx;
    double  argmaxPre = -DBL_MAX, argmaxPost = -DBL_MAX;
    iType   iStar;
    vector  piCbarX;
    
    if ( !(piCbarX= arr_alloc(cell->sigma->cnt, double)) )
        errMsg("Allocation", "refreshSoln", "piCbarX",0);
    
    /* Retrive important lambda */
    for (obs=0; obs<cell->omega->cnt; obs++) {
        iStar = r_computeIstar(prob[1]->num, prob[1]->coord, cell->arima, cell->omega, cell->lambda,
                       cell->gamma, cell->sigma, cell->delta, cell->auxDelta, cell->trimDelta,
                       soln->incumbX, piCbarX, &argmaxPre, &argmaxPost, obs, cell->k, 0);
        cell->lambda->alive[iStar.delta] = 99;
    }
    
    for (lambdaIdx=0; lambdaIdx < cell->lambda->cnt; lambdaIdx++ )
        if (cell->lambda->alive[lambdaIdx] == 99)
            cell->lambda->alive[lambdaIdx] = SURVIVED;
        else
            cell->lambda->alive[lambdaIdx] = DROPPED;
    
    
    /* Re-set simga linkers */
    cell->sigma->aCnt = 0;
    for (sigmaIdx=0; sigmaIdx<cell->sigma->cnt; sigmaIdx++) {
        cell->sigma->linker[sigmaIdx] = DROPPED; //Refresh the linker for safety
        if (cell->lambda->alive[cell->sigma->lambdaIdx[sigmaIdx]] == SURVIVED) {
            cell->sigma->linker[cell->sigma->aCnt] = sigmaIdx;
            cell->sigma->aCnt++;
        }
    }
    mem_free(piCbarX);
#endif
    
    /********** Step 10. :-D **********/
    /* Ready for next time period */
    (*initialPeriod) = FALSE; //Flip the initial period flag
    
    if (run.WARM_START_LB == 1) {
        printf("Validating warm started cust lower bound ... \n");
        /* Conduct bootStrapCheck on warmed up cuts */
        // Currently, the bootStrapCheck check all cuts to detect violation
        if (cell->t > 1) {
            BOOL bsCheck = TRUE;
            status = lowerBoundCheck(prob, cell, soln, &bsCheck);
            if (status) {
                errMsg("validation", "solveCell", "Failed to run bootStarpCheck", 0);
                return 1;
            }
            if (bsCheck==FALSE) {
                printf("\t :: FAILED THE LOWER BOUND TEST !! :-( \n");
            } else {
                printf("\t :: PASSED THE LOWER BOUND TEST !! :-D \n");
            }
        }
    }

    printf("Rolling Completed\n");
    printLine();

#ifdef ROLL_VALID
    printf("Validating warm started stochastic structure ... ");
    
    double  objV, objEst;
    int     stat1;
    vector  rhs, reObs, c_reObs;
    intvec  indices;
    
    sigmaIdx = 0;
    
    reObs = (vector) arr_alloc(cell->arima->subh * cell->arima->var + 1, double);
    c_reObs = (vector) arr_alloc(cell->arima->subh * cell->arima->var + 1, double);
    
    indices = (intvec) arr_alloc(prob[1]->num->rows, int);
    for (i=0; i<prob[1]->num->rows; i++) indices[i] = i;
    
    for (obs=0; obs<cell->omega->cnt; obs++) {
        // Re-constructing observation
        for (i=1; i<=cell->arima->subh * cell->arima->var; i++)
            reObs[i] = cell->omega->vals[obs][i] + cell->auxOmega->vals[obs][i] + cell->gamma->endo->bARMA->val[i];
        
        // Re-simulating scenarios for checking
        generateARIMA(cell->t, cell, c_reObs, NULL, NULL, NULL, cell->omega->vals[obs], NULL);
        for (i=0; i<cell->arima->subh * cell->arima->var; i++) {
            if ( fabs(c_reObs[i] - reObs[i+1]) > run.EPSILON )
                printf("\t :: Re-SIM OBS %d -> mismatch position %d [%lf, %lf]\n", obs, i, c_reObs[i], reObs[i+1]);
        }
        
        
        // Re-consructing subprob rhs
        rhs = computeRHS(prob[1]->num, prob[1]->coord, prob[1]->bBar, prob[1]->Cbar, soln->incumbX, c_reObs-1);
        if (rhs == NULL) {
            errMsg("CALC","setupStoPeriod","failed to calculate right hand side value",0);
            return 1;
        }
        
        // Change subproblem right hand side
        status = changeRHS(cell->subprob->lp, prob[1]->num->rows, indices, rhs+1);
        if (status) {
            errMsg("CPLEX", "setupStoPeriod", "Failed to change subprob right hand side", 0);
            return 1;
        }
        
        // Solve subproblem
        status = solveProblem(cell->subprob->lp, cell->subprob->name, cell->subprob->type, &stat1);
        if ( status ) {
            printf("Subproblem is infeasible: check right hand side.\n");
            printVector(c_reObs-1, prob[1]->num->numRV, NULL);
            return 1;
        }
        
        // Get primal objective value
        objV = getObjective(cell->subprob->lp, PROB_LP);
        
        // Locate sigmaIdx
        sigmaIdx = cell->omega->sigmaIdx[obs];
        
        // Locate lambdaIdx
        lambdaIdx = cell->sigma->lambdaIdx[sigmaIdx];
        
        // Get dual estimate with updated structures, should obtain lower bound (objEst <= objV)
        objEst = cell->sigma->vals[sigmaIdx].b
            + cell->gamma->endo->pibARMA[lambdaIdx]
            - vXv(cell->sigma->vals[sigmaIdx].C, soln->incumbX, prob[1]->coord->colsC, prob[1]->num->cntCcols)
            + cell->delta->val[lambdaIdx][obs].b
            + cell->auxDelta->val[lambdaIdx][obs].b
            + cell->trimDelta->val[lambdaIdx][obs].b
            - vXv(cell->delta->val[lambdaIdx][obs].C, soln->incumbX, prob[1]->coord->rvCols, prob[1]->num->rvColCnt)
            - vXv(cell->auxDelta->val[lambdaIdx][obs].C, soln->incumbX, prob[1]->coord->rvCols, prob[1]->num->rvColCnt);
        
        if ( objV - objEst > - run.EPSILON) {
            printf("\t :: OBS %d (CK %d, SIGMA %d) -> LB(%lf) vs. UB(%lf) \n",
                   obs+1, cell->sigma->ck[sigmaIdx], sigmaIdx, objEst, objV);
        }else {
            errMsg("Validation", "steupStoPeriod", "Validation of warm start fails", 0);
        }
        
        mem_free(rhs);
    }
    
    mem_free(indices);
    mem_free(reObs);
    mem_free(c_reObs);
    
    printf(" DONE\n");
#endif
    
    return 0;
}//END of steupStochastic()

/*
 This subroutine is desinged to re-fresh a cell to get it ready for the next time period. Different from warm-up a cell, the cell
 is cleaned by taking out all cuts, all stochastic information, and all gamma information. The time series model and dyanmic components
 will be updated as usualy after the refresh is conducted.
 */
int refreshPeriod(cellType *cell, stocType *stoc, probType **prob, solnType *soln, BOOL *initialPeriod) {
    
#ifdef ROLL
    printf("~setupCleanStochastic");
#endif
    
    int i, j, length;
    int current = cell->arima->history->open;
    int status;
    
    /* Update the pointers/counters in the cell structure */
    if (*initialPeriod == FALSE) {
        newPeriod = TRUE;
        cell->Kt = 0;                                       // Keep track of the cumulative iterations until now
        cell->Kc += cell->k;
        cell->k = 0;                                        // Refresh the iteration counter
        cell->t++;
        cell->arima->history->open += cell->arima->subh;    // More data is revealed
        cell->arima->history->past += cell->arima->subh;
        cell->arima->noise->open += cell->arima->subh;      // Keep every stream synced
        cell->arima->noise->past += cell->arima->subh;
#ifdef  ROLL
        printf(" :: Adjusting stream pointer open = %d; (VAR1) history[open] = %lf; noise[open] = %lf\n",
               cell->arima->history->open, cell->arima->history->obs[cell->arima->history->open][1],
               cell->arima->noise->obs[cell->arima->noise->open][1]);
#endif
        current = cell->arima->history->open;
    } else {
        /* ---------------------------------------------------
                  Initial Period Explicit Operations
         When a time series model is incorporated, the mean value
         for random rows lose its meaning. Since this mean value 
         will be utilize when calculating omega for RHS, we
         want to clear the value here before even solving the 
         initial time period cell.
         ----------------------------------------------------- */
        cell->t = 1;
        cell->Kt = cell->k;
        cell->Kc = cell->k;
        for (i=1;i<=cell->arima->var;i++){
            for (j=1;j<=cell->arima->subh;j++){
                prob[1]->bBar->val[cell->arima->variate[i]->row[j]] = 0.0;
            }
        }
    }
    
    /* Section II :: Freeing and Reconstructing */
    if (*initialPeriod == FALSE) {
        /* Clean out several algorithm */
        freeDeltaType(cell->delta, cell->lambda->cnt, cell->omega->cnt);
        freeDeltaType(cell->auxDelta, cell->lambda->cnt, cell->omega->cnt);
        freeDeltaType(cell->trimDelta, cell->lambda->cnt, cell->omega->cnt);
        freeLambdaType(cell->lambda);
        freeSigmaType(cell->sigma);
        freeOmegaType(cell->omega);
        freeOmegaType(cell->auxOmega);
        freeOmegaType(cell->trimOmega);
        freeCuts(cell->cuts);
        freeOneProblem(cell->master);
        
        length = (cell->maxIter + cell->maxIter / cell->tau + 1) * run.HORIZON;
        
        /* Re-build clean components of an algorithm cell to get ready for the next time period */
        cell->cuts = newCuts(cell->maxCuts);
        cell->lambda = newLambda(length, 0, prob[1]->num->rvRowCnt, prob[1]->coord);
        cell->sigma = newSigma(length, prob[1]->num->rvColCnt, 0, prob[1]->coord);
        cell->delta = newDelta(length, prob[1]->coord);
        cell->auxDelta = newDelta(length, prob[1]->coord);
        cell->trimDelta = newDelta(length, prob[1]->coord);
        cell->omega = newOmega(cell->maxIter * run.HORIZON + run.HORIZON);
        cell->auxOmega = newOmega(cell->maxIter * run.HORIZON + run.HORIZON);
        cell->trimOmega = newOmega(cell->maxIter * run.HORIZON + run.HORIZON);
        
        /* Update some parameters in cellType */
        cell->quadScalar = run.MIN_QUAD_SCALAR;
        
        /* Refresh gamma type :: cleaning previously known Pi x Deterministic ARMA */
        cell->gamma->endo->pibCnt = 0;
        mem_free(cell->gamma->endo->pibARMA);
        if (!(cell->gamma->endo->pibARMA = (vector) arr_alloc(length, double)))
            errMsg("allocation", "setupCleanPeriod", "refresh gamma->endo->pibMA", 0);
        
        /* Reconstruct the previous time period right hand side :: no updated right hand side from dynamic/stochasticity */
        /* Updates on the one problem structure and probType will be deliveried in warmUpDynamic */
        cell->master = newStoMaster(prob[0]->sp, prob[0]->Bbar[0], cell->cuts, cell->maxCuts);
        if ( cell->master == NULL ) {
            errMsg("setup", "newCell", "failed to setup the master problem", 0);
            return 1;
        }

    }
    
    status = warmUpTimeSeries(cell->arima, current);
    if (status) {
        errMsg("warming", "setupCleanPeriod", "Failed to warm up the stochastic parts", 0);
        return 1;
    }
    
    status = warmUpGamma(cell, cell->gamma, cell->arima, cell->t, cell->lambda, stoc);
    if (status) {
        errMsg("gamma warm-up", "setupCleanPeriod", "Failed to warmup gammaType for next time period", 0);
        return 1;
    }
    
    status = warmUpDynamic(prob, soln, cell, cell->master, cell->subprob);
    if (status) {
        errMsg("dynamic warm-up", "setupCleanPeriod", "failed to warmUpDynamic for next time period", 0);
        return 1;
    }
    
#if defined(ROLL)
    //Write down dynamically updated problem for rolling debug
    char dynMasterfname[30] = "/roll/rollDynMaster  .lp";
    char dynSubprobfname[30] = "/roll/rollDynSubprob  .lp";
    dynMasterfname[19] = '0' + cell->t / 10 % 10;
    dynMasterfname[20] = '0' + cell->t / 1 % 10;
    dynSubprobfname[20] = '0' + cell->t / 10 % 10;
    dynSubprobfname[21] = '0' + cell->t / 1 % 10;
    writeProblem(cell->master->lp, dynMasterfname);
    writeProblem(cell->subprob->lp, dynSubprobfname);
#endif
    
    if (*initialPeriod == FALSE) {
        status = refreshSoln(soln,cell,prob);
        if (status) {
            errMsg("dynamic warm-up", "setupCleanPeriod", "Failed to clean up the solnType", 0);
        }
    }
    
    /* Ready for next time period */
    (*initialPeriod) = FALSE; //Flip the initial period flag
    
    return 0;
}

/*******************************************************************************************************************
 This subroutine is designed for updating the deterministic parts in time series models. When calculating
 time series model, it is possible that some terms can be pre-calculated and stored into the rolling horizon
 structure. These terms are viewed as the gamma type of the implementation. To obtain the gamma type, we need
 to conduct calculation within the time sereis model. These calculation use the values that has been pro-processed
 and hence cannot be directly used for optimization model. Since these values are directly associated with time
 series model. We conduct the calculation within tsType for time series value only and store them in the way that
 time series model would benefit the most.
 -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
 Deterministic Data in tsType is store for easy time series calculation.
 In time series model, matrix calculation Phi * X_t is common
 Hence, detARObs & detMAObs (both are parallel vectors for different time points) are stored in this form
 VAR1    VAR2    ....     VARV
 00:10   [1][1]  [1][2]   ....   [1][var]    --Vector[time][var]
 00:20   [2][1]  [2][2]   ....   [2][var]    --Vector[time][var]
 ...     ....    ....     ....     ....      --....
 01:00   [N][1]  [N][2]   ....   [N][var]    --Vector[time][var]
 You can fetch each time interval really easily. Once you say detARObs[1], you are fetching all variates
 deterministic value at hour/subhour 1.
 Note that everything realted to time series will start its indexing at 1. 0 is empty and no memory is allocated.
 *******************************************************************************************************************/
int warmUpTimeSeries(tsType *arima, int current){
    

#ifdef  ROLL
    printf("\t~warmUpTimeSeries()\n");
#endif
    
    int i, j;
    if (arima->histARObs) {
        for (i=1; i<=arima->subh; i++)
            if (arima->histARObs[i]) mem_free(arima->histARObs[i]);
        mem_free(arima->histARObs);
    }
    if (arima->histMAObs) {
        for (i=1; i<=arima->subh; i++)
            if (arima->histMAObs[i]) mem_free(arima->histMAObs[i]);
        mem_free(arima->histMAObs);
    }
    
    // Memory allocations for pre-calculated information
    if ( !(arima->histARObs = (vector *) arr_alloc(arima->subh+1, vector)) )
        errMsg("allocation","readArma","Failed to allocate memory to detARObs",0);
    if ( !(arima->histMAObs = (vector *) arr_alloc(arima->subh+1, vector)) )
        errMsg("allocation","readArma","Failed to allocate memory to detARObs",0);
    for (i=1; i<=arima->subh; i++) {
        if ( !(arima->histARObs[i] = (vector) arr_alloc(arima->var+1, double)) )
            errMsg("allocation","readArma","detARObs->vectors",0);
        if ( !(arima->histMAObs[i] = (vector) arr_alloc(arima->var+1, double)) )
            errMsg("allocation","readArma","detMAObs->vectors",0);
    }
    
    /* ~~~> ARIMA -> detAR Section: The window slides backward from current time to length p. */
    for (i=1; i<=arima->subh; i++) {
        // Initialize Values in Deterministic Values
        for (j=1; j<=arima->var; j++) arima->histARObs[i][j] = 0.0;
        // Start Calculation
        if ( (arima->p)-i+1 <= 0) {
            // If p sliding window slides out of historical data frame -> get zero vector
            for (j=1; j<=arima->var; j++) arima->histARObs[i][j] = 0.0;
        } else {
            // If p sliding window slides still within some parts of the historical data frame, calculate
            for (j=1; j<=max(0,arima->p-i+1); j++) {
#ifdef ROLL
                // Sanity Check
                if (j+i-1 == 0 || current-j+1 > current) {
                    errMsg("calculation", "warmUpTimeSeries", "wrong indexing when calculating arima->detARObs", 0);
                    return 1;
                }
#endif
                arima->histARObs[i] = MSparsexvAdd(arima->Phi[j+i-1], arima->history->obs[current-j+1], arima->histARObs[i]);
            }
        }
    }
    
    /* ~~~> ARIMA -> detMA Section: The window slides backward from current time to length q. */
    for (i=1; i<=arima->subh; i++) {
        // Initialize Value
        for (j=1; j<=arima->var; j++) arima->histMAObs[i][j] = 0.0;
        // Start Calculation
        if ( (arima->q)-i+1 <= 0) {
            for (j=1; j<=arima->var; j++) arima->histMAObs[i][j] = 0.0;
        } else {
            for (j=1; j<=max(0,arima->q-i+1); j++) {
#ifdef ROLL
                // Sanity Check
                if (j+i-1 == 0 || current-j+1 > current) {
                    errMsg("calculation", "warmUpTimeSeries", "wrong indexing when calculating arima->detMAObs", 0);
                    return 1;
                }
#endif
                arima->histMAObs[i] = MSparsexvAdd(arima->Theta[j+i-1], arima->noise->obs[current-j+1], arima->histMAObs[i]);
            }
        }
    }
    
    printf("Warming Start time series strcutrue ... Done\n");
    return 0;
}//END of warmUpTimeSeries()

/*-------------------------------------------------------------------------------------------------
 Time Series Results are Fetched from tsType. The value fetch back is without reverse-
 pre-process. These values are calculated in this subroutine but stored in tsType for a
 handy-for-checking purpose.
 NOTE HERE, the reverse process is taken care of in here. What we stored into gamma is a
 product depends on what type of pre-process we utilized for time series model. This value
 is used for cut generation and warmUp period.
 
 If no pre-process is applied, this is what we stored:
 tempAR = tAR = sum( j -> p ){ Phi_j * Omega_(t-j) }
 tempMA = tMA = sum( j -> j, q ){ Theta_j * Epsilon_(t-j) }
 Note that if j>=p / j>=q, we don't do any loop.
 Also note that the reverse is correspondingly for each different variate. Hence, we could end
 up looping j and cell->t for different variate.
 
 If EMP pre-process is applied, this is what we stored:
 tempAR = tAR * sigma + mu
 tempMA = tMA * sigma
 Note that if we have hourly model, sigma and mu follows cell->t to pick value.
 If we have sub-hourly model, sigma and mu follows j and cell->t as guidelines for sub-hourly
 intervals and hourly period when picking value from arima->empSD and arima->empMean.
 
 If STL/CLA pre-process is applied, this is what we stored:
 tempAR = tAR + season + trend
 tempMA = tMA
 This is similar to EMP, but easier.
 
 When updating gammaType(warming up), we don't need to worry about updating these terms. Since the
 calculation of (time series contribution x lambda_pi) will be conducted again for all lambda_pi.
 -------------------------------------------------------------------------------------------------*/
int warmUpGamma(cellType *cell, gammaType *gamma, tsType *arima, int t, lambdaType *lambda, stocType *stoc) {
    
#ifdef  ROLL
    printf("\t~warmUpGamma()\n");
#endif
    
    int i, j, k = 1;
    int arj, maj, current;
    double  tempARMA;
    vector *gammaObs;

    
    current = arima->history->open;
    
    // Memory allocations for pre-calculated information
    if ( !(gammaObs = (vector *) arr_alloc(arima->subh+1, vector)) )
        errMsg("allocation","readArma","Failed to allocate memory to detARObs",0);
    for (i=1; i<=arima->subh; i++) {
        if ( !(gammaObs[i] = (vector) arr_alloc(arima->var+1, double)) )
            errMsg("allocation","readArma","detARObs->vectors",0);
    }
    
#ifdef ROLL
    printf("\t :: Validation -> ");
    vector zeroNoise;
    vector ansARMA;
    zeroNoise = (vector) arr_alloc(stoc->numOmega+1, double);
    ansARMA = (vector) arr_alloc(stoc->numOmega+1, double);
    generateARIMA(cell->t, cell, ansARMA, NULL, NULL, NULL, zeroNoise, &run.RUN_SEED);
    printVector(ansARMA-1, stoc->numOmega,NULL);
    mem_free(zeroNoise); mem_free(ansARMA);
#endif
    
    /* Dummy simulation of the inherited deterministic parts and recursive deterministic parts */
    for (i=1; i<=arima->subh; i++) {
        // Initialize Values in Deterministic Values
        for (j=1; j<=arima->var; j++) gammaObs[i][j] = 0.0;
        // Start Calculation from AR section
        for (arj=1; arj<=arima->p; arj++) {
            if (arj < i) { // When AR calculation slides out of historical window, use previously generation observations for new generations
                gammaObs[i] = MSparsexvAdd(arima->Phi[arj], gammaObs[i-arj], gammaObs[i]);
            } else { // When AR calculation slides into the historical window, use historical data for generation
                gammaObs[i] = MSparsexvAdd(arima->Phi[arj], arima->history->obs[current-arj+i], gammaObs[i]);
            }
            
        }
        
        // Add MA Calculation Results Calculation
        for (maj=1; maj<=arima->q; maj++) {
            if ( maj < i) {
                gammaObs[i][maj] += 0.0;
            } else {
                gammaObs[i] = MSparsexvAdd(arima->Theta[maj], arima->noise->obs[current-maj+i], gammaObs[i]);
            }
        }
    }
    
    // Loop through each variate and its each sub-hourly time steps :: ARRANGEMENT IN tsType
    for (j=1; j<=arima->subh; j++) {
        for (i=1; i<=arima->var; i++) {
            tempARMA = reverse(gammaObs[j][i],i,j,t,arima) * run.SCALAR;
            
#ifdef ROLL
            printf("\t :: Calulated temp with var=%d | subh=%d | ARMA(%lf) \n",
                   i,j,tempARMA);
#endif
            /* Storing values into gamma types for ARMA :: Also comput difference value DELTA_ARMA */
            gamma->endo->bARMA->cnt = k;
            gamma->endo->changebARMA->cnt = k;
            gamma->endo->changebARMA->val[k] = tempARMA - gamma->endo->bARMA->val[k];
            gamma->endo->changebARMA->col[k] = arima->variate[i]->row[j];
            gamma->endo->bARMA->val[k] = tempARMA;
            gamma->endo->bARMA->col[k] = arima->variate[i]->row[j];
#ifdef ROLL
            printf("\t :: Change ARMA value with var=%d | subh=%d | changebARMA->val[%d] = %lf \n",
                   i, j, k, gamma->endo->changebARMA->val[k]);
#endif
            k++; //Next random row...
        }
    }
    
    /* Updateing pibAR & pibMA */
    for (i=0; i<lambda->cnt; i++)
        gamma->endo->pibARMA[i] = vXv(lambda->vals[i], gamma->endo->bARMA->val, NULL, stoc->numOmega);
    
    for (i=1; i<=arima->subh; i++)
        if (gammaObs[i]) mem_free(gammaObs[i]);
    mem_free(gammaObs);
    
    printf("Warming start gamma structure ... DONE\n");
    return 0;
}//END of warmUpGamma()


/***********************************************************************************************************************************
 * This subroutine is designed for updating the dynamics between time periods. In the current problem design, the dynamics mainly
 * contains two parts: master stage primal solution dynamics and subprob right hand side. To do the updating work, the difference
 * of corresponding values is calculated and used to update the problem right hand side. This is because the problem right hand
 * side is not solely contributed by the dynamic terms. Besides the updating, we also store the changes of each row so as to update
 * the stochastic components.
 ***********************************************************************************************************************************/
int warmUpDynamic(probType **prob, solnType *soln, cellType *cell, oneProblem *master, oneProblem *subprob){
    
#ifdef ROLL
    printf("\t~warmUpDynamic()\n");
#endif
    
    int     k = 0, i, j, status = 0, masterRowCnt, masterIdx, subprobIdx, stat1, stat2;
    double  memory, change;
    vector  masterRHS = NULL, subprobRHS = NULL, masterRHSChange = NULL, subprobRHSChange = NULL, history = NULL, fullObs = NULL, tempY = NULL;
    intvec  indicesM = NULL, indicesS = NULL;
    BOOL    detFLAG; //It is default to consider this subroutines works for stochastic problem
    
    /* Things to do in this subroutine:
        
        0. Collect x and associated difference
        1. B00 * x
        2. Calculate y and associated difference
        3. B10 * y
        4. B11 * y
        5. Record BDelta
        6. Replace a0
        7. Replace a1
        8. Record aDelta
     
     */
    
    /* Initial Check on the problem type we are working on. If deterministic problem we adjust the offset appropriately. */
    if (subprob == NULL) {
        detFLAG = TRUE;
        masterIdx = 0;
        subprobIdx = 0;
        subprob = master;
    } else {
        detFLAG = FALSE;
        masterIdx = 0;
        subprobIdx = 1;
    }
    
   /*------- Section A: Endogenous Dynamic Connector prob[0]->Cbar for master problem -------*/
    if (cell->t > 1) {
        // Current Solution is in soln->incumbX || Previous Solution is in cell->gamma->endo->optX
        
        masterRowCnt = getNumRows(master->lp); //Size of master problem
        
        /* Memory Allocation and Fetch Right Hand Side from Master Problem */
        if ( !(masterRHS = (vector) arr_alloc(masterRowCnt+1 , double)) )
            errMsg("allocation","warmUpDynamic","masterRHS",0);
        if ( !(masterRHSChange = (vector) arr_alloc(masterRowCnt+1, double)))
            errMsg("allocation", "warmUpDynamic", "masterRHSChange", 0);
        if ( getRhsx(master->lp, 0, masterRowCnt, masterRHS) ) {
            errMsg("solver","warmUpDynamic","Failed to reach current master problem right hand side",0);
            return 1;
        } // Notice that masterRHS starts from 0
        
        /* Construct master indices vector */
        if ( !(indicesM = arr_alloc(masterRowCnt,int)) )
            errMsg("allocation","warmUpDynamic","failed to allocate memory to indices",0);
        for (i=0; i<masterRowCnt; i++) indicesM[i] = i;
        
        /* Track changes in master solution and record new master solution*/
        for (i=1; i<=prob[0]->sp->mac; i++) {
            cell->gamma->endo->optXChange[i] = soln->incumbX[i] - cell->gamma->endo->optX[i];
            cell->gamma->endo->optX[i] = soln->incumbX[i];
        }

        /* Calculate B_{00}: master right hand side changes contributed from previous master problem */
        if (prob[0]->Bbar[0]) {
#if defined(ROLL)
            printf("\t :: Updating B_{00} x X -> %d Rows to update...\n", prob[0]->Bbar[0]->cnt);
#endif
            MSparsexvAdd(prob[0]->Bbar[0], cell->gamma->endo->optXChange, masterRHSChange);
            for (i=1; i<=prob[0]->num->rows; i++) {
                masterRHS[i-1] -= masterRHSChange[i];
#if defined(ROLL) && defined(ROLL_DETAIL)
                if (fabs(masterRHSChange[i]) > run.EPSILON)
                    printf("\t :: master row %d NEW[%f] = OLD[%f] - CHANGE[%f]\n",
                           i, masterRHS[i-1], masterRHS[i-1]+masterRHSChange[i], masterRHSChange[i]);
#endif
            }
            //One Step Solution
            // MSparsexvAdd(prob[0]->Bbar[0], cell->gamma->endo->optX, masterRHS-1);
            
            // Also update the prob[0]->sp->rhsx which is used in many places
            prob[0]->sp->rhsx = MSparsexvSub(prob[0]->Bbar[0], cell->gamma->endo->optXChange, prob[0]->sp->rhsx-1) + 1;
            cell->master->rhsx = MSparsexvSub(prob[0]->Bbar[0], cell->gamma->endo->optXChange, cell->master->rhsx-1) + 1;
        }

        /* Calculate a_{0} */ // Notice that cell->t has already been updated
        if (prob[0]->aBar) {
#if defined(ROLL)
            printf("\t :: Updating aBar x X -> %d Rows to update...\n", prob[0]->aBar->cnt);
#endif
            for (i=1; i<=prob[0]->aBar->cnt; i++) {
                memory = masterRHS[prob[0]->aBar->col[i]-1]; //Record the old coefficient
                //Calculate change value
                change = prob[0]->aBar->val[i] * (cell->gamma->exo->aData[0][cell->t][i] - cell->gamma->exo->aData[0][cell->t-1][i]);
                //Perform change
                masterRHS[prob[0]->aBar->col[i]-1] += change;
                //Record correspondingly in the algorithm structure
                prob[0]->sp->rhsx[prob[0]->aBar->col[i]-1] += change;
                prob[0]->bBar->val[prob[0]->aBar->col[i]] += change;
                cell->master->rhsx[prob[0]->aBar->col[i]-1] += change;
#if defined(ROLL) && defined(ROLL_DETAIL)
                printf("\t :: master row %d NEW[%f] = OLD[%f] - CHANGE[%f]\n",
                       prob[0]->aBar->col[i], masterRHS[prob[0]->aBar->col[i]-1], memory,
                       prob[0]->aBar->val[i] * (cell->gamma->exo->aData[0][cell->t][i] - cell->gamma->exo->aData[0][cell->t-1][i]));
#endif
            }
        }
    
        
        /* Subproblem only exist when solving problme stochastically */
        if (run.DETERMINISTIC == 0) { //Subproblem treatment when solving stochastic problem
            
            //No need to allocate memory to subprobRHS since we compuate it later
            if ( !(subprobRHSChange = (vector) arr_alloc(subprob->mar + 1, double)))
                errMsg("allocation", "warmUpDynamic", "masterRHSChange", 0);
            
            if ( !(indicesS = arr_alloc(subprob->mar + 1,int)) )
                errMsg("allocation","warmUpDynamic","failed to allocate memory to indices",0);
            for (i=0; i<subprob->mar+1; i++) indicesS[i] = i;
            
            /* Memory allocation of a temporary subprob primal solution vector in order to track changes */
            if (!(tempY = (vector) arr_alloc(subprob->mac+1, double)))
                errMsg("allocation", "warmUpDynamic", "tempY", 0);
            
            
            /* Calculate and collect subproblem solution given x_bar and historical information */
            if ( !(history = (vector) arr_alloc(cell->arima->var * cell->arima->subh + 1, double)) ) {
                errMsg("allocation","warmUPDynamic","historyNoise",0);
                return 1;
            }
            
            if ( !(fullObs = (vector) arr_alloc(cell->arima->subh * cell->arima->var + 1, double)) ) {
                errMsg("allocation", "warmUpDynamic", "fullObs", 0);
                return 1;
            }
            
            /* Fetch the previous time period history = what just happened to re-create the historical subprob */
            for (i=1; i<=cell->arima->subh; i++)
                for (j=1; j<=cell->arima->var; j++) {
                    history[k] = cell->arima->history->obs[cell->arima->history->open - cell->arima->subh + i][j];
                    fullObs[k] = max(0, reverse(history[k], j, i, cell->t - 1, cell->arima)) * run.SCALAR;
#if defined(ROLL) && defined(ROLL_DETAIL)
                    printf("\t :: IDX %d :: Fetching historical realization var(%d), subh(%d) :: %lf\n",
                           cell->arima->history->open - cell->arima->subh + i, j, i, fullObs[k]);
#endif
                    k++;
                }
            
            /* Send the right hand side back to a subproblem and solve it */
            subprobRHS = computeRHS(prob[1]->num, prob[1]->coord, prob[1]->bBar, prob[1]->Cbar, soln->incumbX, fullObs-1);
            if ( subprobRHS == NULL ) {
                errMsg("algorithm", "warmUpDynamic", "failed to compute subproblem right-hand side", 0);
                return 1;
            }
            // Different from the masterRHS, which start from 0, this subprobRHS started from 1
            
            /* change the right-hand side in the solver */
            stat1 = changeRHS(cell->subprob->lp, prob[1]->num->rows, indicesS, subprobRHS+1); //Verify if this +1 is leagl or not
            if ( stat1 ) {
                errMsg("solver", "warmUpDynamic", "failed to change the right-hand side in the solver",0);
                return 1;
            }
            
            /* Solve subproblem */
            stat1 = solveProblem(subprob->lp, subprob->name, subprob->type, &stat2);
            if ( stat1 ) {
                if ( stat2 == STAT_INFEASIBLE ) {
                    printf("Subproblem is infeasible with historical information which should not be possible.\n");
                    return 1;
                } else {
                    errMsg("algorithm", "warmUpDynamic", "failed to solve subproblem in solver", 0);
                    return 1;
                }
            }
            
            /* Fetch subprob primal solution */
            status = getPrimal(cell->subprob->lp, tempY, subprob->mac);
            if (status) {
                errMsg("solver", "warmUpDynamic", "failed to fetch the historical subprob primal solution", 0);
            }
            
            /* Track historical subprob primal solution change and record new historical subprob */
            for (i=1; i<=subprob->mac; i++) {
                cell->gamma->endo->optYChange[i] = tempY[i] - cell->gamma->endo->optY[i];
                cell->gamma->endo->optY[i] = tempY[i];
            }
            
            /* Calculate B_{10}: master right hand side changes contributed from previous master problem */

            if (prob[1]->Bbar[0]) {
#if defined(ROLL)
                printf("\t :: Updating B_{10} x X -> %d Rows to update...\n", prob[1]->Bbar[0]->cnt);
#endif
                for (i=0; i<=masterRowCnt; i++)
                    masterRHSChange[i] = 0.0;
                MSparsexvAdd(prob[1]->Bbar[0], cell->gamma->endo->optYChange, masterRHSChange);
                for (i=1; i<=prob[0]->num->rows; i++) {
                    masterRHS[i-1] -= masterRHSChange[i];
#if defined(ROLL) && defined(ROLL_DETAIL)
                    if (fabs(masterRHSChange[i]) > run.EPSILON)
                        printf("\t :: master row %d NEW[%f] = OLD[%f] - CHANGE[%f]\n",
                               i, masterRHS[i-1], masterRHS[i-1] + masterRHSChange[i], masterRHSChange[i]);
#endif
                }
                prob[0]->sp->rhsx = MSparsexvSub(prob[1]->Bbar[0], cell->gamma->endo->optYChange, prob[0]->sp->rhsx-1) + 1;
                cell->master->rhsx = MSparsexvSub(prob[1]->Bbar[0], cell->gamma->endo->optYChange, cell->master->rhsx-1) + 1;
            }
            
            /* Calculate B_{11}: master right hand side changes contributed from previous master problem */
            if (prob[1]->Bbar[1]) {
#if defined(ROLL)
                printf("\t :: Updating B_{11} x X -> %d Rows to update...\n", prob[1]->Bbar[1]->cnt);
#endif
                MSparsexvAdd(prob[1]->Bbar[1], cell->gamma->endo->optYChange, subprobRHSChange);
                for (i=1; i<=prob[1]->Bbar[1]->cnt; i++) {
                    subprobRHS[prob[1]->Bbar[1]->row[i]] -= subprobRHSChange[prob[1]->Bbar[1]->row[i]];
                    cell->gamma->exo->BDelta[i] = -subprobRHSChange[prob[1]->Bbar[1]->row[i]];
#if defined(ROLL) && defined(ROLL_DETAIL)
                    if (fabs(subprobRHSChange[prob[1]->Bbar[1]->row[i]]) > run.EPSILON)
                        printf("\t :: subprob row %d NEW[%f] = OLD[%f] - CHANGE[%f]\n",
                               i, subprobRHS[prob[1]->Bbar[1]->row[i]],
                               subprobRHS[prob[1]->Bbar[1]->row[i]] + subprobRHSChange[prob[1]->Bbar[1]->row[i]],
                               subprobRHSChange[prob[1]->Bbar[1]->row[i]]);
#endif
                }
            }
            
#if defined(ROLL)
            printf("\t :: Updating aBar x X -> %d Rows to update...\n", prob[1]->aBar->cnt);
#endif
            if (prob[1]->aBar) {
                /* Calculate a_{1} */ // Notice that cell->t has already been updated
                for (i=1; i<=prob[1]->aBar->cnt; i++) {
                    
                    memory = subprobRHS[prob[1]->aBar->col[i]]; //Record the old coefficient
                    
                    change = prob[1]->aBar->val[i] * (cell->gamma->exo->aData[1][cell->t][i] - cell->gamma->exo->aData[1][cell->t - 1][i]);
                    
                    subprobRHS[prob[1]->aBar->col[i]] += change;
                    prob[1]->sp->rhsx[prob[1]->aBar->col[i]-1] += change;
                    prob[1]->bBar->val[prob[1]->aBar->col[i]] += change;
                    cell->gamma->exo->aDelta[i] = change;
                    
#if defined(ROLL) && defined(ROLL_DETAIL)
                    printf("\t :: subprob row %d NEW[%f] = OLD[%f] - CHANGE[%f]\n",
                           prob[1]->aBar->col[i], subprobRHS[prob[1]->aBar->col[i]], memory, change);
#endif
                }
                
                //One Step Solution
                //MSparsexvAdd(prob[0]->Bbar[0], cell->gamma->endo->optX, masterRHS);
                
                /* Send ready rhs back to subprob */
                //Notice that subprobRHS is obtained through subroutine computeRHS, index starts from 1
                status = changeRHS(subprob->lp, subprob->mar, indicesS, subprobRHS + 1);
                if (status) {
                    errMsg("solver", "warmUpDynamic", "Failed to change subprob right hand side", 0);
                    return 1;
                }
            }
        }
        
        // Send Ready rhs back to master problem
        status = changeRHS(master->lp, masterRowCnt, indicesM, masterRHS);
        if (status) {
            errMsg("solver","warmUpDynamic","Failed to update the warm-up right hand side",0);
            return 1;
        }
    }
    
    if (subprobRHS)         mem_free(subprobRHS);
    if (masterRHS)          mem_free(masterRHS);
    if (indicesM)           mem_free(indicesM);
    if (indicesS)           mem_free(indicesS);
    if (masterRHSChange)    mem_free(masterRHSChange);
    if (subprobRHSChange)   mem_free(subprobRHSChange);
    if (history)            mem_free(history);
    if (fullObs)            mem_free(fullObs);
    if (tempY)              mem_free(tempY);
    
    printf("Warming start problem dynamically ... DONE\n");
    return 0;
}//End of warmUpDynamic


/**********************************************************************************************************************
 * This subroutine is designed to be the final step of warming up between time periods. All changes to the oneProblem
 * structure is made in here. Both master problem right hand sides, master problem cut coefficients, and subprob right 
 * hand side need to be updated using the information calculated/stored beforehand.
 **********************************************************************************************************************/
int warmStartCell(cellType *cell, probType **prob, solnType *soln){
    
#ifdef  ROLL
    printf("\t~warmUpCell()\n");
#endif
    
    int     status;
    int     i, obs, cutID;
    int     lambdaIdx;
    int     masterRowCnt;
    int     effIn, effOut;
    double  warmupChange, tempChange, pibTrim;
    vector  masterRHS = NULL, reObs = NULL;
    intvec  indices = NULL;
    vector  *trimChange = NULL;
    
    /* Track the trimDelta changes and update the trimOmega section */
    if ( !(reObs = (vector) arr_alloc(cell->arima->subh * cell->arima->var + 1 , double)) )
        errMsg("allocation", "warmUpCell", "reObs", 0);
    if ( !(trimChange = (vector *) arr_alloc(cell->lambda->cnt, vector)) )
    	errMsg("allocation", "warmUpCell", "trimChange", 0);
    /* ***** Update the sigmaType with updated dynamics on the subprob right hand side ***** */
    // What is changing :: exo dynamics on subprob a[1] && endo dynamics on subprob B[1][1]
    if (prob[1]->aBar)
        for (i=0; i<cell->sigma->cnt; i++) {
#ifdef ROLL
            printf("\t\t :: (aDynLambda) Updating sigma->val[%d][Lambda%d].b = %lf", i, cell->sigma->lambdaIdx[i],cell->sigma->vals[i].b);
#endif
            warmupChange = vXv(cell->gamma->exo->aDelta, cell->gamma->exo->aDynLambda->vals[i], NULL, prob[1]->aBar->cnt);
            cell->sigma->vals[i].b += warmupChange;
#ifdef ROLL
            printf(" with value %lf \n", warmupChange);
#endif
        }
    
    if (prob[1]->Bbar[1])
        for (i=0; i<cell->sigma->cnt; i++) {
            warmupChange = vXv(cell->gamma->exo->BDelta, cell->gamma->exo->BDynLambda->vals[i], NULL, cell->gamma->exo->BCnt);
            cell->sigma->vals[i].b += warmupChange;
#ifdef ROLL
            printf("\t\t :: (BDynLambda) Updating sigma->val[%d][Lambda%d].b with value %lf \n", i, cell->sigma->lambdaIdx[i], warmupChange);
#endif
        }
    
    /* **** Warm-start trimDelta structure and track change for cuts updates **** */
    for (i=0; i<cell->lambda->cnt; i++)
        trimChange[i] = (vector) arr_alloc(cell->omega->cnt, double);
    
    for (obs=0; obs<cell->omega->cnt; obs++) {
        
        effIn = 0;
        effOut = 0;
        
        // Re-constructing observations from past epsilons
        for (i=1; i<=cell->arima->subh * cell->arima->var; i++) {
            // Calculate the re-simed observation
            reObs[i] = cell->omega->vals[obs][i] + cell->auxOmega->vals[obs][i] + cell->gamma->endo->bARMA->val[i];
            // Capture the negative trim here
            if (reObs[i] < cell->arima->lb) {
                reObs[i] = cell->arima->lb - reObs[i];
            } else {
                reObs[i] = 0.0;
            }
            // Update the trimOmega structure
            cell->trimOmega->vals[obs][i] = reObs[i];
        }

        // Calculate pi x trimObs for each lambda entry
        for (lambdaIdx=0; lambdaIdx<cell->lambda->cnt; lambdaIdx++) {
            // Calculation pi x trimObs using re-constructed obs
            pibTrim = vXv(cell->lambda->vals[lambdaIdx], reObs, NULL, prob[1]->num->numRV);
            
            if (pibTrim > run.EPSILON)
                effIn++;
            
            if (cell->trimDelta->val[lambdaIdx][obs].b > run.EPSILON)
                effOut++;
            
            // Track changes (ignoring C)
            trimChange[lambdaIdx][obs] = pibTrim - cell->trimDelta->val[lambdaIdx][obs].b;
            
            // Capture new vaue (ignoring C)
            cell->trimDelta->val[lambdaIdx][obs].b = pibTrim;
        }
        
#ifdef ROLL
        printf("\t :: Updating pi x b_trim [OBS %d] = Incurring %d || Resolving %d \n", obs, effIn, effOut);
#endif
    }
    
    /* Some cuts are dropped here base on height value */
    status = transformCuts(cell->master, cell, prob[0]->num->cols, cell->cuts, soln, cell->k);
    if (status) {
        errMsg("transform", "warmUpCell", "Failed to transform cuts", 0);
    }
    

    masterRowCnt = getNumRows(cell->master->lp);
    
#ifdef ROLL
    printf("\t :: Master row count = %d\n", masterRowCnt);
#endif
    
    /* Memory Allocation and Fetch Right Hand Side from Master Problem */
    if ( !(masterRHS = arr_alloc(masterRowCnt+1, double)) )
        errMsg("allocation","warmUpCell","Failed to allocate memory to masterRHS",0);
    if ( getRhsx(cell->master->lp, 0, masterRowCnt, masterRHS) ) {
        errMsg("solver","warmUpCell","Failed to reach current master problem right hand side",0);
        return 1;
    }
    if (!( indices = arr_alloc(cell->master->mar+cell->cuts->cnt+1,int) ))
        errMsg("allocation","warmUpCell","failed to allocate memory to indices",0);
    for (i=0; i<masterRowCnt; i++) indices[i] = i;

#ifdef ROLL
    printf("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n");
#endif
    
    /* **** Updates the RHS only on the surviving optimality Cuts **** */
    for (i=0; i<cell->cuts->cnt; i++) {

        cutID = i;

        /* Prepare the warm-up value :: each optimality is associated one warmup value */
        warmupChange = 0.0;
        
        /* Loop through all distinct observations assocaited with this cut */
        for (obs=0; obs<cell->cuts->val[i]->omegaCnt; obs++) {
            
            /* Fetch the index of a cut's defining dual multiplier for obs */
            lambdaIdx = cell->sigma->lambdaIdx[cell->cuts->val[i]->iStar[obs]];
            
            /* Calculate the change contribution from ARMA */
            tempChange = vXv(cell->lambda->vals[lambdaIdx], cell->gamma->endo->changebARMA->val, NULL, cell->gamma->endo->changebARMA->cnt);
            
            /* Calculate the change contribution from trimDelta */
            tempChange += trimChange[lambdaIdx][obs];
            
            /* Newly added: Calculate the change contribution from subprob dynamics in exougenous rhs */
            if (prob[1]->aBar)
                tempChange += vXv(cell->gamma->exo->aDynLambda->vals[cell->cuts->val[i]->iStar[obs]], cell->gamma->exo->aDelta, NULL, prob[1]->aBar->cnt);
            
            /* Newly added: Calculate the change contribution from subprob dynamic in endogenous decision */
            if (prob[1]->Bbar[1])
                tempChange += vXv(cell->gamma->exo->BDynLambda->vals[cell->cuts->val[i]->iStar[obs]], cell->gamma->exo->BDelta, NULL, cell->gamma->exo->BCnt);
            
            /* Weight the change given associated observation's weight */
            tempChange = ( (double) cell->omega->weight[obs] / (double) cell->k) * tempChange;
            
            /* Accumulate the changes of this observation */
            warmupChange += tempChange;
        }
        
#if defined(ROLL)
        printf("\t :: %d ROW=%d \t CUTOBS=%d \t CHANGE=%f \t AXbar=%f \t NEWAXbar=%f \t INCUMBEND = %d \n",
               i, cell->cuts->val[cutID]->rowNum, cell->cuts->val[cutID]->cutObs, warmupChange,
               cell->cuts->val[cutID]->alphaIncumb, cell->cuts->val[cutID]->alphaIncumb+warmupChange,
               cell->cuts->val[cutID]->isIncumb);
#endif
        
        /* Change the right hand side value in tempRHS */
        masterRHS[cell->cuts->val[i]->rowNum] += warmupChange;
        
        /* Update cuts alpha value and alphaIncumb as well */
        cell->cuts->val[i]->alpha += warmupChange;
        cell->cuts->val[i]->alphaIncumb += warmupChange;
    }

#ifdef ROLL
    printf("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n");
#endif
    
    /* Update the master problem with warm start cuts */
    status = changeRHS(cell->master->lp, prob[0]->num->rows+cell->cuts->cnt, indices, masterRHS);
    if (status) {
        errMsg("warmup","warmUpCell","Failed to update the warm-up right hand side",0);
        return 1;
    }
    
    cell->quadScalar = run.MIN_QUAD_SCALAR;
    
    if (masterRHS) mem_free(masterRHS);
    if (indices) mem_free(indices);
    if (reObs) mem_free(reObs);
    if (trimChange) {
        for (i=0; i<cell->lambda->cnt; i++)
            if (trimChange[i])
                mem_free(trimChange[i]);
        mem_free(trimChange);
    }
    
    printf("Warm start cellStructure ... DONE\n");
    
    return 0;
}//END of warmUpCell()


/***********************************************************************************************************************
 * This subroutine is designed to refresh some of the information in solnType when getting a cell ready for the next
 * time period. Since solnType carries many useful information, we don't want to break and reconstruct again as we
 * warm up. Also, some of the informaiton, mainly flags, will affect how the cell is going to proceed for the next
 * time horizon. Hence, we only neet to update a few members in this structure.
 * This subroutine is only used when we try to solve the problem stochastically.
 
 
 // There are still needs to full go over this subroutine.
 
 
 ***********************************************************************************************************************/
int refreshSoln(solnType *soln, cellType *cell, probType **prob) {
    
    int i;
    double argmaxPre, argmaxPost;
    double argmaxPreSum, argmaxAllSum;
    vector piCbarX;
    
    /*****************************************************************************************************
     * Pi ratio evaluation is with a certain scan length. When running regular rSD, it is always assumed
     * this vector need to be cleared out. However, under rSD with warm start condition, it is okay for
     * user to say don't refresh the pi ratio vector. Then the pi ratio is reconstructed base on the
     * current argmax condition.
     *****************************************************************************************************/
    
    printf("Refreshing solution structure ... \n");
    if (run.PI_EVAL_CLEAR == 1 || run.WARM_UP == 0) {
        mem_free(soln->piSRatio);
        if (!(soln->piSRatio = (vector) arr_alloc(run.SCAN_LEN, double)))
            errMsg("alloation", "cleanUpSoln", "piSRatio", 0);
    }else{
        if ( !(piCbarX= arr_alloc(cell->sigma->cnt, double)) )
            errMsg("Allocation", "refreshSoln", "piCbarX",0);
        for (i = 0; i < run.SCAN_LEN; i++) {
            
            /* Initialize the height indicator */
            argmaxPre = -DBL_MAX;
            argmaxPost = -DBL_MAX;
            argmaxPreSum = 0;
            argmaxAllSum = 0;
            
            /* Search within the observation range (BEGINNING) -> (RATIO MOMENT) */
            r_computeIstar(prob[1]->num, prob[1]->coord, cell->arima, cell->omega, cell->lambda,
                           cell->gamma, cell->sigma, cell->delta, cell->auxDelta, cell->trimDelta,
                           soln->incumbX, piCbarX, &argmaxPre, &argmaxPost, cell->omega->cnt - run.SCAN_LEN + i, cell->k, cell->k - cell->Kt);
            
            argmaxPreSum += max((argmaxPre - 0.0), 0) * cell->omega->weight[cell->omega->cnt - run.SCAN_LEN + i];
            argmaxAllSum += max((max(argmaxPre, argmaxPost) - 0.0), 0.0) * cell->omega->weight[cell->omega->cnt - run.SCAN_LEN + i];
            
            evalDualStability(soln->piSRatio, argmaxPreSum, argmaxAllSum, &soln->dualStableFlag, cell->k - (run.SCAN_LEN - i));
        }
        mem_free(piCbarX);
    }
    
    /* General parameter updates for directing the code in next period */
    soln->gamma = 0.0;
    soln->optFlag = FALSE;
    soln->normDk = 0.0;
    soln->normDk_1 = 0.0;
    soln->subFeasFlag = TRUE;
    soln->optMode = TRUE;
    soln->dualStableFlag = FALSE;
    soln->preCheckFlag = FALSE;
    soln->preCheckEverFlag = FALSE;
    soln->FTError = 0.0;
    soln->repPassed = 0;
    soln->iCutUpdt = cell->k;

    
    /* Diversed scheme base on different running configurations */
    if (run.WARM_UP == 1) {
        soln->incumbEst = maxCutHeight(cell->lbType, cell->lb, cell->cuts, cell->k, soln->incumbX, prob[0]->num->cols, &soln->iCutIdx);
        soln->incumbEst += vXvSparse(soln->incumbX,prob[0]->dBar);
        soln->candidEst = soln->incumbEst;
        soln->cCutIdx = cell->k;
    }else{
        soln->candidEst = cell->lb + vXvSparse(soln->incumbX, prob[0]->dBar);
        soln->incumbEst = soln->candidEst;
        for (i=1; i<=prob[0]->num->cols; i++)
            soln->candidX[i] = soln->incumbX[i];
        soln->iCutIdx = 0;
        soln->cCutIdx = 0;
    }
    
    printf(" ... DONE\n");
    return 0;
}

/* This a auxilary function is used for sorting cuts and arrange their index base on their cuts */
static int compar (const void *a, const void *b){
    int aa = *((int *) a), bb = *((int *) b);
    if (base[aa] < base[bb])
        return -1;
    if (base[aa] == base[bb])
        return 0;
    if (base[aa] > base[bb])
        return 1;
    return 0;
}

/*************************************************************************************************************
 Between two time periods, the warm-up operation focuses on reforming all the cuts in algorithm cell. This 
 may not be necessary since some of the cuts currently in cell is with loose lower bound. It is possible that
 some of the loose cuts may be warmed to be a better one. But it is not necessary to keep all of them since
 the potentially beneficial cut can be regenerated pretty fast in the next period. Hence, we only keep a
 portion of the current cuts for warming up rather than all the cuts. This subroutine first check the
 configuration input and calculate exactly how many cuts should we keep/drop. Then the cuts are dropped base
 on their height. Since the problem is always assume to be an minimize problem, the loose cust are considered
 to be the ones with less height at incumbent solution. The cuts' height value are maintained in the same 
 formation as cuts are dropped to keep operation correct.
 *************************************************************************************************************/
int transformCuts(oneProblem *master, cellType *cell, int mCols, cutType *cuts, solnType *soln, int k) {
    
#ifdef ROLL
    printf("\t\t~transformCuts()\n");
#endif
    
    vector  heights = NULL;
    intvec  sortedCutIdx = NULL;
    int     i,j;
    int     lowCutIdx, warmCutCnt, coldCutCnt, initCutCnt, lastCutIdx;
    
    /* Keep a count number of the total cuts before transformation */
    initCutCnt = cuts->cnt;
    
    /* Allocate memory to a vector that stores all cut's height */
    if (!(heights = (vector) arr_alloc(cuts->cnt, double)))
        errMsg("allocation", "transformCuts", "heights", 0);
    if (!(sortedCutIdx = (intvec) arr_alloc(cuts->cnt, int)))
        errMsg("allocation", "transformCuts", "cutIdx",0);
    
    /* Prepare a vector of every cut's height value */
    for (i=0; i<cuts->cnt; i++) {
        heights[i] = cutHeight(cell->lbType, cell->lb, cell->cuts->val[i], k, soln->incumbX, mCols);
        sortedCutIdx[i] = i;
    }
    
    /* Give augmented array to global sorting base */
    base = heights;
    
    /* Sort the heights values with index kept in place */
    qsort(sortedCutIdx, cuts->cnt, sizeof(int), compar);
    
    /* Calculate the # of cuts that need to be transformed by selecting the maximum height */
    warmCutCnt = (int)(cuts->cnt * run.PERCENT_WARM);
    coldCutCnt = cuts->cnt - warmCutCnt - 1;
    
#ifdef ROLL
    printf("\t\t :: Sorting cuts' heights... Prepare to drop %d cuts [TOTAL %d]\n",coldCutCnt, initCutCnt);
#endif
    
    /* Once the total amount of cuts that need to be dropped is known, drop these cuts */
    for (i=0; i<coldCutCnt; i++) {
        /* Find the lowest cut thourgh searching the active area */
#ifdef ROLL
        printf("\t\t :: Dropping Lowest Cut [Idx] = %d with height %f\n", sortedCutIdx[i], heights[sortedCutIdx[i]]);
#endif
        lowCutIdx = sortedCutIdx[i];
        /* Keep the last cut index for updating cutIdx */
        lastCutIdx = cuts->cnt-1;
        /* Drop the lowest cut */
        dropCut(master->lp, cuts, soln, lowCutIdx);
        /* Update the cutIdx accordingly */
        for (j=i+1; j<initCutCnt; j++) {
            if (sortedCutIdx[j] == lastCutIdx) {
#ifdef ROLL
                printf("\t\t :: Now cut %d has been moved to position %d\n", lastCutIdx,lowCutIdx);
#endif
                sortedCutIdx[j] = lowCutIdx;
                continue;
            }
        }
    }
    
    mem_free(heights); mem_free(sortedCutIdx);
    return 0;
}
