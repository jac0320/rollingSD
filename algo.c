//
//  algo.c
//  algo
//
//  Created by Site Wang on 2/11/16.
//  Copyright © 2016 Site Wang. All rights reserved.
//

#include "rollingSP.h"
#include <time.h>
#include "cell.h"
#include "output.h"
#include "soln.h"
#include "roll.h"
#include "dynamic.h"
#include "algo.h"
#include "evaluate.h"

extern runType      run;
extern string       outputDir;
extern char         inputDir;
extern BOOL         newPeriod;

/****************************************************************************************************************
 This is where everything is. Main algorithm is divided into two major branches.
    1 - Deterministic Solving Rolling Horizon Problem
    2 - Stochastically Solving Rolling Horizon Problem
        |__With Warm-Up Scheme for Optimality Cuts
        |__Without Warm-Up Scheme for Optimality Cuts
 ****************************************************************************************************************/
int algo(oneProblem *orig, stocType *stoc, timeType *tim){
    
#ifdef  TRACE
    printf("~algo()\n");
#endif
    
	probType **prob = NULL;
	cellType *cell = NULL;
	solnType *soln = NULL;
    BOOL     initialPeriod = TRUE;
    int      status = 0;
    clock_t  tic, toc;
    
    FILE *objOUT = NULL;
    FILE *solOUT = NULL;
    FILE *reportOUT = NULL;
    FILE *SubOptVec = NULL;
    FILE *ePtr = NULL;
    
    /* Public File Pointer for both Stochastic and Deterministic */
    SubOptVec = openFile(outputDir, "subOptVec.txt", "w");
    ePtr = openFile(outputDir, "eval.txt","w");

    /* Setup Algorithm */
	cell = setup_rSP(orig, stoc, tim, &prob, &soln);
	if (cell == NULL) {
		errMsg("Setup", "algo", "Failed to construct algorithm cell", 0);
		goto STOPALGO;
	}
    
    /* Print Out Initial Information of the Problem */
    printInitInfo(orig, tim, prob, cell);
    
    if (run.DETERMINISTIC == 1) {
        
        objOUT = openFile(outputDir, "DET_obj.out", "w");
        solOUT = openFile(outputDir, "DET_sol.out", "w");

        /*------------- MAJOR BRANCH :: Problem is solved deterministically in here -------------*/
        char    outputDetMasterDir[BLOCKSIZE];
        strcpy(outputDetMasterDir, outputDir);
        createOutputDir(outputDetMasterDir, "", "detProblem");

        while ( cell->t < run.HORIZON ) {
            
            tic = clock();
            if (run.HORIZON > 1) {
                newPeriod = TRUE;
                /****** Warm Start for Deterministic Problem ******/
                status = setupDetPeriod(cell->master, cell, prob, soln, &initialPeriod);
                if (status) {
                    errMsg("algorithm", "algo", "failed to update the deterministic cell for next time period", 0);
                    return 1;
                }
            } else {
                printf("Solving deterministic cell, cell->t = 1; \n");
                cell->t = 1;
            }
            toc = clock();
            soln->runTime->warmUpTime[cell->t] = ((double) (toc - tic)) / CLOCKS_PER_SEC;
    
            
            tic = clock();
            /* Solve a deterministic cell and collect solution */
            if ( solveDETCell(cell, prob[0], soln) ) {
                errMsg("algorithm", "algo", "failed to solve deterministic cell", 0);
            }
            toc = clock();
            /* Record Overall Running Time */
            soln->runTime->solveTime[cell->t] = ((double) (toc - tic)) / CLOCKS_PER_SEC;
            
            if (run.EVALUATOR > 0) {
                tic = clock();
                /* Read the solution file */
                if (injectSolution(cell, prob[0], soln))
                    errMsg("algorithm", "algo", "failed to inject solution to deterministic cell", 0);
            }
            
            /* Output Periodic Problem Information */
            printPeriodStat(prob, soln, cell, objOUT, solOUT, NULL);
            
            /* Evaluate deterministic solution */
            if ( run.EVAL_FLAG )
            	evalDetOpt(orig, tim, stoc, prob, cell, soln, SubOptVec, ePtr);
            
            if (run.EVALUATOR > 0)
                break;
            
        }
        printLine();
        /*************************************************************************/
        
    } else {
        
        /*--------------- MAJOR BRANCH :: Problem is solved stochastically in here --------------*/
        objOUT = openFile(outputDir, "STO_obj.out", "w");
        solOUT = openFile(outputDir, "STO_sol.out", "w");
        reportOUT = openFile(outputDir, "STO_report.out", "w");

#ifdef MODIFY
        /* ~~~~~~~~~~~~ Prepare Debug Output File Directory ~~~~~~~~~~~~ */
        char    outputMasterDir[BLOCKSIZE];
        char    outputSubprobDir[BLOCKSIZE];
        char    outputIncumbMasterDir[BLOCKSIZE];
        strcpy(outputMasterDir, outputDir);
        strcpy(outputSubprobDir, outputDir);
        strcpy(outputIncumbMasterDir, outputDir);
        createOutputDir(outputMasterDir, "", "master");
        createOutputDir(outputSubprobDir, "", "subprob");
        createOutputDir(outputIncumbMasterDir, "", "incumbMaster");
#endif
        
#ifdef ROLL
        char    outputRollDir[BLOCKSIZE];
        strcpy(outputRollDir, outputDir);
        createOutputDir(outputRollDir, "", "roll");
#endif
        
        /* run through all the time periods in the horizon */
        while ( cell->t < run.HORIZON ) {
            tic = clock();
            // Rolling horizon is only conducted when dealing with such a problem
            if (run.HORIZON > 1) {
                /* Activate this flag for rolling horizon problem */
                newPeriod = TRUE;
                /******* Warm Start :: SETUP NEXT TIME PERIOD PROBLEM & ALGORIHMIC CELL *******/
                if (run.WARM_UP == 1) {
                    // Setup the next time perid s
                    status = setupStoPeriod(cell, stoc, prob, soln, &initialPeriod);
                    if ( status ) {
                        errMsg("algorithm", "algo","failed to update the cell for next time period",0);
                        return 1;
                    }
                } else {
                    status = refreshPeriod(cell, stoc, prob, soln, &initialPeriod);
                    if ( status ) {
                        errMsg("algorithm", "algo", "failed to update the cell for next time period", 0);
                    }
                }
            } else {
                printf("Solving using regular SD, cell->t = 1; \n");
                cell->t = 1;
            }
            toc = clock();
            soln->runTime->warmUpTime[cell->t] = ((double) (toc-tic)) / CLOCKS_PER_SEC;
            
            /* Solve Cell for the Current Time Period */
            tic = clock();
            if ( solveSDCell(stoc, prob, cell, soln) ) {
                errMsg("algorithm", "algo", "failed to solve cell", 0);
                return 1;
            }
            toc = clock();
            soln->runTime->solveTime[cell->t] = ((double) (toc-tic)) / CLOCKS_PER_SEC;
            
            /* Print Periodic Solution */
            printPeriodStat(prob, soln, cell, objOUT, solOUT, reportOUT);
            
            /* Evaluate stochastic solution */
            if ( run.EVAL_FLAG )
            	evalStoOpt(stoc, prob, cell, soln, SubOptVec, ePtr);

        }
        printLine();
        /*************************************************************************/
        
        fclose(reportOUT);

    }
    
    /* Collect Running Time by Printing them */
    printRunningTime(soln);
    
    fclose(ePtr);
    fclose(objOUT); fclose(solOUT); fclose(SubOptVec);
    freeSolnType(soln);
	freeCellType(cell);
    freeProbType(prob, tim->numStages);
	return 0;

STOPALGO:
	fclose(ePtr); fclose(reportOUT);
    fclose(objOUT); fclose(solOUT); fclose(SubOptVec);
    freeSolnType(soln);
	freeCellType(cell);
    freeProbType(prob, tim->numStages);
	return 1;
}//END algo()



cellType *setup_rSP(oneProblem *orig, stocType *stoc, timeType *tim, probType ***prob, solnType **soln) {

    cellType    *cell = NULL;
    vector      xk = NULL, lb = NULL;
    BOOL        arimaFlag = FALSE;
    int         status = 0, stat1, usefulProbIdx = 0, masterProbIdx = 0, subprobIdx = 0;
    
    /* ----- Solve Mean Problem ----- */
    xk = meanProblem(orig, stoc);
    if (xk == NULL) {
        errMsg("solve","setup_rSP","Failed to solve the mean problem.",0);
        return NULL;
    }
    
    /* ----- Get Mean Problem Solution ---- */
    lb = calcLowerBound(orig, tim);
    if (lb == NULL) {
        errMsg("Solving","setup_rSP","Failed to obtain the mean problem solution",0);
        return NULL;
    }
    
    
    /******************************************************************************************
     
                                    rhs
     -----------------------     ----------  <---- MasterStart              -|
       First Stage Problem          DYN      <---- Actual Master Problem     | Deterministic
     -----------------------     ----------  <---- SubprobStart              | Problem
       Second Stage Problem       DYN + TS   <---- Actual Subprob           _|
     -----------------------     ----------
     
     ******************************************************************************************/
    
    /* --------------------- Decomposed Problem Navigation --------------------- */
    if (run.DETERMINISTIC == 1) {
        tim->numStages--;
        usefulProbIdx   = 0;  //The deterministic problem is stored in prob[0];
        masterProbIdx   = 0;
        subprobIdx      = 0; // This is disabled
    } else {
        usefulProbIdx   = 1; // The recourse problem is store in prob[1]
        masterProbIdx   = 0; // Master Problem
        subprobIdx      = 1; // SubProblem
    }

    /* ------------------------- Problem Decomposition ------------------------- */
    (*prob) = newProb(orig, stoc, tim, lb, 0.1);
    if (prob == NULL) {
        errMsg("decompose","setup_rSP","Failed to decompose the problem",0);
        goto SETUP_FAIL;
    }
    
    
    /* ----- Construct Algorithm Cell ----- */
    cell = newCell((*prob), stoc);
    if (cell == NULL) {
        errMsg("Setup", "setup_rSP", "Failed to construct algorithm cell", 0);
        goto SETUP_FAIL;
    }

    /******************************************************************************************
     Assumed that run.HORIZON will be greater than 1 if deallsing with rolling horizon problem. 
     When dealing with rolling Horizon problems, time series model is required. As we are testing
     rolling horizon problem on ENERGY application, we assume that a powerCurve is need for 
     translation between wind speed and wind power genreation.
     ******************************************************************************************/
    if (run.HORIZON > 1) {
        /* Read Time Series Model */ /* Useful index is assigned during decomposition */
        cell->arima = readARIMA(&inputDir, orig->name, (*prob)[usefulProbIdx]->sp, &arimaFlag, &stat1);
        if ( cell->arima == NULL && stat1 == 1 ) {
            errMsg("read","main","Failed to read arma file",0);
            goto SETUP_FAIL;
        } else {
            // Indicate in stoc type that the stoc type is ARIMA
            printf(".arma file detected and read. Changing stoc type to ARIMA\n");
            strcpy(stoc->type, "ARIMA");
        }
        
        /* Read Dynamic Information and Initial File */
        status = readDyn
        (&inputDir, orig->name, cell, (*prob), tim->numStages);
        if (status) {
            errMsg("read", "setup_rSP", "Failed to read dynmic file", 0);
            goto SETUP_FAIL;
        } else {
            printf(".dyn file detected and read.\n");
        }
    }
    
    (*soln) = newSoln((*prob), cell->maxCuts, xk, lb[0], usefulProbIdx);
    if ((*soln) == NULL) {
        errMsg("setup", "setup_rSP", "failed to setup the solution used by the algorithm", 0);
        goto SETUP_FAIL;
    }

    if (xk) mem_free(xk);
    if (lb) mem_free(lb);
    return cell;

SETUP_FAIL:
    if (xk) mem_free(xk);
    if (lb) mem_free(lb);
    return NULL;
    
}//END setup_rSP()



void freeSingleProbtype (probType *prob, int stageCnt) {
    int i;
    
    for (i=0; i<stageCnt; i++) {
        if (prob->Bbar[i]) freeSparseMatrix(prob->Bbar[i]);
    }
    if ( prob->Bbar ) mem_free(prob->Bbar);
    if ( prob->Abar ) freeSparseMatrix(prob->Abar);
    if ( prob->Cbar ) freeSparseMatrix(prob->Cbar);
    if ( prob->Dbar ) freeSparseMatrix(prob->Dbar);
    if ( prob->aBar ) freeSparseVector(prob->aBar);
    if ( prob->bBar ) freeSparseVector(prob->bBar);
    if ( prob->cBar ) freeSparseVector(prob->cBar);
    if ( prob->dBar ) freeSparseVector(prob->dBar);
    if ( prob->sp )   freeOneProblem(prob->sp);
    if ( prob->name)  mem_free(prob->name);
    if ( prob->num )  mem_free(prob->num);
    if ( prob->coord) freeCoordType(prob->coord);
    if ( prob->omegas)freeOmegastuff(prob->omegas);
    mem_free(prob);
}
