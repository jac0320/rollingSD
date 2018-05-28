//
//  out.c
//  rollingSP
//
//


#include "output.h"

extern runType run;
extern string outputDir;

/*******************************************************************************
 This subroutine is designed to write output information into a file. It is 
 currently in its initial shape. More modifications will be made on this later.
 *******************************************************************************/
void printPeriodStat (probType **prob, solnType *soln, cellType *cell, FILE *objV, FILE *solution, FILE *report) {

    if (run.DETERMINISTIC == 1) {
        printf(" :: Time Period \t %d \n", cell->t);
        printf(" :: Deterministic Problem Solution \t %f \n",soln->optValM);
        printf(" :: Solution one-Norm = %f \n :: ",soln->incumbX[0]);
    } else {
        printf(" :: Time Period   \t %d   \n", cell->t);
        printf(" :: Optimization Running Time \t %lf\n", soln->runTime->solveTime[cell->t]);
        printf(" :: Incumb Est    \t %lf  \n", soln->incumbEst);
        if (run.WARM_UP == 1) {
            printf(" :: Periodic Run Iterations = %d, Cumulative Iteration = %d \n", cell->k - cell->Kc, cell->k);
        } else {
            printf(" :: Periodic Run Iterations = %d, Cumulative Iteration = %d \n", cell->k, cell->Kc + cell->k);
        }
        printf(" :: Solution one-Norm = %f \n ::",soln->incumbX[0]);
    }
    
    if (objV != NULL) {
        if (run.DETERMINISTIC == 1) {
            fprintf(objV, "%f\n",soln->optValM);
            printVector(soln->incumbX-1, prob[0]->sp->mac+1, solution);
        } else {
            fprintf(objV, "%f\n",soln->incumbEst);
            printVector(soln->incumbX-1, prob[0]->sp->mac+1, solution);
            if (cell->t == 1) {
                fprintf(report, "t \t k \t omegaCnt \t lambdaCnt \t sigmaCnt \t argmaxCnt \t LPCnt \n");
            }
            fprintf(report, "%d \t %d \t %d \t %d \t %d \t %d \t %d \n",
                    cell->t, cell->k, cell->omega->cnt, cell->lambda->cnt, cell->sigma->cnt
                    ,soln->runTime->argmaxCnt, cell->LPcnt);
        }
    }
}

void printInitInfo(oneProblem *orig, timeType *tim, probType **prob, cellType *cell) {
    
    int i;
    
    FILE *logF;
    logF = openFile(outputDir, "log.txt", "w");
    
    /* Print Initial Problem Information onto the Screen */
    printf("\n\n########################### Starting Rolling Horizon Algorithm ############################\n");
    fprintf(logF, "\n\n########################### Starting Rolling Horizon Algorithm ############################\n");
    printf("Problem Name        :: %s \n", orig->name);
    fprintf(logF, "Problem Name        :: %s \n", orig->name);
    printf("Horizon             :: %d \n", run.HORIZON);
    fprintf(logF, "Horizon             :: %d \n", run.HORIZON);
    printf("Original Problem    :: Cols = %d; Rows = %d;\n", orig->mac, orig->mar);
    fprintf(logF, "Original Problem    :: Cols = %d; Rows = %d;\n", orig->mac, orig->mar);
    printf("Decomposed Stages   :: %d \n", tim->numStages);
    fprintf(logF, "Decomposed Stages   :: %d \n", tim->numStages);
    for ( i=0; i<tim->numStages; i++) {
        if (prob[i]) {
            // A double check here is conducted
            if (prob[i]->sp) {
                printf("\t Stage %d :: %d Cols, %d Rows, %d RVs \n",
                        i+1, prob[i]->num->cols, prob[i]->num->rows, prob[i]->num->numRV);
                fprintf(logF, "\t Stage %d :: %d Cols, %d Rows, %d RVs \n",
                       i+1, prob[i]->num->cols, prob[i]->num->rows, prob[i]->num->numRV);
            }
        }
    }
    if (cell->master->type == PROB_LP) {
        printf("Master Type         :: LP \n");
        fprintf(logF, "Master Type         :: LP \n");
    } else if (cell->master->type == PROB_QP) {
        printf("Master Type         :: QP \n");
        fprintf(logF, "Master Type         :: QP \n");
    } else {
        printf("Master Type         :: Unknown %d \n", cell->master->type);
        fprintf(logF, "Master Type         :: Unknown %d \n", cell->master->type);
    }
    if(run.DETERMINISTIC == 1) {
        printf("Solving Problem Deterministically...\n");
        fprintf(logF, "Solving Problem Deterministically...\n");
    } else {
        printf("Solving Problem Stochasticallly...\n");
        fprintf(logF, "Solving Problem Stochasticallly...\n");
        if (run.WARM_UP == 1) {
            printf("Warm Start Scheme is ON...\n");
            fprintf(logF, "Warm Start Scheme is ON...\n");
            printf("Warm Start %.0f %% cuts...\n", run.PERCENT_WARM * 100);
            fprintf(logF, "Warm Start %.0f %% cuts...\n", run.PERCENT_WARM * 100);
            printf("EPSILON value is %lf\n", run.EPSILON);
            fprintf(logF, "EPSILON value is %lf\n", run.EPSILON);
            printf("RUN SEED == %llu\n", run.RUN_SEED);
            fprintf(logF, "RUN SEED == %llu\n", run.RUN_SEED);
            printf("Evaluation SEED == %llu\n", run.RESAMPLE_SEED);
            fprintf(logF, "Evaluation SEED == %llu\n", run.RESAMPLE_SEED);
            printf("Full Check SEED == %llu\n", run.EVAL_SEED);
            fprintf(logF, "Full Check SEED == %llu\n", run.EVAL_SEED);
            printf("Scalar of RV = %lf\n",run.SCALAR);
            fprintf(logF, "Scalar of RV = %lf\n", run.SCALAR);
            if (run.PI_EVAL_CLEAR == 1) {
                printf("Clearing piRatios... Minimum Run Iteration is Automatically set to SCAN_LENGTH%d\n", run.SCAN_LEN);
                printf("Allowed error level is %lf\n", run.DUAL_EPSILON);
                fprintf(logF, "Clearing piRatios... Minimum Run Iteration is Automatically set to SCAN_LENGTH%d\n", run.SCAN_LEN);
                fprintf(logF, "Allowed error level is %lf\n",run.DUAL_EPSILON);
            } else {
                printf("Warm start piRatios as well... Minimum run Length could be less than SCAN_LENGTH[%d]\n", run.SCAN_LEN);
                fprintf(logF, "Warm start piRatios as well... Minimum run Length could be less than SCAN_LENGTH[%d]\n", run.SCAN_LEN);
            }
        } else {
            printf("Warm Start Scheme is OFF...\n");
            fprintf(logF, "Warm Start Scheme is OFF...\n");
        }
    }
    
    fclose(logF);
    printLine();
}

/* Customized Utility Printing Subroutine */
void rSPprintf(FILE *out, const char *s) {
    
    if (out == NULL) {
        printf("%s",s);
    } else {
        fprintf(out, "%s", s);
    }
}

/***********************************************************************\
 ** This function prints the relevant information in a cut.
 ** It is meant to be used for debugging.
 \***********************************************************************/
void printCut(cutType *cuts, numType *num, int idx) {
    
    int cnt;
    
    printf("\t\t :: Cut #%d:: c:%d o:%d \t\t :: a:%f B:",
           idx, cuts->val[idx]->cutObs, cuts->val[idx]->omegaCnt, cuts->val[idx]->alpha);
    for (cnt = 0; cnt <= num->prevCols; cnt++)
        printf("%f ", cuts->val[idx]->beta[cnt]);
    printf("\nistar: ");
    for (cnt = 0; cnt < cuts->val[idx]->omegaCnt; cnt++)
        printf("%d ", cuts->val[idx]->iStar[cnt]);
    printf("\n");
   
}

void printRunningTime(solnType *soln) {
    
    FILE *pT, *iterT, *cutT;
    int i;

    pT = openFile(outputDir, "PERIOD_runtime.txt", "w");
    iterT = openFile(outputDir, "ITER_runtime.txt", "w");
    cutT = openFile(outputDir, "CUT_runtime.txt","w");
    
    /* Periodic Level Running Time Output */
    fprintf(pT, "Period -> %d\n", run.HORIZON);
    fprintf(pT, "periodTime \t warmUpTime \t evalOptTime \n");
    for (i=1; i<=run.HORIZON; i++) {
        fprintf(pT, "%lf \t %lf \t %lf \n",
                soln->runTime->solveTime[i], soln->runTime->warmUpTime[i], soln->runTime->evalOptTime[i]);
    }
    
    /* Iteration Level Running Time Output */
    fprintf(iterT, "Iter -> %d\n", soln->runTime->masterCnt);
    fprintf(iterT, "t \t masterTime \t optTime \t fC_sample \t fC_reform \t fC_lb \n");
    for (i=1; i<=soln->runTime->masterCnt; i++) {
        fprintf(iterT, "%d \t%lf \t %lf \t %lf \t %lf \t %lf\n",
                soln->runTime->iterIndex[i], soln->runTime->masterTime[i], soln->runTime->optTime[i],
                soln->runTime->fullCheck_sample[i], soln->runTime->fullCheck_reform[i], soln->runTime->fullCheck_lb[i]);
    }
    
    /* Solving Level Running Time Output */
    fprintf(cutT, "SubLPCnt -> %d\n", soln->runTime->subprobCnt);
    fprintf(cutT, "t \t subprobTime \t cutGenTime \t argmaxTime \t stoStrucTime \n");
    for (i=1; i<=soln->runTime->subprobCnt; i++) {
        fprintf(cutT, "%d \t %lf \t %lf \t %lf \t %lf \n",
                soln->runTime->lpCntIndex[i], soln->runTime->subprobTime[i], soln->runTime->cutGenTime[i], soln->runTime->argmaxTime[i], soln->runTime->stoStrucTime[i]);
    }
    
    fclose(pT); fclose(iterT); fclose(cutT);
}
