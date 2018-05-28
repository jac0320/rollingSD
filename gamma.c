//
//  gamma.c
//  rollingSP
//
//  Created by Site Wang on 4/9/16.
//  Copyright © 2016 Site Wang. All rights reserved.
//

#include "rollingSP.h"
#include "gamma.h"

extern runType run;

/* This subroutine updates the gammaType structure when a new dual multiplier is detected. There are two main things to update:
 * updating the stochastic information and updating the dynamics components. To update the stochastic information, we need to
 * calculate the pi x bAR/bMA, where bAR/bMA are the deterministic components from the time series where historical information
 * is used. In the mean time, we update the dynamics components by storing the dual solution that is associated with the dynamics
 * rows. The calculation using these information is later conducted when warming up the cell for the next period.
 * One need to note that the dual solution is with oneNorm and associated row coordination in gammaType for dynamics is indexed from
 * 0. Hence, when storing coressponding dual solution, we need +1 to seek for the right dual solution we need here. The information
 * stored start reading from index 1.
 */
int calcGamma(int lambdaIdx, int iter, gammaType *gamma, sparseVector *aBar, sparseMatrix **Bbar, vector pi, BOOL newSigmaFlag, BOOL newLambdaFlag) {
#ifdef TRACE
    printf("\t\t\t~CalcGamma()\n");
#endif
    int i;
    double pibARMA;
    
    if (run.HORIZON == 1) {
        /* Don't Do Anything here when no gamma is required */
        pibARMA = 0.0;
    } else {
        // Dynamic structure is synced with simgaType
        if (newSigmaFlag == TRUE) {
            //Dual solutions starts from one and the index is originated from 0
            /* Record the dual solutions that is associated with the subprob dynamics RHS */
            if (aBar) {
                /* Dynamics part: Allocate memory for dual multiplier storage */
                if ( !(gamma->exo->aDynLambda->vals[gamma->exo->aDynLambda->cnt] = (vector) arr_alloc(gamma->exo->aCnt+1, double)))
                    errMsg("allocation", "calcGamma", "new vals for dynLambda", 0);
                for (i=1; i<=gamma->exo->aCnt; i++)
                    gamma->exo->aDynLambda->vals[gamma->exo->aDynLambda->cnt][i] = pi[aBar->col[i]];
                /* A new entry in the gamma type for exo dynamics part */
                gamma->exo->aDynLambda->cnt++;
            }
            
            /* Record the dual solutions that is assocaited with the subprob dynamic decisions */
            if (Bbar[1]) {
                /* Dynamics part: Allocate memory for dual multiplier storage */
                if ( !(gamma->exo->BDynLambda->vals[gamma->exo->BDynLambda->cnt] = (vector) arr_alloc(gamma->exo->BCnt+1, double)))
                    errMsg("allocation", "calcGamma", "new vals for dynLambda", 0);
                for (i=1; i<=gamma->exo->BCnt; i++)
                    gamma->exo->BDynLambda->vals[gamma->exo->BDynLambda->cnt][i] = pi[Bbar[1]->row[i]];
                gamma->exo->BDynLambda->cnt++;
            }
        }
        
        // Time-series Stochastic structure is synced with lambdaType
        if (newLambdaFlag == TRUE) {
            /* Calculate Pi * gamma.AR & Pi gamma.MA */
            pibARMA = vXvSparse(pi, gamma->endo->bARMA);
            
            /* Store pi x bAR/bMA into the gamma structure  */
            gamma->endo->pibARMA[gamma->endo->pibCnt] = pibARMA;
            gamma->endo->pibARMA_ck[gamma->endo->pibCnt] = iter;
            gamma->endo->pibCnt++;
        }
        
    }
    
    /* This flag is not used but we are keeping it for no good reason, just being careful ... */
    return gamma->endo->pibCnt - 1;
}

/*
 This subroutine generates a new gamma structure in cell that stores all related dynamic and warm start information
 */
gammaType *newGamma(int numPeriod, int masterCols, int subprobCols, int subbBarLength, int lambdaCnt) {
    
    gammaType *gamma = NULL;
    
    if ( !(gamma = (gammaType *) mem_malloc(sizeof(gammaType))) ) {
        errMsg("allocation","newGamma","gamma",0);
        return NULL;
    }
    
    /* Main member of gammaType 
        exo :: exogenous stochasticity (supported by the external data)
        endo :: endogenous stochasticity (time-series related and dynamic components)
     */
    gamma->exo = NULL;
    gamma->endo = NULL;
    
    // Allocate Memory to gamma->exdo
    if ( !(gamma->exo = (exoType *) mem_malloc(sizeof(exoType))) ) {
        errMsg("allocation", "newGamma", "gamma->exdo", 0);
        return NULL;
    }
    gamma->exo->aData = NULL;        // Allocation is conducted when reading dynamic initial file
    gamma->exo->aCnt = 0;            // Assignment is conducted when reading dynamic initial file
    gamma->exo->BCnt = 0;            // Assignment is conducted when reading dynamic initial file
    gamma->exo->aDelta = NULL;       // Improve this part for later
    gamma->exo->BDelta = NULL;       // Improve this part for later
    gamma->exo->aDynLambda = NULL;   //
    gamma->exo->BDynLambda = NULL;   //
    
    if ( !(gamma->exo->aDynLambda = (dynLambdaType *) mem_malloc(sizeof(dynLambdaType)))) {
        errMsg("allocation", "newGamma", "aDynLambda", 0);
        return NULL;
    }
    if (!(gamma->exo->BDynLambda = (dynLambdaType *) mem_malloc(sizeof(dynLambdaType)))) {
        errMsg("allocation", "newGamma", "BDynLambda" , 0);
        return NULL;
    }
    
    gamma->exo->aDynLambda->cnt = 0;
    if ( !(gamma->exo->aDynLambda->vals = (vector *) arr_alloc(lambdaCnt, vector)) ) {
        errMsg("allocation", "newGamma", "aDynLambda->vals", 0);
        return NULL;
    }
    gamma->exo->BDynLambda->cnt = 0;
    if ( !(gamma->exo->BDynLambda->vals = (vector *) arr_alloc(lambdaCnt, vector)) ){
        errMsg("allocation", "newGamma", "BDynLambda->vals", 0);
        return NULL;
    }
    
    // Allocate Memory to gamma->endo
    if ( !(gamma->endo = (endoType *) mem_malloc(sizeof(endoType))) ) {
        errMsg("allocation","newGamma","gamma->indo",0);
        return NULL;
    }
    gamma->endo->pibCnt = 0;
    gamma->endo->bARMA = NULL;
    gamma->endo->changebARMA = NULL;
    gamma->endo->pibARMA = NULL;
    gamma->endo->pibARMA_ck = NULL;
    gamma->endo->optX = NULL;
    gamma->endo->optXChange = NULL;
    gamma->endo->optY = NULL;
    gamma->endo->optYChange = NULL;
    
    // (Sparse Vector) Allocate memory to bARMA
    if ( !(gamma->endo->bARMA = (sparseVector *) mem_malloc(sizeof(sparseVector))) ) {
        errMsg("allocation","newGamma","tsbBar",0);
        return NULL;
    } else {
        gamma->endo->bARMA->cnt = 0;
        if (!(gamma->endo->bARMA->val = (vector) arr_alloc(subbBarLength+1, double))) {
            errMsg("allocation","newGamma","tsbBar->val",0);
            return NULL;
        }
        if (!(gamma->endo->bARMA->col = (intvec) arr_alloc(subbBarLength+1, int))) {
            errMsg("allocation","newGamma","tsbBar->col",0);
            return NULL;
        }
    }
    
    // Allocate memory to pibma
    if ( !(gamma->endo->pibARMA = (vector) arr_alloc(lambdaCnt, double)) ) {
        errMsg("allocation","newGamma","pibar",0);
        return NULL;
    }
    // Allocate memory to ck
    if ( !(gamma->endo->pibARMA_ck = (intvec) arr_alloc(lambdaCnt, int)) ) {
        errMsg("allocation","newGamma","pibar",0);
        return NULL;
    }
    
    /* (Sparse Vector) Allocate Memory to Offset Values Space for late warm-up update */
    if(!(gamma->endo->changebARMA = (sparseVector *) mem_malloc(sizeof(sparseVector)))){
        errMsg("allocation","rolling_update","offsetb",0);
    } else {
        gamma->endo->changebARMA->cnt = 0;
        if ( !(gamma->endo->changebARMA->col = (intvec) arr_alloc(subbBarLength+1, int)) ) {
            errMsg("allocation","rollingUpdate","offsetb->col",0);
            return NULL;
        }
        if ( !(gamma->endo->changebARMA->val = (vector) arr_alloc(subbBarLength+1, double)) ) {
            errMsg("allocation","rollingUpdate","offsetb->val",0);
            return NULL;
        }
    }
    
    /* (Vector) Allocate memory to store best incumbent solution from the just solved time period */
    if (!(gamma->endo->optX = (vector) arr_alloc(masterCols + 1, double))) {
        errMsg("allocation", "newGamma", "endo->optX", 0);
        return NULL;
    }
    if (!(gamma->endo->optXChange = (vector) arr_alloc(masterCols + 1, double))) {
        errMsg("allocation", "newGamma", "endo->optXChange", 0);
        return NULL;
    }
    if (!(gamma->endo->optY = (vector) arr_alloc(subprobCols + 1, double))) {
        errMsg("allocation", "newGamma", "endo->optY", 0);
        return NULL;
    }
    if (!(gamma->endo->optYChange = (vector) arr_alloc(subprobCols + 1, double))) {
        errMsg("allocation", "newGamma", "endo->optYChange", 0);
        return NULL;
    }
    
    return gamma;
}

void freeExoType(exoType *exo){
    
    int i, n, stageCnt=2;
    
    // Releasing memorty that hold external data
    if (exo->aData) {
        if (run.DETERMINISTIC == 1)
            stageCnt = 1;
        for (n=0; n<stageCnt; n++) {
            if (exo->aData[n]) {
                for (i=1; i<=run.HORIZON; i++)
                    if (exo->aData[n][i])
                        mem_free(exo->aData[n][i]);
                mem_free(exo->aData[n]);
            }
        }
        mem_free(exo->aData);
    }
    if (exo->aDelta) mem_free(exo->aDelta);
    if (exo->BDelta) mem_free(exo->BDelta);
    if (exo->aDynLambda) {
        for (i=0; i<exo->aDynLambda->cnt; i++)
            if (exo->aDynLambda->vals[i]) mem_free(exo->aDynLambda->vals[i]);
        if (exo->aDynLambda->vals) mem_free(exo->aDynLambda->vals);
        mem_free(exo->aDynLambda);
    }
    if (exo->BDynLambda) {
        for (i=0; i<exo->BDynLambda->cnt; i++)
            if (exo->BDynLambda->vals[i]) mem_free(exo->BDynLambda->vals[i]);
        if (exo->BDynLambda->vals) mem_free(exo->BDynLambda->vals);
        mem_free(exo->BDynLambda);
    }
    
    mem_free(exo);
}//End of freeExoType

void freeEndoType(endoType *endo){
    if (endo->bARMA) freeSparseVector(endo->bARMA);
    if (endo->changebARMA) freeSparseVector(endo->changebARMA);
    if (endo->pibARMA_ck) mem_free(endo->pibARMA_ck);
    if (endo->pibARMA) mem_free(endo->pibARMA);
    if (endo->optX) mem_free(endo->optX);
    if (endo->optXChange) mem_free(endo->optXChange);
    if (endo->optY) mem_free(endo->optY);
    if (endo->optYChange) mem_free(endo->optYChange);
    mem_free(endo);
}//End of freeEndoType

void freeGammaType(gammaType *gamma) {
    if (gamma->endo) freeEndoType(gamma->endo);
    if (gamma->exo) freeExoType(gamma->exo);
    mem_free(gamma);
}//End of freeGammaType
