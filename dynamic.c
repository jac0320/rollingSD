//
//  dynamic.c
//  rollingSP
//
//

#include "dynamic.h"

extern runType run;

int readDyn(string inputDir, string probName, cellType *cell, probType **prob, int stageCnt) {
    
    char    probpath[2*BLOCKSIZE], line[BLOCKSIZE], **fields, fieldType;
    FILE    *fptr;
    int     status, i, n, numFields, masterStage = 0, subprobStage = 0, aSeg = 0;
    BOOL    initFlag = FALSE;
    
    if (run.DETERMINISTIC == 0) {
        masterStage = 0;
        subprobStage = 1;
    }
    
    /* -------------------- Initialization of Pointers in probType -------------------- */
    for (n=0; n<stageCnt; n++) {
        if (!(prob[n]->Bbar = (sparseMatrix **) arr_alloc(stageCnt, sparseMatrix *))) {
            errMsg("allocation", "readDyn", "prob->Bbar", 0);
            return 1;
        }
        // Initial NULL pointer
        for (i=0; i<stageCnt; i++) {
            prob[n]->Bbar[i] = NULL;
        }
        prob[n]->aBar = NULL;
    }
    /* --------------------------------------------------------------------------------- */
    
    sprintf(probpath, "%s%s/%s.dyn", inputDir, probName, probName);
    
    /* Open the dynamic file */
    fptr = fopen(probpath, "r");
    if (fptr == NULL) {
        errMsg("read", "readDyn", "failed to open dynamic file :: .dyn missing", 0);
        return 1;
    }
    
    /* allocate memory to seven field locations */
    if (!(fields = (string *) arr_alloc(7, string)) )
        errMsg("allocation", "readDyn", "field locations", 0);
    for (n = 0; n < 7; n++ )
        if ( !(fields[n] = (string) arr_alloc(NAMESIZE, char)) )
            errMsg("allocation", "readDyn", "individual field location", 0);
    
    /* Read and check problem name  */
    if (fgets(line, sizeof line, fptr) != NULL)
        sscanf(line, "%s %s", fields[0], fields[1]);
    else {
        errMsg("read", "readDyn", "failed to read the problem name", 0);
        return 1;
    }
    if ( !(strncmp(fields[0], "NAME", 4)))
        if ( strcmp(probName, fields[1]) ) {
            errMsg("read", "readDyn", "Problem name doesn't match", 0);
            return 1;
        }
    
    /* Read and check problem name  */
    if (fgets(line, sizeof line, fptr) != NULL)
        sscanf(line, "%s %s", fields[0], fields[1]);
    else {
        errMsg("read", "readDyn", "failed to read the problem name", 0);
        return 1;
    }
    if ( !(strncmp(fields[0], "TYPE", 4)))
        if ( strcmp("DYNAMIC", fields[1]) ) {
            errMsg("read", "readDyn", "Problem name doesn't match", 0);
            return 1;
        }
    
    while ( !(getLine(&fptr, fields, &fieldType, &numFields)) ) {
    RESUME_READ:
        if ( !(strcmp(fields[0], "BM")) ) {
            status = readBM(fptr, fields, prob);
            if (status)
                return 1;
            goto RESUME_READ;
        } else if ( !(strcmp(fields[0], "AV")) ) {
            status = readAV(fptr, fields, prob, &aSeg);
            if (status)
                return 1;
            goto RESUME_READ;
        } else if ( !(strcmp(fields[0], "INIT")) ) {
            status = readDynInit(fptr, inputDir, fields, probName, &initFlag, prob, cell, stageCnt, aSeg);
            if (status)
                return 1;
            initFlag = TRUE;
            goto RESUME_READ;
        } else if ( !(strcmp(fields[0], "ENDATA"))) {
            break;
        } else
            continue;
    }
    
    
    for (n=0; n<7; n++)
        if (fields[n]) mem_free(fields[n]);
    mem_free(fields);
    
    fclose(fptr);
    
    return 0; //Suppose to return dynType
}//End of _readDyn

int readBM(FILE *fptr, string *fields, probType **prob) {
    
    char    strType;
    int     n, i,BMCnt = 0;
    int     numFields, expOffset = 0, expSize = 0, rowBufCnt, colBufCnt, coefBufCnt, targetPeriod = 0, dynPeriod = 0;
    intvec  rowBuf, colBuf, newCol, newRow; //Buffer space
    vector  coefBuf, newVal; //Buffer space
    
    /*  This is what BM portion looks like
    |BM    PERIOD00
    | ...COL    ...ROW      ...COEF
    |BM    PERIOD10
    | ...COL    ...ROW      ...COEF
    |BM    PERIOD11
    | ...COL    ...ROW      ...COEF
     */
    
    if (run.DETERMINISTIC == 0) {
        dynPeriod = (int) (fields[1][0] - '0');
        targetPeriod = (int) (fields[1][1] - '0');
    } else {
        // Manual override when solving thing deterministically
        dynPeriod = 0;
        targetPeriod = 0;
    }
    
    rowBufCnt = max(prob[dynPeriod]->sp->mar,prob[targetPeriod]->sp->mar);
    colBufCnt = max(prob[dynPeriod]->sp->mac,prob[targetPeriod]->sp->mac);
    coefBufCnt = max(rowBufCnt, colBufCnt);
    if (!(rowBuf = (intvec) arr_alloc(rowBufCnt, int)))
        errMsg("allocation", "_readBM", "buffer for rows", 0);
    if (!(colBuf = (intvec) arr_alloc(colBufCnt, int)))
        errMsg("allocation", "_readBM", "buffer for rows", 0);
    if (!(coefBuf = (vector) arr_alloc(coefBufCnt, double)))
        errMsg("allocation", "_readBM", "buffer for rows", 0);
    
    if ( prob[dynPeriod]->Bbar[targetPeriod] == NULL ){
        if ( !(prob[dynPeriod]->Bbar[targetPeriod] = (sparseMatrix *) mem_malloc(sizeof(sparseMatrix))) )
            errMsg("allocation", "_readBM", "give a sparseMatrix structure", 0);
    } else {
        //Need to expand sparseMatrix in this case
        expOffset = prob[dynPeriod]->Bbar[targetPeriod]->cnt;
    }
    
    while (TRUE) {
        getLine(&fptr, fields, &strType, &numFields);
        if ( strType == 't' )
            break;
        if ( !strncmp("C",fields[0], 1) ) {
            BMCnt++;
            // Search for coordination of this column in dynamic oneProblem: column source
            n = 0;
            while ( n < prob[dynPeriod]->sp->mac ) {
                if (!(strncmp(fields[0],prob[dynPeriod]->sp->cname[n],NAMESIZE))) {
                    break;
                } else {
                    n++;
                }
            }
            colBuf[BMCnt] = n + 1; //Convert to algorithm index mode
            if ( !strncmp("R", fields[1], 1) ) {
                // Search for coordination of this rows in target oneProblem: row source
                n = 0;
                while ( n < prob[targetPeriod]->sp->mar ) {
                    if (!(strncmp(fields[1], prob[targetPeriod]->sp->rname[n], NAMESIZE))){
                        break;
                    } else {
                        n++;
                    }
                }
                rowBuf[BMCnt] = n + 1; //Convert to algorithm index mode
                coefBuf[BMCnt] = str2float(fields[2]);
            }
        } else {
            errMsg("reading", "_readBM", "false column indication", 0);
            goto ERROR_BM;
        }
    }
    
    //Just finished reading a block of matrix
    if (expOffset > 0) {
        // Case when attaching
        expSize = expOffset + BMCnt;
        // Allocate memory to new vectors
        if ( !(newVal = (vector) arr_alloc(expSize+1, double)) )
            errMsg("allocation", "_readBM", "Bbar->val", 0);
        if ( !(newCol = (intvec) arr_alloc(expSize+1, int)) )
            errMsg("allocation", "_readBM", "Bbar->col", 0);
        if ( !(newRow = (intvec) arr_alloc(expSize+1, int)) )
            errMsg("allocation", "_readBM", "Bbar->row", 0);
        // Copy the existing information
        for (i=1; i<=expOffset; i++) {
            newVal[i] = prob[dynPeriod]->Bbar[targetPeriod]->val[i];
            newCol[i] = prob[dynPeriod]->Bbar[targetPeriod]->col[i];
            newRow[i] = prob[dynPeriod]->Bbar[targetPeriod]->row[i];
        }
        // Free the old vectors
        mem_free(prob[dynPeriod]->Bbar[targetPeriod]->val);
        mem_free(prob[dynPeriod]->Bbar[targetPeriod]->col);
        mem_free(prob[dynPeriod]->Bbar[targetPeriod]->row);
        // Attach the new vectors
        prob[dynPeriod]->Bbar[targetPeriod]->val = newVal;
        prob[dynPeriod]->Bbar[targetPeriod]->row = newRow;
        prob[dynPeriod]->Bbar[targetPeriod]->col = newCol;
    } else {
        // Case when not attaching
        if ( !(prob[dynPeriod]->Bbar[targetPeriod]->val = (vector) arr_alloc(BMCnt+1, double)) )
            errMsg("allocation", "_readBM", "Bbar->val", 0);
        if ( !(prob[dynPeriod]->Bbar[targetPeriod]->col = (intvec) arr_alloc(BMCnt+1, int)) )
            errMsg("allocation", "_readBM", "Bbar->col", 0);
        if ( !(prob[dynPeriod]->Bbar[targetPeriod]->row = (intvec) arr_alloc(BMCnt+1, int)) )
            errMsg("allocation", "_readBM", "Bbar->row", 0);
    }
    
    prob[dynPeriod]->Bbar[targetPeriod]->cnt = expOffset + BMCnt;
    for (i=1; i<=BMCnt; i++) {
        prob[dynPeriod]->Bbar[targetPeriod]->col[expOffset+i] = colBuf[i];
        prob[dynPeriod]->Bbar[targetPeriod]->row[expOffset+i] = rowBuf[i];
        prob[dynPeriod]->Bbar[targetPeriod]->val[expOffset+i] = coefBuf[i];
    }
    
    mem_free(colBuf);
    mem_free(rowBuf);
    mem_free(coefBuf);
    return 0;
    
ERROR_BM:
    mem_free(colBuf);
    mem_free(rowBuf);
    mem_free(coefBuf);
    return 1;
}

int readAV(FILE *fptr, string *fields, probType **prob, int *aSeg) {
    
    char    strType;
    int     numFields, n, i, BMCnt = 0;
    int     colBufCnt, coefBufCnt, expOffset = 0, expSize = 0, dynPeriod = 0;
    intvec  colBuf, newCol; //Buffer space
    vector  coefBuf, newVal; //Buffer space
    
    /*  This is what BM portion looks like
     |AV    PERIOD00
     | ...ROW      ...COEF
     |AV    PERIOD01
     | ...ROW      ...COEF
     */
    
    if (run.DETERMINISTIC == 0) {
        dynPeriod = (int) (fields[1][0] - '0');
    } else {
        // Manual override when solving thing deterministically
        dynPeriod = 0;
    }
    
    colBufCnt = max(prob[dynPeriod]->sp->mac,prob[dynPeriod]->sp->mar);
    coefBufCnt = colBufCnt;
    if (!(colBuf = (intvec) arr_alloc(colBufCnt, int)))
        errMsg("allocation", "_readBM", "buffer for rows", 0);
    if (!(coefBuf = (vector) arr_alloc(coefBufCnt, double)))
        errMsg("allocation", "_readBM", "buffer for rows", 0);

    if ( prob[dynPeriod]->aBar == NULL) {
        if ( !(prob[dynPeriod]->aBar = (sparseVector *) mem_malloc(sizeof(sparseVector))) )
            errMsg("allocation", "_readBM", "give a sparseMatrix structure", 0);
    }else{
        //Need to expand sparseMatrix in this case
        expOffset = prob[dynPeriod]->aBar->cnt;
    }

    
    while (TRUE) {
        getLine(&fptr, fields, &strType, &numFields);
        if ( strType == 't')
            break;
        if ( !strncmp("R",fields[0], 1) ) {
            BMCnt++;
            // Search for coordination of this column in dynamic oneProblem: column source
            n = 0;
            while ( n < prob[dynPeriod]->sp->mar ) {
                if (!(strncmp(fields[0],prob[dynPeriod]->sp->rname[n],NAMESIZE))) {
                    break;
                } else {
                    n++;
                }
            }
            colBuf[BMCnt] = n + 1; //Convert to algorithm index mode
            coefBuf[BMCnt] = str2float(fields[1]);
        } else {
            errMsg("reading", "_readBM", "false column indication", 0);
            goto ERROR_AV;
        }
        
    }
    
    //Just finished reading a block of matrix
    if (expOffset > 0) {
        // Case when attaching
        expSize = expOffset + BMCnt;
        *aSeg = expOffset;
        // Allocate memory to new vectors
        if ( !(newVal = (vector) arr_alloc(expSize+1, double)) )
            errMsg("allocation", "_readBM", "Bbar->val", 0);
        if ( !(newCol = (intvec) arr_alloc(expSize+1, int)) )
            errMsg("allocation", "_readBM", "Bbar->col", 0);
        // Copy the existing information
        for (i=1; i<=expOffset; i++) {
            newVal[i] = prob[dynPeriod]->aBar->val[i];
            newCol[i] = prob[dynPeriod]->aBar->col[i];
        }
        // Free the old vectors
        mem_free(prob[dynPeriod]->aBar->val);
        mem_free(prob[dynPeriod]->aBar->col);
        // Attach the new vectors
        prob[dynPeriod]->aBar->val = newVal;
        prob[dynPeriod]->aBar->col = newCol;
    } else {
        if ( !(prob[dynPeriod]->aBar->val = (vector) arr_alloc(BMCnt+1, double)) )
            errMsg("allocation", "_readBM", "Bbar->val", 0);
        if ( !(prob[dynPeriod]->aBar->col = (intvec) arr_alloc(BMCnt+1, int)) )
            errMsg("allocation", "_readBM", "Bbar->col", 0);
    }

    prob[dynPeriod]->aBar->cnt = expOffset + BMCnt;
    for (i=1; i<=BMCnt; i++) {
        prob[dynPeriod]->aBar->col[expOffset+i] = colBuf[i];
        prob[dynPeriod]->aBar->val[expOffset+i] = coefBuf[i];
    }
    
    mem_free(colBuf);
    mem_free(coefBuf);
    return 0;
    
ERROR_AV:
    mem_free(colBuf);
    mem_free(coefBuf);
    return 1;
}

int readDynInit(FILE *fptr, string inputDir, string *fields, string probName, BOOL *initFlag,
                 probType **prob, cellType *cell, int stageCnt, int aSeg) {
    
    char    strType, line[BLOCKSIZE], field1[NAMESIZE], initFilePath[BLOCKSIZE];
    int     numFields, n, i, j, dynPeriod, targetPeriod;
    double  param;
    string  *longFields;
    FILE    *initFptr = NULL;
    
    if (!(longFields = (string *) arr_alloc(7, string)) )
        errMsg("allocation", "readDyn", "field locations", 0);
    for (n = 0; n < 7; n++ )
        if ( !(longFields[n] = (string) arr_alloc(2*BLOCKSIZE, char)) )
            errMsg("allocation", "readDyn", "individual field location", 0);
    
    /* Memory Initialization on the cellType :: Move this later to cell */
    cell->gamma->exo->aData = (vector **) arr_alloc(stageCnt, vector *);
    for (i=0; i<stageCnt; i++) {
        cell->gamma->exo->aData[i] = (vector *) arr_alloc(run.HORIZON+1, vector);
        for (j=1; j<=run.HORIZON; j++)
            cell->gamma->exo->aData[i][j] = (vector) arr_alloc(prob[i]->aBar->cnt+1, double);
    }
    if (run.DETERMINISTIC == 0) {
        if (prob[1]->aBar) {
            cell->gamma->exo->aCnt = prob[1]->aBar->cnt;
            if (!(cell->gamma->exo->aDelta = (vector) arr_alloc(prob[1]->aBar->cnt+1, double)))
                errMsg("allocation", "readDynInit", "gamma->exo->aDelta", 0);
        }
        if (prob[1]->Bbar[1]) {
            cell->gamma->exo->BCnt = prob[1]->Bbar[1]->cnt;
            if (!(cell->gamma->exo->BDelta = (vector) arr_alloc(prob[1]->Bbar[1]->cnt+1, double)))
                errMsg("allocation", "readDynInit", "gamma->exo->BDelta", 0);
        }
    }


    /*  This is what INIT file portion looks like
     |NAME      ****
     |TYPE      DYNINIT
     |BM        0
     | ### ### ### ### ...
     |BM        1
     | ### ### ### ### ...
     |AV        0
     | T ### ### ### ###
     | T ### ### ###
     | ...
     |AV        1
     | T ### ### ### ###
     | T ### ### ### ###
     | ...
     |ENDATA
    */
    
    if ( !strcmp("ENDATA", fields[0]) )  { //Enter the ENDATA, no init path indicated
        initFlag = FALSE;
        return 0;
    }
    // Fetch the path line
    while (TRUE) {
        getLine(&fptr, longFields, &strType, &numFields);
        if ( !(strcmp("ENDATA", longFields[0])) )
            break; //Exit point of this while loop: when finish reading the init file
        
        /* Open the dynamic file */
        if (!(strcmp("PATH", longFields[0]))) {
            if ( !(strncmp(longFields[1], "./", 2))) {
                strcpy(initFilePath, inputDir);
                strcat(initFilePath, probName);
                strcat(initFilePath, "/");
                strcat(initFilePath, probName);
                strcat(initFilePath, ".init");
                initFptr = fopen(initFilePath, "r");
            } else {
                initFptr = fopen(longFields[1], "r");
            }
            if (initFptr == NULL) {
                errMsg("read", "readDyn", "failed to open dynamic file :: .init missing", 0);
                return 1;
            }
        } else {
            errMsg("read", "readDynInit", "wrong indication of init file path", 0);
            return 1;
        }

        /* Read and check problem name  */
        if (fgets(line, sizeof line, initFptr) != NULL)
            sscanf(line, "%s %s", longFields[0], longFields[1]);
        else {
            errMsg("read", "readDyn", "failed to read the problem name", 0);
            return 1;
        }
        if ( !(strncmp(longFields[0], "NAME", 4)))
            if ( strcmp(probName, longFields[1])) {
                errMsg("read", "readDyn", "Problem name doesn't match", 0);
                return 1;
            }
        
        /* Read and check the type name */
        if (fgets(line, sizeof line, initFptr) != NULL)
            sscanf(line, "%s %s", longFields[0], longFields[1]);
        else {
            errMsg("read", "_readDynInit", "failed to read file type", 0);
        }
        if ( !(strncmp(longFields[0], "TYPE", 4)))
            if ( strncmp(longFields[1], "DYNINIT", 4) ) {
                errMsg("read", "_readDynInit", "wrong type line", 0);
            }
        
        
        while (TRUE) {
            getLine(&initFptr, fields, &strType, &numFields);
        RESUME_INIT_READING:
            if (strType == 't') {
                if (!(strcmp(fields[0], "BM"))) {
                    if (run.DETERMINISTIC == 0) {
                        dynPeriod = (int) (fields[1][6] - '0');
                        targetPeriod = (int) (fields[1][7] - '0');
                    }else{
                        dynPeriod = 0;
                        targetPeriod = 0;
                    }
                    for (n=1; n<=prob[dynPeriod]->Bbar[targetPeriod]->cnt; n++) {
                        fscanf(initFptr, "%s", field1);
                        while (field1[0] < '0' || field1[0] > '9')
                            fscanf(initFptr, "%s", &field1[0]);
                        param = str2float(field1);
                        if (dynPeriod == 0)
                            cell->gamma->endo->optX[prob[dynPeriod]->Bbar[targetPeriod]->col[n]] = param;
                        if (dynPeriod == 1)
                            cell->gamma->endo->optY[prob[dynPeriod]->Bbar[targetPeriod]->col[n]] = param;
                    }
                    getLine(&initFptr, fields, &strType, &numFields);
                    goto RESUME_INIT_READING;
                } else if ( !(strcmp(fields[0], "AV")) ) {
                    if (run.DETERMINISTIC == 0) {
                        dynPeriod = (int) (fields[1][0] - '0');
                        for (n=1; n<=run.HORIZON; n++) {
                            fscanf(initFptr, "%s", field1); //This is the row/horizon index
                            for (j=1; j<=prob[dynPeriod]->aBar->cnt; j++) {
                                fscanf(initFptr, "%s", field1);
                                param = str2float(field1);
                                cell->gamma->exo->aData[dynPeriod][n][j] = param;
                            }
                        }
                    }else{
                        // Special treatment for deterministic case
                        dynPeriod = 0;
                        for (n=1; n<=run.HORIZON; n++) {
                            fscanf(initFptr, "%s", field1); //This is the row/horizon index
                            for (j=1; j<=aSeg; j++) {
                                fscanf(initFptr, "%s", field1);
                                param = str2float(field1);
                                cell->gamma->exo->aData[dynPeriod][n][j] = param;
                            }
                        }
                        fscanf(initFptr, "%s", &field1[0]); //Expected  encounter
                        while (field1[0] < '0' || field1[0] > '9')
                            fscanf(initFptr, "%s", &field1[0]);
                        fscanf(initFptr, "%s", &field1[0]); //Expected encounter of index
                        for (n=1; n<=run.HORIZON; n++) {
                            for (j=aSeg + 1; j<=prob[dynPeriod]->aBar->cnt; j++) {
                                fscanf(initFptr, "%s", field1);
                                param = str2float(field1);
                                cell->gamma->exo->aData[dynPeriod][n][j] = param;
                            }
                            if (n < run.HORIZON) //Skip the last ENDATA
                                fscanf(initFptr, "%s", field1); //This is the row/horizon index
                        }
                    }
                    getLine(&initFptr, fields, &strType, &numFields);
                    goto RESUME_INIT_READING;
                } else if (!(strcmp(fields[0], "ENDATA"))) {
                    goto FINISHED_INITFILE;
                } else {
                    errMsg("read", "_readDynInit", "Unkown input line in dynInit file", 0);
                    return 1;
                }
            }
        }
    }
FINISHED_INITFILE:
    for (n=0; n<7; n++)
        if (longFields[n]) mem_free(longFields[n]);
    mem_free(longFields);
    
    fclose(initFptr);
    
    return 0;
}
