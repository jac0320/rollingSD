
//
//  rollingARMA.c
//  rollingSP
//


#include "arma.h"

extern runType run;
extern string outputDir;

tsType *readARIMA(string inputDir, string probName, oneProblem *orig, BOOL *arimaFlag, int *stat1) {
    
    char        probpath[BLOCKSIZE], line[BLOCKSIZE], field1[NAMESIZE], field2[NAMESIZE];
    int         status, i;
    double      temp;
    FILE        *fptr = NULL;
    tsType      *arima = NULL;
    
    /* Locate the problem arma file */
    sprintf(probpath, "%s%s/%s.arma", inputDir, probName, probName);
    fptr = fopen(probpath, "r");
    if ( fptr == NULL ) {
        printf("No .arma file detected...\n");
        *arimaFlag = FALSE;
        *stat1 = 1; return NULL;
    } else {
        *arimaFlag = TRUE;
    }
    
    /*************************** NAME SECTION: Check time series file problem name *****************************/
    if ( fgets(line, sizeof line, fptr) != NULL ) {
        sscanf(line, "%s %s", field1, field2);
    }
    else {
        errMsg("read", "readARMA", "Failed to read name in ARMA file", 0);
        *stat1 = 1; return NULL;
    } // Ready to go
    if ( !(strncmp(field1, "NAME", 4)) ) {
        if( strcmp(probName, field2) ) {
            errMsg("read", "read_arma", "problem name does not match, in NAME section", 0);
            *stat1 = 1; return NULL;
        }
    }
    
    /* Allocate to initialize */
    if ( !(arima = (tsType *) mem_malloc(sizeof(tsType)) ) ) {
        errMsg("allocation", "read_arma", "Failed to allocate memory to tsType",0);
        *stat1 = 1;
        return NULL;
    }
    /************************************************************************************************************/
    
    /* Initialization of All Variates */
    arima->lb = 0.0;
    arima->var = 0;
    arima->subh = 0;
    arima->p = 0;
    arima->q = 0;
    arima->d = 0;
    arima->window = 0;
    arima->histARObs = NULL;
    arima->histMAObs = NULL;
    arima->decomp = NULL;
    arima->history = NULL;
    arima->mean = NULL;
    arima->stdev = NULL;
    arima->variate = NULL;
    arima->Theta = NULL;
    arima->Phi = NULL;
    
    /* Read the main coefficients of the time series model */
    status = readARIMACore(fptr, arima, stat1);
    if (status) {
        errMsg("read", "readARMA", "Failed to read core ARMA information", 0);
        return NULL;
    }
    
    /* Read Phi and Theta of the time series model */
    status = readARIMACoef(fptr, arima, stat1);
    if (status) {
        errMsg("read", "readARMA", "ARIMA coefficients", 0);
        return NULL;
    }
    
    /* Read each variate information */
    status = readARIMAVarBlock(fptr, arima, orig, stat1);
    if (status) {
        errMsg("read", "readARIMA", "ARIMA variable blocks", 0);
        return NULL;
    }
    
    /* Read historical time series observation and noise */
    status = readARIMAStream(fptr, arima, stat1);
    if (status) {
        errMsg("read", "readARIMA", "ARIMA stream data", 0);
        return NULL;
    }
    
    /* Read trend and seasonality or other decomposed deterministical information */
    status = readARIMAAux(fptr, arima, stat1);
    if (status) {
        errMsg("read", "readARIMA", "ARIMA aux data for reverse pre-process", 0);
        return NULL;
    }
    
    /************************ ENDATA Section: Detecting the end of the arma file ************************/
    if ( fgets(line,sizeof line, fptr) != NULL ) {
        sscanf(line, "%s", field1);
    } else {
        errMsg("read","read_arma","Failed to read Ending signal",0);
        *stat1 = 1; return NULL;
    } // Ready to go
    if ( strncmp(field1, "ENDATA", 6)) {
        errMsg("read", "read_arma", "Failed to detect the end of a file.", 0);
        *stat1 = 1; return NULL;
    }
    /****************************************************************************************************/
    
    // Close .arma file
    fclose(fptr);
    
#ifdef REPORT
    printf("Successfully Read Time Series Model for %s\n", probName);
    printf("Time Series Model Type -> %d\n", arima->transType);
    printf("Total Sub-hourly Intervals -> %d\n", arima->subh);
    printf("Model Order -> (p=%d, d=%d, q=%d)\n", arima->p, arima->d, arima->q);
    printf("Total Variate -> %d\n", arima->var);
    printf("Historical Observation Stream Length ->%d\n", arima->history->total);
    printf("Historical Observation Head Values -> \n \t ::");
    for (i=1; i<=5; i++) printf("[%d]%f ",arima->history->open+i, arima->history->obs[arima->history->open+i][1]);
    printf("\nHistorical Noise Stream -> \n \t ::");
    for (i=1; i<=5; i++) printf("[%d]%f ",arima->noise->open+i, arima->noise->obs[arima->noise->open+i][1]);
    printf("\nAux Data Stream (Trend) -> \n \t ::");
    for (i=1; i<=5; i++) printf("[%d]%f ",i,arima->decomp[1]->trend[i]);
    printf("\nAux Data Stream (Seasonal) -> \n \t ::");
    for (i=1; i<=5; i++) printf("[%d]%f ",i,arima->decomp[1]->seasonal[i]);
#endif
    
    /* Some very basic validation work is conducted here to see if we have a valid ARIMA model input. */
    printf("\n\n");
    printf("-----Validating ARIMA File Input\n");
    // Checking Phi Coefficients Amount
    printf("Matching Phi sparseMatrix Length with model...");
    temp = 0;
    for (i=1; i<=arima->p; i++)
        temp += arima->Phi[i]->cnt;
    if ( (arima->var * arima->var * arima->p) == temp) {
        printf("Passed\n");
    } else {
        printf("Failed\n");
        errMsg("validate","readARMA","Not enough Phi coefficients.",0);
        *stat1 = 1; return NULL;
    }
    // Checking Theta Coefficient Amount
    printf("Matching Theta sparseMatrix Length with model...");
    temp = 0;
    for (i=1; i<=arima->q; i++)
        temp += arima->Theta[i]->cnt;
    if ((arima->var*arima->var*arima->q) == temp) {
        printf("Passed\n");
    } else {
        printf("Failed\n");
        errMsg("validate","readARMA","Not enough theta coefficents.",0);
        *stat1 = 1; return NULL;
    }
    // Checking if historical data is sufficient for rolling horizon problem using run.HORIZON

    printf("-----Finish Validation\n\n");
    
    return arima;
}

int readARIMACore(FILE *fptr, tsType *arima, int *stat1) {
    
    char        line[BLOCKSIZE], field1[NAMESIZE], field2[NAMESIZE];
    int         int_param;
    
    /********************************* TYPE SECTION: Collect time series types **********************************/
    if (fgets(line, sizeof line, fptr) != NULL) {
        sscanf(line, "%s %s", field1, field2);
    }
    else {
        errMsg("read","read_arma","Failed to read time series type",0);
        *stat1 = 1;
        return 1;
    } // Ready to read in
    if ((!strncmp(field1, "TYPE", 4))) {
        if (!(strncmp(field2, "UNIVAR", 6)) && !(strncmp(field2,"MULTIVAR",8))) {
            errMsg("read", "read_arma", "Unknown time series type",0);
            *stat1 = 1;
            return 1;
        }
    }
    if ( !(strncmp(field2, "EMP",3)) ) {
        arima->transType = 1;
    } else if (  !(strncmp(field2, "STL",3)) ) {
        arima->transType = 2;
    } else if (  !(strncmp(field2, "CLA",3)) ) {
        arima->transType = 3;
    } else {
        errMsg("read","readArma","Unkown timer series inversion type",0);
        *stat1 = 1;
        return 1;
    }
    /************************************************************************************************************/
    
    /**************************** PARAMETER SECTION: Collect time series parameters *****************************/
    //------Parameter SUBHOUR
    if (fgets(line, sizeof line, fptr) != NULL) {
        sscanf(line, "%s %d", field1, &int_param);
    }
    else {
        errMsg("read", "read_arma","Failed to read parameter P",0);
        *stat1 = 1;
        return 1;
    } // Ready to read in
    
    if ( strncmp(field1, "SUBHOUR", 7) ) {
        errMsg("read", "read_arma", "Order P not detected",0);
        *stat1 = 1;
        return 1;
    }
    arima->subh = int_param;
    
    //------Parameter P
    if (fgets(line, sizeof line, fptr) != NULL) {
        sscanf(line, "%s %d", field1, &int_param);
    }
    else {
        errMsg("read", "read_arma","Failed to read parameter P",0);
        *stat1 = 1;
        return 1;
    } // Ready to read in
    
    if ( strncmp(field1, "P", 1) ) {
        errMsg("read", "read_arma", "Order P not detected",0);
        *stat1 = 1;
        return 1;
    }
    arima->p = int_param;
    
    //------Parameter D
    if (fgets(line, sizeof line, fptr) != NULL) {
        sscanf(line, "%s %d", field1, &int_param);
    }
    else {
        errMsg("read", "read_arma","Failed to read parameter D",0);
        *stat1 = 1;
        return 1;
    } // Ready to read in
    if ( strncmp(field1, "D", 1) ) {
        errMsg("read", "read_arma", "Difference Parameter",0);
        *stat1 = 1;
        return 1;
    }
    arima->d = int_param;
    
    //------Parameter Q
    if (fgets(line, sizeof line, fptr) != NULL) {
        sscanf(line, "%s %d", field1, &int_param);
    }
    else {
        errMsg("read", "read_arma","Failed to read parameter Q",0);
        *stat1 = 1;
        return 1;
    } // Ready to read in
    if ( strncmp(field1, "Q", 1) ) {
        errMsg("read", "read_arma", "Order Q not detected",0);
        *stat1 = 1;
        return 1;
    }
    arima->q = int_param;
    
    //------Parameter T
    // UNIVAR -> several time series with several vectors of coefficients in multiple lines
    // MULTIVAR -> one time series with a sparse matrix of coefficients
    // This really don't matter any more...
    if (fgets(line, sizeof line, fptr) != NULL) {
        sscanf(line, "%s %d", field1, &int_param);
    }
    else {
        errMsg("read", "read_arma","Failed to read parameter T",0);
        *stat1 = 1;
        return 1;
    } // Read to read in
    if ( strncmp(field1, "VAR", 3) ) {
        errMsg("read", "read_arma", "Instrument amount T not detected",0);
        *stat1 = 1;
        return 1;
    }
    arima->var = int_param;
    /**************************************************************************************************************/
    
    return 0;
}

int readARIMACoef(FILE *fptr, tsType *arima, int *stat1) {
    char        line[BLOCKSIZE], field1[NAMESIZE];
    int         int_param, int_param2;
    double      param;
    int         i, k;
    
    /************************** COEFFICIENT SECTION: Collect time series coefficients *****************************/
    if ( !(arima->Phi = (sparseMatrix **) arr_alloc(arima->p+1, sparseMatrix *))){
        errMsg("allocation","read_arma","Failed to allocate memory to Phi coefficients space.",0);
        *stat1 = 1;
        return 1;
    }
    for ( k=1; k<=arima->p; k++) {
        //----->Phi Section
        if (fgets(line, sizeof line, fptr) != NULL) {
            sscanf(line, "%s %d %d", field1, &int_param2, &int_param);
        }
        else {
            errMsg("read","read_arma","Failed to read AR coefficient section",0);
            *stat1 = 1;
            return 1;
        }// Ready to read in
        if ( strncmp(field1, "AR", 2) ) {
            errMsg("read", "read_arma", "Unknown time series type",0);
            *stat1 = 1;
            return 1;
        } // Ready to go
        //---->Prepare Memory Space for Coefficients SparseMatrix
        if (!(arima->Phi[k] = (sparseMatrix *) mem_malloc(sizeof(sparseMatrix)))) {
            errMsg("allocation", "read_arma","Failed allocate memory to Phi sparseMatrix",0);
            *stat1 = 1;
            return 1;
        }
        arima->Phi[k]->cnt = int_param;
        // Columns goes first
        if (!(arima->Phi[k]->col = (intvec) arr_alloc(arima->Phi[k]->cnt+1,int)))
            errMsg("allocation", "read_arma", "Failed to allocate memory to column space of Phi",0);
        if (!(arima->Phi[k]->row = (intvec) arr_alloc(arima->Phi[k]->cnt+1,int)))
            errMsg("allocation", "read_arma", "Failed to allocate memory to row space of Phi",0);
        if (!(arima->Phi[k]->val = (vector) arr_alloc(arima->Phi[k]->cnt+1,double)))
            errMsg("allocation","read_arma","Fialed to allocate memory to value space of Phi",0);
        for ( i=1; i<=arima->Phi[k]->cnt; i++) {
            fscanf(fptr, "%d", &int_param); arima->Phi[k]->col[i] = int_param;
        }
        // Detect Ending Mark, there is a double check here.
        if ( !(fscanf(fptr, "%s", field1)) ) {
            errMsg("read","read_arma","Failed to detect sparseMatrix end.",0);
            *stat1 = 1;
            return 1;
        } else {
            if (strncmp(field1, "E", 1)) {
                errMsg("read","read_arma","Wrong sparseMatrix end indicator.",0);
                *stat1 = 1;
                return 1;
            }
        }
        for ( i=1; i<=arima->Phi[k]->cnt; i++) {
            fscanf(fptr, "%d", &int_param); arima->Phi[k]->row[i] = int_param;
        }
        // Detect Ending Mark, there is a double check here.
        if ( !(fscanf(fptr, "%s", field1)) ) {
            errMsg("read","read_arma","Failed to detect sparseMatrix end.",0);
            *stat1 = 1;
            return 1;
        } else {
            if (strncmp(field1, "E", 1)) {
                errMsg("read","read_arma","Wrong sparseMatrix end indicator.",0);
                *stat1 = 1;
                return 1;
            }
        }
        for ( i=1; i<=arima->Phi[k]->cnt; i++) {
            fscanf(fptr, "%lf", &param); arima->Phi[k]->val[i] = param;
        }
        // Detect Ending Mark, there is a double check here.
        if ( !(fscanf(fptr, "%s", field1)) ) {
            errMsg("read","read_arma","Failed to detect sparseMatrix end.",0);
            *stat1 = 1;
            return 1;
        } else {
            if (strncmp(field1, "E", 1)) {
                errMsg("read","read_arma","Wrong sparseMatrix end indicator.",0);
                *stat1 = 1;
                return 1;
            }
        }
        // Skip this lines
        fgets(line, sizeof line, fptr);
    }
    
    //----->Theta Section
    if (!( arima->Theta = (sparseMatrix **) arr_alloc(arima->q+1, sparseMatrix *))){
        errMsg("allocation","read_arma","Failed to allocate memory to Phi coefficients space.",0);
        return 1;
    }
    for ( k=1; k<=arima->q; k++) {
        if (fgets(line, sizeof line, fptr) != NULL) {
            sscanf(line, "%s %d %d", field1, &int_param2, &int_param);
        } else {
            errMsg("read","read_arma","Failed to read MA coefficient section",0);
            *stat1 = 1;
            return 1;
        } // Ready to read in
        if ( strncmp(field1, "MA", 2) ) {
            errMsg("read", "read_arma", "Unknown time series type",0);
            *stat1 = 1;
            return 1;
        } // Ready to go
        if (!(arima->Theta[k] = (sparseMatrix *) mem_malloc(sizeof(sparseMatrix)))) {
            errMsg("allocation","read_arma","Failed to allocate memtory to Theta sparseMatrix",0);
            *stat1 = 1;
            return 1;
        }
        arima->Theta[k]->cnt = int_param;
        if (!(arima->Theta[k]->col = (intvec) arr_alloc(arima->Theta[k]->cnt+1,int)))
            errMsg("allocation", "read_arma", "Failed to allocate memory to column space of Phi",0);
        if (!(arima->Theta[k]->row = (intvec) arr_alloc(arima->Theta[k]->cnt+1,int)))
            errMsg("allocation", "read_arma", "Failed to allocate memory to row space of Phi",0);
        if (!(arima->Theta[k]->val = (vector) arr_alloc(arima->Theta[k]->cnt+1,double)))
            errMsg("allocation","read_arma","Fialed to allocate memory to value space of Phi",0);
        for ( i=1; i<=arima->Theta[k]->cnt; i++) {
            fscanf(fptr, "%d", &int_param); arima->Theta[k]->col[i] = int_param;
        }
        // Detect Ending Mark, there is a double check here.
        if ( !(fscanf(fptr, "%s", field1)) ) {
            errMsg("read","read_arma","Failed to detect sparseMatrix end.",0);
            *stat1 = 1;
            return 1;
        } else {
            if (strncmp(field1, "E", 1)) {
                errMsg("read","read_arma","Wrong sparseMatrix end indicator.",0);
                *stat1 = 1;
                return 1;
            }
        }
        for ( i=1; i<=arima->Theta[k]->cnt; i++) {
            fscanf(fptr, "%d", &int_param); arima->Theta[k]->row[i] = int_param;
        }
        // Detect Ending Mark, there is a double check here.
        if ( !(fscanf(fptr, "%s", field1)) ) {
            errMsg("read","read_arma","Failed to detect sparseMatrix end.",0);
            *stat1 = 1;
            return 1;
        } else {
            if (strncmp(field1, "E", 1)) {
                errMsg("read","read_arma","Wrong sparseMatrix end indicator.",0);
                *stat1 = 1;
                return 1;
            }
        }
        for ( i=1; i<=arima->Theta[k]->cnt; i++) {
            fscanf(fptr, "%lf", &param); arima->Theta[k]->val[i] = param;
        }
        // Detect Ending Mark, there is a double check here.
        if ( !(fscanf(fptr, "%s", field1)) ) {
            errMsg("read","read_arma","Failed to detect stream end.",0);
            *stat1 = 1;
            return 1;
        } else {
            if (strncmp(field1, "E", 1)) {
                errMsg("read","read_arma","Wrong stream end indicator.",0);
                *stat1 = 1;
                return 1;
            }
        }
        // Skip this lines
        fgets(line, sizeof line, fptr);
    }
    /**************************************************************************************************************/
    return 0;
}

int readARIMAVarBlock(FILE *fptr, tsType *arima, oneProblem *orig, int *stat1) {
    
    char        line[BLOCKSIZE], field1[NAMESIZE], field2[NAMESIZE];
    int         int_param;
    double      param, param2;
    int         i, j, n;
    
    /*********************** VAR BLOCK SECTION: Find out time series associated row names *************************/
    //------>Allocate memory to time series blocks
    if ( !(arima->variate = (varBlkType **) arr_alloc(arima->var+1, varBlkType *))) {
        errMsg("allocation","read_arma","Failed to allocate memeroty to variate blocks",0);
        *stat1 = 1;
        return 1;
    }
    if ( !(arima->mean = (vector) arr_alloc(arima->var, double)) ) {
        errMsg("allocation","readArma","Failed to allocate memory to arima->mean",0);
        *stat1 = 1;
        return 1;
    }
    if ( !(arima->stdev = (vector) arr_alloc(arima->var, double)) ) {
        errMsg("allocation","readArma","Failed to allocate memory to arima->stdev",0);
        *stat1 = 1;
        return 1;
    }
    //------>Get Row/Column Information for every variate
    for ( i = 1; i <=arima->var; i++) {
        // Read the head of a block
        if (fgets(line, sizeof line, fptr) != NULL) {
            sscanf(line, "%s %s %d %lf %lf", field1, field2, &int_param, &param, &param2);
        } else {
            errMsg("read","read_arma","failed to read the next line",0);
            *stat1 = 1;
            return 1;
        } // Ready to check content
        if ( strncmp(field1,"BL",2) ) {
            errMsg("read","read_arma","failed to find a block in right position",0);
            *stat1 = 1;
            return 1;
        }
        // Allocate memory when encounter block
        if (!( arima->variate[i] = (varBlkType *) mem_malloc(sizeof(varBlkType)))) {
            errMsg("allocation","read_arma","failed to allocate memory to variate block",0);
            *stat1 = 1;
            return 1;
        }
        // Record block names, record block rows count
        arima->mean[i-1] = param;
        arima->stdev[i-1] = param2;
        arima->variate[i]->cnt = int_param;
        // Prepare to record names fo random variables, first allocate memory to integer index vectors
        if ( !( arima->variate[i]->col = (intvec) arr_alloc(arima->variate[i]->cnt+1,int) ) ) {
            errMsg("allocation","read_arma","Failed to allocate memory to variate block column index int vector",0);
            *stat1 = 1;
            return 1;
        }
        if ( !( arima->variate[i]->row = (intvec) arr_alloc(arima->variate[i]->cnt+1,int) ) ) {
            errMsg("allocation","read_arma","Failed to allocate memory to variate block column index int vector",0);
            *stat1 = 1;
            return 1;
        }
        if ( !( arima->variate[i]->seq = (intvec) arr_alloc(arima->variate[i]->cnt+1,int) ) ) {
            errMsg("allocation","read_arma","Failed to allocate memory to variate block column index int vector",0);
            *stat1 = 1; return 1;
        }
        for ( j = 1; j<=arima->variate[i]->cnt; j++) {
            // Keep reading lines
            if (fgets(line, sizeof line, fptr) != NULL) {
                // field1->column; field->row
                sscanf(line, "%s  %s  %d", field1, field2, &int_param);
            } else {
                errMsg("read","read_arma","Failed to read variate block line",0);
                *stat1 = 1; return 1;
            } // Ready to check content
            //------>Record the sequence position
            arima->variate[i]->seq[j] = int_param;
            //------>Search for column coordinates
            if ( !(strncmp(field1, "RHS", 3)) ) {
                n = -1;
            } else {
                n = 1;
                while ( n < orig->mac ) {
                    if ( !(strncmp(field1, orig->cname[n],NAMESIZE)) ) {
                        break;
                    } else {
                        n++; //This is the results of the CPLEX problem row
                    }
                }
            }
            // Record Column Coordinates: The future column set starts from one in algorithm
            arima->variate[i]->col[j] = n+1;
            //------>Search for row coordinates
            if ( !(strncmp(field2, orig->objname, NAMESIZE))) {
                n = -1;
            } else {
                n = 1;
                while ( n < orig->mar ) {
                    if (!(strncmp(field2, orig->rname[n], NAMESIZE))) {
                        break;
                    } else {
                        n++; //This is the results of the CPLEX problem row
                    }
                }
            }
            // Record Column Coordinates: The future rhs starts from one in algorithm
            arima->variate[i]->row[j] = n+1;
        }
    }
    /**************************************************************************************************************/
    
    return 0;
}

/* Description goes here... */
int readARIMAStream(FILE *fptr, tsType *arima, int *stat1) {
    
    
    char        line[BLOCKSIZE], field1[NAMESIZE];
    int         int_param, int_param2;
    double      param;
    int         i, j;
    
    /********** STREAM SECTION: This section stores all the historical data from input time series model **********/
    //------>Allocate memory to time series blocks
    if ( !(arima->history = (streamType *) mem_malloc(sizeof(streamType)))) {
        errMsg("allocation","readARIMAStream","Failed to allocate memeroty to variate blocks",0);
        *stat1 = 1; return 1;
    }
    //------>Collect Data Stream for different variate
    if (fgets(line, sizeof line, fptr) != NULL) {
        sscanf(line, "%s %d %d", field1, &int_param, &int_param2);
    } else {
        errMsg("read","readARIMAStream","Failed to read the next line",0);
        *stat1 = 1; return 1;
    } // Ready to check
    if (strncmp(field1, "STREAM", 6)) {
        errMsg("read","readARIMAStream","Failed to find a bloc in right position", 0);
        *stat1 = 1; return 1;
    }
    printf("%d Historical Data Points Available. Starting index = %d \n", int_param, int_param2);
    arima->history->start = int_param2;
    arima->history->past  = int_param2-1;
    arima->history->diff  = 0;
    arima->history->open  = int_param2;
    arima->history->total = int_param;
    //------>Allocate memory to data stream structure
    if ( !(arima->history->obs = (vector *) arr_alloc(arima->history->total+1, vector) )) {
        errMsg("allocation","readARIMAStream","failed to allocate memory to data stream.",0);
        *stat1 = 1; return 1;
    }
    //------>Allocate memory to store historical data at each time point for all variates
    for (i=1; i<=arima->history->total; i++) {
        if (!((arima->history->obs[i] = arr_alloc(arima->var+1, double)))){
            errMsg("allocation","readARIMAStream","Failed to allocate memory to obs space",0);
            *stat1 = 1; return 1;
        }
    }
    for (i=1; i<=arima->var; i++) {
        // Scan the line to read in data stream
        for ( j=1; j<=arima->history->total; j++) {
            fscanf(fptr, "%lf", &param);
            arima->history->obs[j][i] = param;
        }
        // Detect Ending Mark, there is a double check here.
        if ( !(fscanf(fptr, "%s", field1)) ) {
            errMsg("read","readARIMAStream","Failed to detect stream end.",0); *stat1 = 1; return 1;
        } else {
            if (strncmp(field1, "E", 1)) {
                errMsg("read","readARIMAStream","Wrong stream end indicator.",0); *stat1 = 1; return 1;
            }
        }
    }
    /*************************************************************************************************************/
    
    // Skip this line
    fgets(line, sizeof line, fptr);
    
    /*************** NOISE SECTION: This sections stores all corresponding historical noise data *****************/
    //------>Allocate memory to time series blocks
    if ( !(arima->noise = (streamType *) mem_malloc(sizeof(streamType)))) {
        errMsg("allocation","readARIMAStream","Failed to allocate memeroty to variate blocks",0);
        *stat1 = 1; return 1;
    }
    //------>Collect Data Stream for different variate
    if (fgets(line, sizeof line, fptr) != NULL) {
        sscanf(line, "%s %d %d", field1, &int_param, &int_param2);
    } else {
        errMsg("read","readARIMAStream","Failed to read the next line",0);
        *stat1 = 1; return 1;
    } // Ready to check
    if (strncmp(field1, "NOISE", 5)) {
        errMsg("read","readARIMAStream","Failed to find a bloc in right position", 0);
        *stat1 = 1; return 1;
    }
    printf("%d Historical Noise Data Points Available. Starting index = %d \n", int_param, int_param2);
    arima->noise->start = int_param2;
    arima->noise->past  = int_param2-1;
    arima->noise->diff  = 0;
    arima->noise->open  = int_param2;
    arima->noise->total = int_param;

    //------>Allocate memory to data stream structure
    if ( !(arima->noise->obs = (vector *) arr_alloc(arima->noise->total+1, vector) )) {
        errMsg("allocation","readARIMAStream","failed to allocate memory to data stream.",0);
        *stat1 = 1; return 1;
    }
    //------>Allocate memory to store historical data at each time point for all variates
    for (i=1; i<=arima->noise->total; i++) {
        if (!((arima->noise->obs[i] = arr_alloc(arima->var+1, double)))){
            errMsg("allocation","readARIMAStream","Failed to allocate memory to obs space",0);
            *stat1 = 1; return 1;
        }
    }
    for (i=1; i<=arima->var; i++) {
        // Scan the line to read in data stream
        for ( j=1; j<=arima->noise->total; j++) {
            fscanf(fptr, "%lf", &param);
            arima->noise->obs[j][i] = param;
        }
        // Detect Ending Mark, there is a double check here.
        if ( !(fscanf(fptr, "%s", field1)) ) {
            errMsg("read","readARIMAStream","Failed to detect stream end.",0); *stat1 = 1; return 1;
        } else {
            if (strncmp(field1, "E", 1)) {
                errMsg("read","readARIMAStream","Wrong stream end indicator.",0); *stat1 = 1; return 1;
            }
        }
    }
    /************************************************************************************************************/
    
    // Skip this line
    fgets(line, sizeof line, fptr);
    
    return 0;
}

int readARIMAAux(FILE *fptr, tsType *arima, int *stat1) {
    
    char        line[BLOCKSIZE], field1[NAMESIZE];
    int         int_param, decompLength;
    double      param;
    int         i, j;
    
    /************** Inversion Data Section: Read emperical data hourly mean and standard deviation **************/
    if (arima->transType == 1) {
        errMsg("no implementation", "readARIMAAux", "no implementation for EMP method", 0);
    } else if (arima->transType == 2 || arima->transType == 3) {
        /****** CLASSIC DECOMPOSITION Section: Read emperical data hourly mean and standard deviation ******/
        if (fgets(line, sizeof line, fptr) != NULL) {
            sscanf(line, "%s %d", field1, &int_param);
        } else {
            errMsg("read","read_arma","Failed to read the next line",0);
            *stat1 = 1; return 1;
        } // Ready to check
        // Read Head Line
        if (strncmp(field1, "STL", 3) ) {
            errMsg("read","read_arma","Failed to find EMP SD head line", 0);
            *stat1 = 1; return 1;
        }
        if ( !(arima->decomp = (decomposeType **) arr_alloc(arima->var+1, decomposeType *)) ) {
            errMsg("allocation","readArma","Failed to allocate memory to cla array",0);
            *stat1 = 1; return 1;
        }
        
        decompLength = int_param;
        
        for (i=1; i<=arima->var; i++) {
            if (!(arima->decomp[i] = (decomposeType *) mem_malloc(sizeof(decomposeType)))) {
                errMsg("allocation","readArma","Failed to allocate memory to individual decomposeType for classic",0);
                *stat1 = 1;
                return 1;
            }
            if ( !(arima->decomp[i]->trend = arr_alloc(decompLength + 1, double) ) ) {
                errMsg("allocation","readArma","Failed to allocate memory to individual decomposeType->trend array",0);
                *stat1 = 1;
                return 1;
            }
            if ( !(arima->decomp[i]->seasonal = arr_alloc(decompLength + 1, double) )) {
                errMsg("allocation","readArma","Failed to allocate memory to individual decomposeType->trend array",0);
                *stat1 = 1;
                return 1;
            }
            
            arima->decomp[i]->cnt = decompLength;
            
            if (arima->decomp[i]->cnt > arima->subh*run.HORIZON) {
                printf("VAR[%d] :: More STL/CLA information is provided. \n NEEDED %d. \n Reading the all %d in. Setting cnt to %d\n",
                       i, arima->subh*run.HORIZON, arima->subh*run.HORIZON, arima->decomp[i]->cnt);
            } else {
                printf("VAR[%d] :: Not sufficient STL/CLA information is provided. NEEDED %d. Setting cnt to %d\n",
                       i, arima->subh*run.HORIZON, arima->decomp[i]->cnt);
            }
            
            // Check variate ID for trend data
            fscanf(fptr, "%d", &int_param);
            if (int_param != i) {
                errMsg("read","readArma","wrong stl/cla trend title",0);
                *stat1 = 1;
                return 1;
            }
            // Read Trend Data
            for (j=1; j<=arima->decomp[i]->cnt; j++) {
                fscanf(fptr, "%lf", &param);
                arima->decomp[i]->trend[j] = param;
            }
            // Detect Ending Mark, there is a double check here.
            if ( !(fscanf(fptr, "%s", field1)) ) {
                errMsg("read","read_arma","Failed to detect stream end.",0); *stat1 = 1;
                return 1;
            } else {
                if (strncmp(field1, "E", 1)) {
                    errMsg("read","read_arma","Wrong stream end indicator.",0); *stat1 = 1;
                    return 1;
                }
            }
            // Skip this "\n"
            fgets(line, sizeof line, fptr);
            
            // Check variate ID for seasonality data
            fscanf(fptr, "%d", &int_param);
            if (int_param != i) {
                errMsg("read","readArma","wrong stl/cla trend title",0);
                *stat1 = 1; return 1;
            }
            // Read Seasonal Data
            for (j=1; j<=arima->decomp[i]->cnt; j++) {
                fscanf(fptr, "%lf", &param);
                arima->decomp[i]->seasonal[j] = param;
            }
            // Detect Ending Mark, there is a double check here.
            if ( !(fscanf(fptr, "%s", field1)) ) {
                errMsg("read","read_arma","Failed to detect stream end.",0); *stat1 = 1;
                return 1;
            } else {
                if (strncmp(field1, "E", 1)) {
                    errMsg("read","read_arma","Wrong stream end indicator.",0); *stat1 = 1;
                    return 1;
                }
            }
        }
    } else {
        errMsg("no implmentation","readARIMA","no implementation yet",0);
        *stat1 = 1; return 1;
    }
    
    // Skip this ending mark "\n"
    fgets(line, sizeof line, fptr);
    
    /************************************************************************************************************/
    return 0;
}

/*******************************************************************************
 The noise vector is indexed using time first then variate.
 t = 1 -> x_1, y_1, z_1, ... -> this vector can be used to multiply Phi matrix
 t = 2 -> x_2, y_2, z_2, ... -> x,y,z are different variates
 t = 2 -> x_3, y_3, z_3, ...
 ...
 t = N -> x_N, y_N, z_N, ...
 *******************************************************************************/
int simulateNoise (tsType *arima, vector *white, long long *seed){
    
    int totalN = arima->subh;
    int totalVar = arima->var;
    int i, status = 1;
    
    for (i=1; i<=totalN; i++) {
        status = normal(arima->mean, arima->stdev, totalVar, white[i]+1, seed);
        if ( !status ) {
            errMsg("noise","simulate_noise","Failed to generate noise vector",0);
            return 1;
        }
    }

    return 0;
}

/************************************************************************
 This subroutine translate a observation vector into a real right hand 
 side vector. This vector is assume to be under the order as
 var_1_1, var_1_2, var_1_3, ... , var_1_subh, var_2_1, var_2_2, ..., 
                                            ,..., var_var_subh....
 Make sure subproblem random rows follow the same order.
 ************************************************************************/
void reverseVector(vector observ, tsType *arima, int t) {
    
    int i,j,k=0;
    
    // Notice: Variable/Row Line-up matters ...
    for (j=1; j<=arima->subh; j++) {
        for (i=1; i<=arima->var; i++) {
            observ[k] = reverse(observ[k],i,j,t,arima);
            k++;
        }
    }
    
}


/************************************************************************
 This transfer the time series model output back to real observation of 
 the original data. The method here used is directly emperical 
 distribution method. minPos and hourPos need to be indicated. 
 This subroutine indicates a horizon of a day.
 ------------------------------------------------------------------------
 In case of EMP method...
 Given 10 minutes sub-hourly interval, 24 hour horizon---
 minPos = 1->XX:10; 2->XX:20; 3->XX:30; 4->XX:40; 5->XX:50; 6->XX+1:00;
 hourPos = 1->00:00; 2->01:00; 3->02:00; ... ; 23->22:00; 24->23:00;
 varID => which variate (wind farm);
 ------------------------------------------------------------------------
 In case of STL/CLA method...
 Given 10 minutes sub-hourly interval, 24 hour horizon---
 minPos     = 1->XX:10; 2->XX:20; 3->XX:30; 4->XX:40; 5->XX:50; 6->XX+1:00;
 hourPos    = 
 ************************************************************************/
double reverse(double tsResult, int varID, int minPos, int hourPos, tsType *arima) {
        
    double inv;
    int pos;
    
    if (arima->transType == 0) {
        return tsResult;
    }

    // Check model type: hourly or sub-hourly
    if (arima->subh == 1) {
        minPos = 0;
        pos = hourPos * arima->subh;
    } else {
        pos = (hourPos-1) * arima->subh + minPos;
    }
    
    // Fetch information from time series type for corresponding trends and seasonality
    if (arima->transType == 1) {
        errMsg("no Implementation", "reverse", "no implementation for EMP method", 0);
    } else if (arima->transType == 2 || arima->transType == 3) {
        inv = tsResult + arima->decomp[varID]->trend[pos] + arima->decomp[varID]->seasonal[pos];
        return inv;
    }
    
    return -1;
}


void freeVarBlkType(varBlkType *variate) {
    mem_free(variate->col);
    mem_free(variate->row);
    mem_free(variate->seq);
    mem_free(variate);
}//End of freeVarBlkType

void freeStreamType(streamType *stream, int totalTime) {
    int i;
    for (i=1;i<=totalTime;i++) {
        if (stream->obs[i]) mem_free(stream->obs[i]);
    }
    mem_free(stream->obs);
    mem_free(stream);
}//End of freeStreamType

void freeDecomposeType(decomposeType *decompose) {
    if (decompose->seasonal) mem_free(decompose->seasonal);
    if (decompose->trend)    mem_free(decompose->trend);
    mem_free(decompose);
}//End of freeDecomposeType

void freeTsType(tsType *arima) {
    int i;
    
    
    mem_free(arima->mean);
    mem_free(arima->stdev);

    for (i=1;i<=arima->p;i++)
        if (arima->Phi[i]) freeSparseMatrix(arima->Phi[i]);
    mem_free(arima->Phi);

    
    for (i=1;i<=arima->q;i++)
        if (arima->Theta[i]) freeSparseMatrix(arima->Theta[i]);
    mem_free(arima->Theta);

    
    for (i=1;i<=arima->var;i++)
        if (arima->variate[i]) freeVarBlkType(arima->variate[i]);
    mem_free(arima->variate);
    
    if (arima->history) freeStreamType(arima->history, arima->history->total);
    if (arima->noise) freeStreamType(arima->noise, arima->noise->total);
    
    for (i=1; i<=arima->subh; i++) {
        if (arima->histARObs[i]) mem_free(arima->histARObs[i]);
        if (arima->histMAObs[i]) mem_free(arima->histMAObs[i]);
    }
    mem_free(arima->histARObs);
    mem_free(arima->histMAObs);
    
    if (arima->transType == 2 || arima->transType == 3) {
        for (i=1; i<=arima->var; i++) {
            if (arima->decomp[i]) freeDecomposeType(arima->decomp[i]);
        }
        mem_free(arima->decomp);
    }
    
    mem_free(arima);
}//End of freeTsType
