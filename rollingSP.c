//
//  rollingSP.c
//
//

#include "rollingSP.h"
#include "algo.h"

// Basic Global Variables
long        MEM_USED = 0;
string      outputDir;
string      evalsolDir;
char        inputDir[BLOCKSIZE];
BOOL        newPeriod = FALSE;

// Configuration Type :: Controls and measures many things in rolling stochastic programming
runType     run;

// Please find more information about rollingSP with the following address
// www.www.www

/*--------------------------------------Main Subroutine-----------------------------------------*/
int main(int argc, char *argv[]){
    
    int status = 0;
    oneProblem 	*orig = NULL;
    timeType 	*tim = NULL;
    stocType 	*stoc = NULL;
    char		probName[NAMESIZE];
    
    /* Print Title */
    printLine(); printf("Build 0731\n");
    
	/* ~~~~~~~~~~~~ Open solver environment ~~~~~~~~~~~~ */
	openSolver();

	/* ~~~~~~~~~~~~ Read algorithm configuration files ~~~~~~~~~~~~ */
    status = readRun(inputDir);
    if (status) {
        errMsg("read","main","failed to read run file", 0);
        goto QUIT;
    }

	/* ~~~~~~~~~~~~ Read problem information ~~~~~~~~~~~~ */
	/* request for problem name to be solved, the path should be provided in the configuration file */
	if ( argc < 5 )
		parseCmdLine(argv, probName, &run.HORIZON, &run.DETERMINISTIC, &run.WARM_UP);
    else {
		strcpy(probName, argv[1]);
        run.HORIZON = (int) str2float(argv[2]);
        run.DETERMINISTIC = (int) str2float(argv[3]);
        run.WARM_UP = (int) str2float(argv[4]);
    }
	/* ~~~~~~~~~~~~ Setup the output directory ~~~~~~~~~~~~ */
	createOutputDir(outputDir, "rollingSP", probName);
    
    /* ~~~~~~~~~~~~ Read problem files ~~~~~~~~~~~~ */
    status = readFiles(inputDir, probName, &orig, &tim, &stoc);
    if (status) {
        errMsg("read","main", "failed to read in the problem file", 0);
        goto QUIT;
    }

    /* ~~~~~~~~~~~~ Begin the algorithm ~~~~~~~~~~~~ */
    //status = setupAlgo(orig,stoc,tim);
    status = algo(orig, stoc, tim);
    if (status) {
        errMsg("solve","main","Failed to solve the problem.",0);
        goto QUIT;
    }

    printf("\n\n################################# Successfully completed ################################## \n");


QUIT:
    freeOneProblem(orig);
    freeStocType(stoc);
    freeTimeType(tim);
    mem_free(outputDir);
    closeSolver();

    return 0;
}//END main()

void parseCmdLine(string *argv, string probName, int *horizon, int *deterministic, int *warmStart) {
    
    printf("Please enter the name of the problem: ");
    scanf("%s", probName);
    printf("Please enter the total horizon: ");
    scanf("%d", horizon);
    printf("Run deterministic problem? (0-NO || 1-YES) ");
    scanf("%d", deterministic);
    printf("Apply warm start techniques? (0-NO || 1-YES) ");
    scanf("%d", warmStart);

}//END parseCmdLine

int readRun(string inputDir) {
    
	int 	status;
	char    param[BLOCKSIZE], comment[2*BLOCKSIZE];
	FILE    *fptr = NULL;

	outputDir = (string) mem_malloc(BLOCKSIZE*sizeof(char));
    evalsolDir = (string) mem_malloc(BLOCKSIZE*sizeof(char));

	/* Assign default values for configuration/run parameters */
    run.DETERMINISTIC = 0;
    run.HORIZON = 1;
    run.RUN_SEED = 3554548844580680;
    run.MIN_ITER = 1;
    run.MAX_ITER = 100;          // Default Maximum Iteration is designed for debugging
    run.CUT_MULT = 5;
    run.TOLERANCE = 0.001;       // Used for condition checking and 0 determination // Nominal Tolerance
    run.R1 = 0.2;
    run.R2 = 0.95;
    run.R3 = 2.0;
    run.MIN_QUAD_SCALAR = 0.001;
    run.MAX_QUAD_SCALAR = 10000.0;
    run.MASTER_TYPE = 5;         // Default MASTERTYPE is QP
    run.TAU = 2;
    run.TEST_TYPE = 1;           // Default setting allows full test
    run.EPSILON = 0.001;         // Nominal Tolerance
    run.SCAN_LEN = 256;          // Newly added, must a power of 2!!
    run.PRE_EPSILON = 0.01;      // Ratio used in preTestI for optimality
    run.PERCENT_PASS = 0.95;
    run.RESAMPLE_SEED = 9495518635394380;
    run.RESERVE_SEED = 9495518635394380;
    run.M = 50;
    run.EVAL_RUN_FLAG = 0;
    run.EVAL_FLAG = 0;
    run.EVAL_MIN_ITER = 1;
    run.EVAL_SEED = 2668655841019641;
    run.EVAL_ERROR = 0.05;
    run.PI_EVAL_START = 1;
    run.BOOTSTRAP_TEST = 1;
    run.CONFID_HI = 1.0;
    run.CONFID_LO = 1.45;
    run.AUTO_SEED = 1;
    run.PRINT_CYCLE = 100;
    run.PI_CYCLE = 1;
    run.PI_EVAL_CLEAR = 1;
    run.PI_EVAL_DELAY = 1;
    run.SOLVER = 0;
    run.DUAL_EPSILON = 0.000002;
    run.WARM_UP = 1;
    run.WARM_START_LB = 0; //This will take some time
    run.SCALAR = 1.0;
    run.EVALUATOR = 0;
	strcpy(inputDir, "./");
	strcpy(outputDir, "./");
    strcpy(evalsolDir, "./");

	/* open the configuration/run file */
	fptr = fopen("./rollingSP.run", "r");
	if (fptr == NULL) {
		errMsg("read","read_run","Failed to open run/configuration file",0);
		return 1;
	}
    
    // Configuration File Provides 39 Input Parameters
	while ((status = (fscanf(fptr, "%s", param) != EOF))) {
		if (!(strcmp(param, "INPUTDIR")))               //0
			fscanf(fptr, "%s", inputDir);
		else if (!(strcmp(param, "OUTPUTDIR")))         //1
			fscanf(fptr, "%s", outputDir);
        else if (!(strcmp(param, "EVALSOLFILE")))        // Additional solution path
            fscanf(fptr, "%s", evalsolDir);
        else if (!(strcmp(param, "EVALUATOR")))        // Additional solution path
            fscanf(fptr, "%d", &run.EVALUATOR);
		else if (!(strcmp(param, "RUN_SEED")))          //2
			fscanf(fptr, "%lld", &run.RUN_SEED);
        else if (!(strcmp(param, "DETERMINISTIC")))     //3
            fscanf(fptr, "%d", &run.DETERMINISTIC);
		else if (!(strcmp(param, "MAX_ITER")))          //4
			fscanf(fptr, "%d", &run.MAX_ITER);
		else if (!(strcmp(param, "MIN_ITER")))          //5
			fscanf(fptr, "%d", &run.MIN_ITER);
		else if (!(strcmp(param, "CUT_MULT")))          //6
			fscanf(fptr, "%d", &run.CUT_MULT);
		else if (!(strcmp(param, "TOLERANCE")))         //7
			fscanf(fptr, "%lf", &run.TOLERANCE);
		else if (!(strcmp(param, "R1")))                //8
			fscanf(fptr, "%lf", &run.R1);
		else if (!(strcmp(param, "R2")))                //9
			fscanf(fptr, "%lf", &run.R2);
		else if (!(strcmp(param, "R3")))                //10
			fscanf(fptr, "%lf", &run.R3);
		else if (!(strcmp(param, "MIN_QUAD_SCALAR")))   //11
			fscanf(fptr, "%lf", &run.MIN_QUAD_SCALAR);
		else if (!(strcmp(param, "MAX_QUAD_SCALAR")))   //12
			fscanf(fptr, "%lf", &run.MAX_QUAD_SCALAR);
		else if (!(strcmp(param, "TAU")))               //13
			fscanf(fptr, "%d", &run.TAU);
		else if (!(strcmp(param, "MASTER_TYPE")))       //14
			fscanf(fptr, "%d", &run.MASTER_TYPE);
		else if (!(strcmp(param, "TEST_TYPE")))         //15
			fscanf(fptr, "%d", &run.TEST_TYPE);
		else if (!(strcmp(param, "EPSILON")))           //16
			fscanf(fptr, "%lf", &run.EPSILON);
		else if (!(strcmp(param, "PRE_EPSILON")))       //17
			fscanf(fptr, "%lf", &run.PRE_EPSILON);
		else if (!(strcmp(param, "PERCENT_PASS")))      //18
			fscanf(fptr, "%lf", &run.PERCENT_PASS);
		else if (!(strcmp(param, "RESAMPLE_SEED")))     //19
			fscanf(fptr, "%lld", &run.RESAMPLE_SEED);
		else if (!(strcmp(param, "M")))                 //20
			fscanf(fptr, "%d", &run.M);
		else if (!(strcmp(param, "EVAL_FLAG")))         //21
			fscanf(fptr, "%d", &run.EVAL_FLAG);
		else if (!(strcmp(param, "EVAL_MIN_ITER")))     //22
			fscanf(fptr, "%d", &run.EVAL_MIN_ITER);
		else if (!(strcmp(param, "EVAL_SEED")))         //23
			fscanf(fptr, "%lld", &run.EVAL_SEED);
		else if (!(strcmp(param, "EVAL_ERROR")))        //24
			fscanf(fptr, "%lf", &run.EVAL_ERROR);
        else if (!(strcmp(param, "SCAN_LEN")))          //25
            fscanf(fptr, "%d", &run.SCAN_LEN);
        else if (!(strcmp(param, "EVAL_RUN_FLAG")))     //26
            fscanf(fptr, "%d", &run.EVAL_RUN_FLAG);
        else if (!(strcmp(param, "PI_EVAL_START")))     //27
            fscanf(fptr, "%d", &run.PI_EVAL_START);
        else if (!(strcmp(param, "PRINT_CYCLE")))       //28
            fscanf(fptr, "%d", &run.PRINT_CYCLE);
        else if (!(strcmp(param, "BOOTSTRAP_TEST")))    //29
            fscanf(fptr, "%d", &run.BOOTSTRAP_TEST);
        else if (!(strcmp(param, "CONFID_LO")))         //30
            fscanf(fptr, "%lf", &run.CONFID_LO);
        else if (!(strcmp(param, "CONFID_HI")))         //31
            fscanf(fptr, "%lf", &run.CONFID_HI);
        else if (!(strcmp(param, "AUTO_SEED")))         //32
            fscanf(fptr, "%d", &run.AUTO_SEED);
        else if (!(strcmp(param, "PI_CYCLE")))          //33
            fscanf(fptr, "%d", &run.PI_CYCLE);
        else if (!(strcmp(param, "SOLVER")))            //34
            fscanf(fptr, "%d", &run.SOLVER);
        else if (!(strcmp(param, "DUAL_EPSILON")))      //35
            fscanf(fptr, "%lf", &run.DUAL_EPSILON);
        else if (!(strcmp(param, "PI_EVAL_CLEAR")))     //36
            fscanf(fptr, "%d", &run.PI_EVAL_CLEAR);
        else if (!(strcmp(param, "PERCENT_WARM")))      //37
            fscanf(fptr, "%lf", &run.PERCENT_WARM);
        else if (!(strcmp(param, "BOOTSTRAP_VALID")))   //38
            fscanf(fptr, "%d", &run.BOOTSTRAP_VALID);
        else if (!(strcmp(param, "WARM_UP")))           //39
            fscanf(fptr, "%d", &run.WARM_UP);
        else if (!(strcmp(param, "PI_EVAL_DELAY")))     //40
            fscanf(fptr, "%d", &run.PI_EVAL_DELAY);
        else if (!(strcmp(param, "WARM_START_LB")))     //41
            fscanf(fptr, "%d", &run.WARM_START_LB);
        else if (!(strcmp(param, "SCALAR")))            //42
            fscanf(fptr, "%lf", &run.SCALAR);
		else if (!strcmp(param, "//"))
			fgets(comment, 2*BLOCKSIZE, fptr);
		else {
			printf ("%s\n", param);
			errMsg("read", "readConfig", "unrecognized parameter in configuration file", 1);
		}
	}

	fclose(fptr);
	return 0;
}//END read_run()
