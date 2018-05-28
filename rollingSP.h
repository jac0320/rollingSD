//
//  rollingSP.h
//
//

#ifndef rollingSP_h
#define rollingSP_h

#include <stdio.h>
#include "utils.h"
#include "solver.h"
#include "smps.h"
#include "prob.h"

#define DROPPED     0
#define SURVIVED    1

#undef  SETUP_CHECK         //Print all setup information
#undef  SOLVE_CHECK         //Check all things related to CPLEX
#undef  DETAIL_CHECK        //Infomration generated at each iteration, literaly everything

#undef  OFF
#undef  INPUT
#undef  VALID
#undef  LB_VALID
#undef  ROLL_VALID

#undef  TRACE               //Tells every single detailed calculation of the process
#undef  OMEGA               //Record evaluated scenarios
#undef  REPORT

#undef  ENFORCE             //Speed up the warm-start process
#undef  CUTS_SIGMA

#undef  LAMBDA
#undef  TUNE
#undef  EVAL
#undef  TSVALID
#undef  SUBOBJ
#undef  CUTS
#undef  ROLL
#undef  ROLL_DETAIL
#undef  DETAIL              //(Implementation not completed. Code will run very slow)Further detials of the code
#undef  OPTIMAL             //Write down important information in optimality check section
#undef  DEBUGGING           //Turn the debugging check before solving master problem
#undef  MODIFY              //Write all generated problems during the algorithm
#undef  DETAIL_OUT          //Write detailed conditions of each iteration
#define NOW                 //Debug check, this should always on

/* Determines the maximum number of cuts allowed, given the number of primal variables in the master problem. */
//#define		MAX_CUTS(c)		(3 * (c) + 3)
#define MAX_CUTS(c)     (1 *(c) + 3)
#define TRIVIAL			0
#define	NONTRIVIAL		1

//TODO: [SITE] (Comments on the structure)
typedef struct {
    int			HORIZON;
    int         DETERMINISTIC;
	long long 	RUN_SEED;
	int			MIN_ITER;
	int			MAX_ITER;
	int			CUT_MULT;
	double		TOLERANCE;
	double		R1;
	double		R2;
	double		R3;
	double		MIN_QUAD_SCALAR;
	double		MAX_QUAD_SCALAR;
	int 		MASTER_TYPE;
	int 		TAU;
	int         TEST_TYPE;          // 0->PRE-TEST ONLY || 1->FULL-TEST ALLOWED
	double		EPSILON;
    int         SCAN_LEN;
	double		PRE_EPSILON;
	double		PERCENT_PASS;
	long long	RESAMPLE_SEED;
    long long   RESERVE_SEED;
	int			M;
    int         EVAL_RUN_FLAG;
	int         EVAL_FLAG;
	long long	EVAL_SEED;
	int			EVAL_MIN_ITER;
	double		EVAL_ERROR;
    int         PRINT_CYCLE;
    int         BOOTSTRAP_TEST;
    double      CONFID_LO;
    double      CONFID_HI;
    int         AUTO_SEED;
    int         PI_EVAL_CLEAR;
    int         PI_EVAL_START;
    int         PI_EVAL_DELAY;
    int         PI_CYCLE;
    int         SOLVER;
    double      DUAL_EPSILON;
    double      PERCENT_WARM;
    int         BOOTSTRAP_VALID;
    int         WARM_UP;
    int         WARM_START_LB;
    int         EVALUATOR;
    double      SCALAR;
}runType;

/* To save time and space, Pi x b and Pi x C are calculated as soon as possible and stored in structures like sigma and delta.  Toward
 * this end, pixbCType represents a single calculation of pi X b (which is a scalar) and pi X C (which is a vector).*/
typedef struct{
	double 	b;
	vector 	C;
}pixbCType;

/* The lambda structure stores some of the dual variable values from every distinct dual vector obtained during the program.  Each vector contains
 * only those dual variables whose corresponding rows in the subproblem constraint matrix contain random elements.  _val_ is an array of
 * these dual vectors (thus it is 2-D). _row_ gives the corresponding row number for a given dual variable in _val_.  _cnt_ represents
 * the number of dual vectors currently stored in lambda. */
typedef struct {
	int 	cnt;
    intvec  alive;
	vector 	*vals;
}lambdaType;

/* The sigma matrix contains the values of Pi x bBar and Pi x Cbar  for all values of pi obtained so far (note it does not depend on
 * observations of omega).  _col_ gives the column number of each non-zero element in pi X Cbar.  _val_ is an array of values
 * for pi X bBar and pi X Cbar, one entry for each pi.  Note that values  which are always zero (because Rbar or Cbar is zero there) are not
 * stored.  The _lamb_ array is the same size as the _val_ array, and for each element in _val_ the corresponding element in _lamb_ references
 * the dual vector in lambda that was used to calculate that entry in sigma. */
typedef struct {
	int 		cnt;
    int         aCnt;               /* Actual count of the active sigma values */
    intvec      linker;
	pixbCType 	*vals;
	intvec		lambdaIdx;
	intvec		ck; 				/* record the iteration # when sigma was created */
}sigmaType;

/* When calculating istar for a cut, it is useful to have two separate references into the sigma and delta structures, since each dual vector
 * is stored in two places -- part in sigma and part in delta.  The final entry in cut->istar[] will just be the _sigma_ field of this structure. */
typedef struct {
	int     delta;
	int     sigma;
    double  height;
}iType;

/* Each cut consists of a scalar value for the right-hand side, stored in _alpha_, a vector of coefficients for the master program's primal
 * variables, stored in _beta_, and an array of indices to the maximal pi for each observation of omega, stored in _istar_ (these are references
 * into sigma).  In order to weight the cuts properly, _cut_obs_ gives the number of samples on which the given cut was based.  In contrast,
 * _omega_cnt_ gives the number of observations used to generate the cut. Cuts which are "loose" when the master is solved will be dropped periodically,
 * based on the _slack_cnt_ field, which counts the number of consecutive iterations the cut's constraint was slack.  Finally, in order to interface
 * with the LP solver, each cut should know what its row number is in the master constraint matrix, so this is stored in _row_num_
 * Note that all cuts in a given cell will have the same "num_obs" (from Sen's earlier implementation paper) so this is stored in theta. */
typedef struct {
	double 	alpha;					/* scalar value for the right-hand side */
	vector 	beta;					/* coefficients for the master program's primal variables */
    int     cutPeriod;
	int 	cutObs;					/* number of samples on which the given cut was based */
	int 	omegaCnt;				/* number of *distinct* observations on which the cut was based (this is also the length of istar) */
	intvec	iStar;					/* indices of maximal pi for each distinct observation */
	BOOL	isIncumb;				/* indicates if the cut is an incumbent cut */
	double 	alphaIncumb;			/* right-hand side when using QP master, this is useful for quick updates */
	int 	slackCnt;				/* number of times a cut has been slack, used in deciding when the cut needs to be dropped */
	int 	rowNum;					/* row number for master problem in solver */
}oneCut;

/* A collection of the single cuts described above is stored here. The _val_ array holds pointers to cut structures, while the _cnt_ field tells
 * how many cuts are currently stored in _val_. */
typedef struct {
	int 	cnt;
	oneCut 	**val;
}cutType;

/* The delta matrix contains the values of lambda_pi X bOmega and lambda_pi X Comega for all values of pi and all observations of omega.
 * _col_ gives the column number of the each non-zero element in the multiplication of lambda_pi X Comega (the same elements are non-zero
 * each time).  _val_ is an array of vectors of (lambda_pi X bOmega, lambda_pi X Comega) pairs of calculations.  A row in _val_ corresponds
 * to a distinct dual vector, and a column in _val_ corresponds to a distinct observation of omega.  Thus, every pi-omega combination is
 * represented here, and the size of the delta matrix can be determine ffrom lambda->cnt and omega->cnt.
 * Note that when elements of omega get dropped, vacant columns appear in delta.  This is ok, but be sure to loop carefully! */
typedef struct {
	pixbCType 	**val;
}deltaType;

/* TODO: Update this part */
typedef struct {
	int     cnt;
	intvec  weight;
    intvec  sigmaIdx;
	vector	*vals;
}omegaType;

/* Aux structure for ARIMA models. Similar to pixbCType but without pi. Can be used both by ARIMA model and gammaType. */
typedef struct{
    double b;
    vector C;
}bCType;

typedef struct{
    double  allTime;
    int     detSolveCnt;
    int     masterCnt;
    int     subprobCnt;
    int     cutGenCnt;
    int     dropCutCnt;
    int     optCnt;
    int     argmaxCnt;
    // Periodic Level
    vector  warmUpTime;  //Time for warming up
    vector  solveTime;     //Time for each time period
    vector  evalOptTime; //Time for optimality evaluation
    // Iteration Level
    intvec  iterIndex;
    vector  masterTime;         //Time for solving master problem
    vector  optTime;            //Time for optimality check
    vector  fullCheck_sample;   //Time for fullCheck to sample
    vector  fullCheck_reform;   //Time for fullCheck to reform cuts
    vector  fullCheck_lb;       //Time for fullCheck to evaluate lower bounds
    // Solving Cnt Level
    intvec  lpCntIndex;
    vector  subprobTime; //Time for solving subprobTime
    vector  cutGenTime;  //Time for cut generation
    vector  argmaxTime;  //Time for conudcting the argmax operation
    vector  stoStrucTime;//Time for updating the stochastic structures
}runTimeType;

/* Whereas the prob structure contain relatively permanent data, the cell structure keeps track of the structural changes in the problem and the soln
 * structure contains temporary data necessary for solving a single SD problem. To avoid allocating and freeing vectors every iteration, the soln
 * structure also contains a number of solution vectors (_Pi_, _candid_x_ and _incumb_x_) allocated once for the duration of the solution.
 * The _candidEst_, _incumbEst_, and _incumbStdev_ fields store info about each of the solution vectors, for use in controlling the algorithm.
 * _iCutupdt_ stores the iteration of the last time the incumbent cut was updated.  _gamma_ stores the expected (in the previous iteration)
 * improvement of the objective function, and is used in deciding whether or not to update the incumbent solution. The structure also holds information
 * (normDk_1 and normDk) used for updating the scaling factor in QP. */
typedef struct {
	double      optValM;
    vector      meanX;
	vector      piS;           		/* dual solutions for subproblem */
    vector      piSRatio;            /* TODO: [SITE] newly added part for calculate dual stability */
	double      mubBar;				/* dual slack information for subproblem */
	vector      piM;           		/* dual solutions for master */
	vector      djM;				/* reduced cost vector for master */
	vector      candidX;			/* candidate master solution */
	double      candidEst;			/* estimate at the candidate solution */
	vector      incumbX;			/* incumbent master solution */
	double      incumbEst;			/* estimate at incumbent solution */
	double      incumbStdev;		/* standard deviation of incumbent estimate */
	BOOL        incumbChg;			/* set to be true if the incumbent solution has changed in an iteration */
	int         iCutIdx;			/* index of incumbent cut in cell->cuts structure */
    int         cCutIdx;            /* */
	int         iCutUpdt;			/* iteration number when incumbent cut is updated */
	double      gamma;				/* improvement in objective function value */
	BOOL        optFlag;			/* optimality flag */
	double      normDk_1;			/* (\Delta x^{k-1})^2 */
	double      normDk;				/* (\Delta x^k)^2 */
	BOOL        subFeasFlag;		/* flag indicates infeasible subproblem */
	BOOL        infeasIncumb;		/* flag indicates if incumbent solution is infeasible */
	BOOL        optMode;
    BOOL        dualStableFlag;     /* Indicator whether dual multiplier is stable or not */
    int         dualStableIter;     /* Indicator of which iteration dual solution is stablized */
    BOOL        preCheckFlag;       /* Indicator of pre test I flag */
    BOOL        preCheckEverFlag;   /* Indicator of pre test I flag ever passed */
    double      FTError;            /* Aerage error of failed replication in full test */ //TODO: Need to be initialized
    int         repPassed;          /* This record the replcation passed on optimal solution when conducting full test */
    runTimeType *runTime;
} solnType;

void parseCmdLine(string *argv, string probName, int *horizon, int *deterministic, int *warmStart);
int readRun(string inputDir);

// Setup.c
int setupAlgo(oneProblem *orig, stocType *stoc, timeType *tim);

#endif /* rollingSP_h */
