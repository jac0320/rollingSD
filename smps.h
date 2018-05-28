//
//  smps.h
//  rollingSP
//


#ifndef SMPS_H_
#define SMPS_H_

#include "solver.h"

typedef struct{
    int		type;			/* type of problem: LP, QP, MIP or MIQP */
    void	*lp;			/* problem pointer to be used by solver */
    string	name;			/* name of the problem */
    int		objsen;			/* sense of the objective: 1 for minimization and -1 for maximization */
    int		mac;			/* number of columns */
    int 	mar;			/* number of rows */
    int		numBin;			/* number of binary variables in the problem */
    int		numInt;			/* number of integer variables in the problem */
    int		numnz;			/* number of non-zero elements in constraint matrix */
    vector	objx;			/* objective function coefficients */
    vector	rhsx;			/* right-hand side */
    string	senx;			/* constraint sense */
    intvec	matbeg;			/* sparse matrix representation: column beginning */
    intvec	matcnt;			/* sparse matrix representation: number of non-zero entries in a column */
    intvec	matind;			/* sparse matrix representation: rows with non-zero entries */
    vector	matval;			/* sparse matrix representation: non-zero coefficients of the matrix */
    vector	bdl;			/* lower bound */
    vector	bdu;			/* upper bound */
    string	ctype;			/* type of decision variables: 'C' continuous, 'B' binary, 'I' general integer, 'S' semi-continuous, 'N' semi-integer */
    string	objname;		/* objective function name */
    int		rstorsz;		/* memory size for storing row names */
    string	*rname;			/* vector of row names */
    string	rstore;			/* row names string */
    int		cstorsz;		/* memory size for storing column names */
    string	*cname;			/* vector of column names */
    string	cstore;			/* column name string */
    int		macsz;			/* extended column size */
    int		marsz;			/* extended row size */
    int		matsz;			/* extended matrix size */
}oneProblem;

typedef struct {
    int		type;			/* type of time file declaration, 0 for implicit and 1 for explicit */
    string	probName;		/* name of the problem as read from time file */
    int		numStages;	    /* number of stages in the problem */
    string	*stgNames;		/* unique strings to identify stages*/
    intvec	row;			/* a list of row names which mark the beginning of a new stage */
    intvec	col;			/* a list of column names which mark the beginning of a new stage */
    int		numRows;		/* used with explicit time file declaration only, set to numStages in implicit declaration */
    intvec	rowStg;			/* used with explicit time file declaration only */
    int		numCols;		/* used with explicit time file declaration only, set to numStages in implicit declaration */
    intvec	colStg;			/* used with explicit time file declaration only */
}timeType;

typedef struct {
    string	type;
    BOOL	sim;
    int		numOmega; 			/* number of stochastic elements stored in structure */
    int 	numCipher; 			/* number of ints needed to encode an observation */
    int		numGroups;
    intvec	row; 				/* row number array in the original problem; -1 indicates objective function */
    intvec	col; 				/* column number array in the original problem; -1 indicates right-hand side */
    intvec	numVals;			/* number of realization for each random variable */
    vector	*vals; 				/* indexed array of discrete realizations of random variable */
    vector	*probs;				/* indexed array of probabilities associated with discrete realizations*/
    intvec	numPerGroup;        /* TODO: (HG) What is this? When is this used? */
    intvec	groupBeg;
    vector	mean;         		/* mean of each rv */
}stocType;

/* subroutines in input.c */
int readFiles(string inputDir, string probName, oneProblem **orig, timeType **tim, stocType **stoc);
oneProblem *readCore(string inputDir, string probName);
timeType *readTime(string inputDir, string probName, oneProblem *orig);
stocType *readStoc(string inputDir, string probName, oneProblem *orig);
int readIndep(FILE *fptr, string *fields, oneProblem *orig, int maxOmegas, int maxVals, stocType *stoc);
int readBlocks(FILE *fptr, string *fields, oneProblem *orig, int maxOmegas, int maxVals, stocType *stoc);
int readBlk(FILE *fptr, string *fields, oneProblem *orig, int maxOmegas, int maxVals, BOOL origRV, stocType *stoc);

void freeOneProblem(oneProblem *p);
void freeTimeType(timeType *tim);
void freeStocType(stocType *stoc);

#endif /* SMPS_H_ */
