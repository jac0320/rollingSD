//
//  subprob.h
//  rollingSP
//
//  Created by Site Wang on 4/27/16.
//  Copyright © 2016 Site Wang. All rights reserved.
//

#ifndef subprob_h
#define subprob_h

#include "rollingSP.h"
#include "stocUpdt.h"

/* subroutines in subprob.c */
int solveSubprob(probType *prob, oneProblem *subprob, solnType *soln, vector Xvect, vector observ, double *mubBar, BOOL *subFeasFlag);
vector computeRHS(numType *num, coordType *coord, sparseVector *bBar, sparseMatrix *Cbar, vector X, vector obs);
oneProblem *newSubprob(oneProblem *subprob);
void freeSubprob(oneProblem *subprob);


#endif /* subprob_h */
