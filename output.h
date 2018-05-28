//
//  output.h
//  rollingSP
//
//  Created by Site Wang on 4/23/16.
//  Copyright © 2016 Site Wang. All rights reserved.
//

#ifndef out_h
#define out_h


#include "rollingSP.h"
#include "cell.h"

void printPeriodStat(probType **prob, solnType *soln, cellType *cell, FILE *out, FILE *outII, FILE *outIII);
void printInitInfo(oneProblem *orig, timeType *tim, probType **prob, cellType *cell);
void rSPprint(FILE *out, const char *s);
void printCut(cutType *cuts, numType *num, int idx);
void printRunningTime(solnType *soln);

#endif /* out_h */
