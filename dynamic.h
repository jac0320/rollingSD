//
//  dynamic.h
//  rollingSP
//
//  Created by Site Wang on 4/27/16.
//  Copyright © 2016 Site Wang. All rights reserved.
//

#ifndef dynamic_h
#define dynamic_h

#include "rollingSP.h"
#include "cell.h"

/* Subroutines in dynamic.c */
int readDyn(string inputDir, string probName, cellType *cell, probType **prob, int stageCnt);
int readBM(FILE *fptr, string *fields, probType **prob);
int readAV(FILE *fptr, string *fields, probType **prob, int *aSeg);
int readDynInit(FILE *fptr, string inputDir, string *fields, string probName, BOOL *initFlag,
                probType **prob, cellType *cell, int stageCnt, int aSeg);

#endif /* dynamic_h */
