//
//  algo.h
//  rollingSP
//
//  Created by Site Wang on 7/18/16.
//  Copyright © 2016 Site Wang. All rights reserved.
//

#ifndef algo_h
#define algo_h

#include "rollingSP.h"
#include "cell.h"

cellType *setup_rSP(oneProblem *orig, stocType *stoc, timeType *tim, probType ***prob, solnType **soln);
int algo(oneProblem *orig, stocType *stoc, timeType *tim);
void freeSingleProbtype (probType *prob, int stageCnt);

#endif /* algo_h */
