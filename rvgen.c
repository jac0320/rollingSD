//
//  rvgen.c
//  rollingSP
//
//

#include <stdio.h>
#include "utils.h"
#include "rollingSP.h"
#include "cell.h"


void generateOmega(cellType *cell, stocType *stoc, vector residual, vector observ, vector auxObs, vector trimObs, vector inputNoise, long long *seed) {

#ifdef TRACE
    printf("\t~generateOmega()\n");
#endif

    // resolve numtiple numGroups issue
    int n, m, offset = 0;
    
    for ( n = 0; n < stoc->numGroups; n++ ) {
        /* -------------- CATEGORY A : INDEPENDENT STOC FILE -------------- */
        if ( strstr(stoc->type, "INDEP") != NULL ) {
            if ( strstr(stoc->type, "DISCRETE") != NULL )
                generateIndep(stoc, observ+offset, n, seed);
            else if ( strstr(stoc->type, "NORMAL") != NULL )
                normal(stoc->mean+offset, stoc->vals[0]+offset, stoc->numPerGroup[n], observ+offset, seed);
            else
                errMsg("rvGeneration", "generateOmega", "random number generation for input type is missing",0);
            for ( m = 0 ; m < stoc->numOmega ; m++ ) {
                residual[m+offset] = observ[m+offset] - stoc->mean[m];
                observ[m+offset] = observ[m+offset] - stoc->mean[m];
                auxObs[m] = 0.0;
            }
            offset += stoc->numPerGroup[n];
        }
        /* -------------- CATEGORY B : BLOCK DISCRETE STOC FILE -------------- */
        else if ( strstr(stoc->type, "BLOCKS") != NULL) {
            if ( strstr(stoc->type, "DISCRETE") != NULL)
                generateBlocks(stoc, observ+offset, n, seed);
            else
                printf("Not yet\n");
            for ( m = 0; m < stoc->numOmega; m++ ) {
                residual[m+offset] = observ[m+offset] - stoc->mean[m];
                observ[m+offset] = observ[m+offset] - stoc->mean[m];
                auxObs[m] = 0.0;
            }
            offset += stoc->numPerGroup[n];
        }
        /* -------------- CATEGORY C : DISTRIBUTION GENERATION -------------- */
        else if ( strstr(stoc->type, "DISTRIB") != NULL ) {
            //~~~~>generateDistrib(observ+offset, n, seed);
            offset += stoc->numPerGroup[n];
            errMsg("no implementation","generateOmega","Not implemented feature.",0);
        }
        /* -------------- CATEGORY D : TIME SERIES MODEL -------------- */
        else if ( strstr(stoc->type, "ARIMA") != NULL) {
            //~~~~>Generate ARIMA observation vector, residual vector, and auxiluary observation for subhourly models
            generateARIMA(cell->t, cell, observ, residual, auxObs, trimObs, inputNoise, seed);
        }
        else
            errMsg("rvgen", "generateOmega", "unknown section type in omegastuff", 0);
    }
    
}//END generateOmega()

void generateBlocks(stocType *stoc, vector observ, int groupID, long long *seed) {
    int 	n, m;
    double 	val, cumm;
    
    /* select which block */
    val = scalit(0,1, seed);
    cumm = 0;
    for ( m = 0; val > cumm; m++ )
        cumm += stoc->probs[groupID][m];
    
    /* read block realizations */
    for (n = 0; n < stoc->numPerGroup[groupID]; n++ )
        observ[n] = stoc->vals[stoc->groupBeg[groupID]+n][m-1];
    
}//END generateBlocks()

void generateIndep(stocType *stoc, vector observ, int groupID, long long *seed) {
    double 	val, cumm;
    int		n, m;
    
    for ( n = 0; n < stoc->numPerGroup[groupID]; n++) {
        val = scalit(0, 1, seed);
        cumm = 0;
        for ( m = 0; val > cumm; ++m )
            cumm += stoc->probs[stoc->groupBeg[groupID] + n][m];
        observ[stoc->groupBeg[groupID]+n] = stoc->vals[stoc->groupBeg[groupID]+n][m-1];
    }
    
}//END generateIndep()

/* The following inverse normal variate generator was published by Micheal J. Wichura, University of Chicago
 * in Applied Statistics, as Algorithm AS 241.  The C function normal() was converted from the Fortran function
 * PPND7 and produces normal random variates for the lower tail of a normal distribution accurate to approx.
 * 7 significant figures. */

int normal(vector mu, vector stdev, int numOmega, vector observ, long long *seed) {
    int i;
    float zero, one, half, split1, split2, const1, const2, a0, a1, a2, a3, b1;
    float b2, b3, c0, c1, c2, c3, d1, d2, e0, e1, e2, e3, f1, f2, p, q, r;
    float endval;
    
    for (i = 0; i < numOmega; i++) {
        p = scalit(0, 1, seed);
        
        zero = 0.0;
        one = 1.0;
        half = one / 2.0;
        split1 = 0.425;
        split2 = 5.0;
        const1 = 0.180625;
        const2 = 1.6;
        
        /* coefficients for p close to 1/2 */
        a0 = 3.3871327179;
        a1 = 50.434271938;
        a2 = 159.29113202;
        a3 = 59.109374720;
        b1 = 17.895169469;
        b2 = 78.775757664;
        b3 = 67.18756360;
        
        /* coefficients for p neither close to 1/2 nor 0 or 1 */
        c0 = 1.4234372777;
        c1 = 2.7568153900;
        c2 = 1.3067284816;
        c3 = .17023821103;
        d1 = .73700164250;
        d2 = .12021132975;
        
        /* coefficients for p near 0 or 1 */
        e0 = 6.6579051150;
        e1 = 3.0812263860;
        e2 = .42868294337;
        e3 = .017337203997;
        f1 = .24197894225;
        f2 = .012258202635;
        
        q = p - half;
        
        if (fabs(q) <= split1) {
            r = const1 - q * q;
            endval = q * (((a3 * r + a2) * r + a1) * r + a0)
            / (((b3 * r + b2) * r + b1) * r + one);
            observ[i] = mu[i] + stdev[i] * endval;
            continue;
        }
        
        if (q < 0.0)
            r = p;
        else
            r = one - p;
        
        if (r <= zero)
            return 0;
        
        r = sqrt(-log(r));
        
        if (r <= split2) {
            r = r - const2;
            endval = (((c3 * r + c2) * r + c1) * r + c0)
            / ((d2 * r + d1) * r + one);
            observ[i] = endval;
        }
        else {
            r = r - split2;
            endval = (((e3 * r + e2) * r + e1) * r + e0)
            / ((f2 * r + f1) * r + one);
            observ[i] = endval;
        }
        if (q < 0)
            observ[i] = -1 * observ[i];
        observ[i] = mu[i] + stdev[i] * observ[i];
    }
    
    return (1);
}//normal()

float scalit(float lower, float upper, long long *seed) {
    float val, wide;
    
    wide = upper - lower;
    val = randUniform(seed);
    
    return ((wide * val) + lower);
}//END scalit()

float randUniform(long long *SEED) {
    static int lo_bits, hi_bits;
    
    lo_bits = ((*SEED) & 0xFFFFL) * 16807;
    hi_bits = (int) (((*SEED) >> 16) * 16807) + (lo_bits >> 16);
    *SEED = ((lo_bits & 0xFFFFL) - 0x7FFFFFFFL) + ((hi_bits & 0x7FFFL) << 16)
    + (hi_bits >> 15);
    return ((*SEED) < 0 ? ((*SEED) += 0x7FFFFFFFL) : (*SEED)) * 4.656612875E-10;
}//END randUniform()

int randInteger(long long *SEED, int iMax) {
    
    return (int) (randUniform(SEED)*iMax);
    
}//END randInteger()
