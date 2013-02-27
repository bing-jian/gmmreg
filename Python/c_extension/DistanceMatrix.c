/* ----------------------------------------------------------
 * $RCSfile: DistanceMatrix.c,v $
 * Language:  C
 * $Author: bing.jian $
 * $Date: 2008-11-06 01:51:00 -0500 (Thu, 06 Nov 2008) $
 * $Revision: 96 $
 * ---------------------------------------------------------- */

#include <stdlib.h>
#include <math.h>

/* PLEASE READ THIS HEADER, IMPORTANT INFORMATION INSIDE */
#include "memory_layout_note.h"

#ifdef WIN32
__declspec( dllexport )
#endif
void squared_distance_matrix(const double*X, const double*Y, const double*g, int m, int n, int d, double* dist)
{
    int i,c,r,p,k;
    double *temp, *dist_ij;
    
    temp = (double*)malloc(d*sizeof(double));
    /* memset(temp,0,d*sizeof(double)); */
    dist_ij = (double*)malloc(d*sizeof(double));

    for(i=0, r=0;r<m;++r) {
	    for(c=0;c<n;++c,++i){
            for(p=0;p<d;++p){
                temp[p] = 0;
            }
            for(p=0;p<d;++p){
                dist_ij[p] = Y[c*d+p]-X[r*d+p];
                for (k=0;k<d;++k){
                    temp[k] += dist_ij[p]*g[p*d+k];
                }
            }
            for (p=0;p<d;++p){
                dist[i] += temp[p]*dist_ij[p];
            }
        }
    }
    /* clean up */
    free(dist_ij);
    free(temp);
}

#ifdef WIN32
__declspec( dllexport )
#endif
void distance_matrix(const double*X, const double*Y, const double*g, int m, int n, int d, double* dist)
{
    int i;
    squared_distance_matrix(X, Y, g, m, n, d, dist);
    for (i=0;i<m*n;++i){
        dist[i] = sqrt(dist[i]);
    }
}
