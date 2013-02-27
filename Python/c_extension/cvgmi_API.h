#ifndef CVGMI_API_H
#define CVGMI_API_H


#ifdef __cplusplus
extern "C" {
#endif

double GaussTransform(const double* A, const double* B, int m, int n, int dim, double scale, double* grad);
void squared_distance_matrix(const double*X, const double*Y, const double*g, int m, int n, int d, double* dist);
void distance_matrix(const double*X, const double*Y, const double*g, int m, int n, int d, double* dist);

#ifdef __cplusplus
}
#endif

#endif
