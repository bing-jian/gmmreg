#ifndef vnl_gmmreg_h_
#define vnl_gmmreg_h_

#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>
#include <vector>

#ifdef WIN32
#include <windows.h>
#endif

#define BUFSIZE 4096

int LoadMatrixFromTxt(const char* filename, vnl_matrix<double>& matrix);

double GaussTransform(const vnl_matrix<double>& A,
                      const vnl_matrix<double>& B,
                      double scale);
double GaussTransform(const vnl_matrix<double>& A,
                      const vnl_matrix<double>& B,
                      double scale,
                      vnl_matrix<double>& gradient);
//double GaussTransform(const vnl_matrix<double>& A, const vnl_matrix<double>& B, double scale, double* gradient);
//double GaussTransform(const double* A, const double* B, int m, int n, int dim, double scale, double* grad);
void ComputeTPSKernel(const vnl_matrix<double>& model,
                      const vnl_matrix<double>& ctrl_pts,
                      vnl_matrix<double>& U,
                      vnl_matrix<double>& K);
void ComputeGaussianKernel(const vnl_matrix<double>& model,
                           const vnl_matrix<double>& ctrl_pts,
                           vnl_matrix<double>& G,
                           vnl_matrix<double>& K,
                           double beta);

void parse_tokens(char* str, const char delims[],
                  std::vector<double>& v_tokens);
void parse_tokens(char* str, const char delims[],
                  std::vector<int>& v_tokens);

void ComputeSquaredDistanceMatrix(const vnl_matrix<double>& A,
                                  const vnl_matrix<double>& B,
                                  vnl_matrix<double>& D);
void pick_indices(const vnl_matrix<double>& dist,
                  vnl_matrix<int>& pairs,
                  const double threshold);
int find_working_pair(const vnl_matrix<double>& M,
                      const vnl_matrix<double>& S,
                      const vnl_matrix<double>& Transformed_M,
                      const double threshold,
                      vnl_matrix<double>& working_M,
                      vnl_matrix<double>& working_S);

int get_config_fullpath(const char* input_config,char* f_config);
void save_matrix( const char * filename, const vnl_matrix<double>& x);
void save_vector( const char * filename, const vnl_vector<double>& x);

//void normalize(const vnl_matrix<double>& x, vnl_vector<double>& centroid, double& scale, vnl_matrix<double>& normalized_x);
//void denormalize(const vnl_matrix<double>& x, const vnl_vector<double>& centroid, const double scale, vnl_matrix<double>& denormalized_x);
void normalize(vnl_matrix<double>& x,
               vnl_vector<double>& centroid,
               double& scale);
void denormalize(vnl_matrix<double>& x,
                 const vnl_vector<double>& centroid,
                 const double scale);
void compute_P(const vnl_matrix<double>& x,
               const vnl_matrix<double>& y,
               vnl_matrix<double>& P, double &E, double sigma, int outliers);

void quaternion2rotation(vnl_vector<double> q,
                         vnl_matrix<double>& R,
                         vnl_matrix<double>& g1,
                         vnl_matrix<double>& g2,
                         vnl_matrix<double>& g3,
                         vnl_matrix<double>& g4);

void quaternion2rotation(vnl_vector<double> q, vnl_matrix<double>& R);

#ifndef WIN32
char *strupr(char *string);
char *strlwr(char *string);
#endif

#endif //#ifndef vnl_gmmreg_h_
