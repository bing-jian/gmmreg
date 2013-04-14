#include <vcl_iostream.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

#define SQR(X)  ((X)*(X))
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>

#include "gmmreg_utils.h"

int LoadMatrixFromTxt(const char* filename, vnl_matrix<double>& matrix) {
  std::ifstream infile(filename, std::ios_base::in);
  if (infile.is_open()) {
    if (matrix.read_ascii(infile)) {
      return matrix.rows();
    } else {
      std::cerr << "unable to parse input file " << filename
                << " as a matrix." << std::endl;
      return -1;
    }
  } else {
    std::cerr << "unable to open model file " << filename << std::endl;
    return -1;
  }
}

/*
 *  Note: The input point set containing 'n' points in 'd'-dimensional
 *  space should be arranged in the memory such that the j-th coordinate of i-th point
 *  is at location  (i*d + j), i.e., the input matrix in MATLAB should
 *  be of size (d, n) since MATLAB is column-major while the input array
 *  in NumPy should be of size (n, d) since NumPy is row-major.
 */

double GaussTransform(const double* A, const double* B,
    int m, int n, int dim, double scale) {
  double cross_term = 0;
  scale = SQR(scale);
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      double dist_ij = 0;
      for (int d = 0; d < dim; ++d) {
        dist_ij += SQR(A[i * dim + d] - B[j * dim + d]);
      }
      double cost_ij = exp(-1.0 * dist_ij / scale);
      cross_term += cost_ij;
    }
    /* printf("cross_term = %.3f\n", cross_term);  */
  }
  return cross_term / (m * n);
}

void GaussianAffinityMatrix(const double* A, const double* B,
    int m, int n, int dim, double scale, double* dist) {
  scale = -2.0*SQR(scale);
  for (int i = 0, k = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j, ++k) {
      double dist_ij = 0;
      for (int d = 0; d < dim; ++d) {
        dist_ij +=  SQR(A[i * dim + d] - B[j * dim + d]);
      }
      dist[k] = exp(dist_ij / scale);
    }
  }
}

double GaussTransform(const double* A, const double* B,
    int m, int n, int dim, double scale, double* grad) {
  double cross_term = 0;
  for (int i = 0; i < m * dim; ++i) {
    grad[i] = 0;
  }

  scale = SQR(scale);
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      double dist_ij = 0;
      for (int d = 0; d < dim; ++d) {
        dist_ij +=  SQR(A[i * dim + d] - B[j * dim + d]);
      }
      double cost_ij = exp(-1.0 * dist_ij / scale);
      for (int d = 0; d < dim; ++d){
        grad[i * dim + d] -= cost_ij * 2.0 * (A[i * dim + d] - B[j * dim + d]);
      }
      cross_term += cost_ij;
    }
    /* printf("cross_term = %.3f\n", cross_term);  */
  }
  scale *= m * n;
  for (int i = 0; i < m * dim; ++i) {
    grad[i] /= scale;
  }
  return cross_term / (m * n);
}

double GaussTransform(const vnl_matrix<double>& A,
    const vnl_matrix<double>& B, double scale) {
  // assert A.cols() == B.cols()
  return GaussTransform(A.data_block(), B.data_block(),
      A.rows(), B.rows(), A.cols(), scale);
}

double GaussTransform(const vnl_matrix<double>& A,
    const vnl_matrix<double>& B, double scale,
    vnl_matrix<double>& gradient) {
  // assert A.cols() == B.cols()
  return GaussTransform(A.data_block(), B.data_block(),
      A.rows(), B.rows(), A.cols(), scale,
      gradient.data_block());
}

// todo: add one more version when the model is same as ctrl_pts
// reference:  Landmark-based Image Analysis, Karl Rohr, p195
void ComputeTPSKernel(const vnl_matrix<double>& model,
    const vnl_matrix<double>& ctrl_pts,
    vnl_matrix<double>& U, vnl_matrix<double>& K) {
  int m = model.rows();
  int n = ctrl_pts.rows();
  int d = ctrl_pts.cols();
  //asssert(model.cols()==d==(2|3));
  K.set_size(n, n);
  K.fill(0);
  U.set_size(m, n);
  U.fill(0);
  double eps = 1e-006;

  vnl_vector<double> v_ij;
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      v_ij = model.get_row(i) - ctrl_pts.get_row(j);
      if (d == 2) {
        double r = v_ij.squared_magnitude();
        if (r > eps) {
          U(i, j) = r * log(r) / 2;
        }
      } else if (d == 3) {
        double r = v_ij.two_norm();
        U(i, j) = -r;
      }
    }
  }
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      v_ij = ctrl_pts.get_row(i) - ctrl_pts.get_row(j);
      if (d == 2) {
        double r = v_ij.squared_magnitude();
        if (r > eps) {
          K(i, j) = r * log(r) / 2;
        }
      }
      else if (d == 3) {
        double r = v_ij.two_norm();
        K(i, j) = -r;
      }
    }
  }
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      K(i, j) = K(j, i);
    }
  }

}

/*
   Matlab code in cpd_G.m:
   k=-2*beta^2;
   [n, d]=size(x);  [m, d]=size(y);

   G=repmat(x,[1 1 m])-permute(repmat(y,[1 1 n]),[3 2 1]);
   G=squeeze(sum(G.^2,2));
   G=G/k;
   G=exp(G);
   */
void ComputeGaussianKernel(const vnl_matrix<double>& model,
    const vnl_matrix<double>& ctrl_pts,
    vnl_matrix<double>& G, vnl_matrix<double>& K,
    double lambda) {
  int m,n,d;
  m = model.rows();
  n = ctrl_pts.rows();
  d = ctrl_pts.cols();
  //asssert(model.cols()==d);
  //assert(lambda>0);

  G.set_size(m,n);
  GaussianAffinityMatrix(model.data_block(), ctrl_pts.data_block(),
      m, n, d, lambda, G.data_block());

  if (model == ctrl_pts) {
    K = G;
  } else {
    K.set_size(n,n);
    GaussianAffinityMatrix(ctrl_pts.data_block(), ctrl_pts.data_block(),
        n, n, d, lambda, K.data_block());
  }
}

void ComputeSquaredDistanceMatrix(const vnl_matrix<double>& A,
    const vnl_matrix<double>& B, vnl_matrix<double>& D) {
  int m = A.rows();
  int n = B.rows();
  //asssert(A.cols()==B.cols());
  D.set_size(m,n);
  vnl_vector<double> v_ij;
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      v_ij = A.get_row(i) - B.get_row(j);
      D(i,j) = v_ij.squared_magnitude();
    }
  }
}

void parse_tokens(char* str, const char delims[],
    std::vector<double>& v_tokens) {
  char* pch = strtok (str, delims);
  while (pch != NULL) {
    v_tokens.push_back(atof(pch));
    pch = strtok (NULL, delims);
  }
}

void parse_tokens(char* str, const char delims[],
    std::vector<int>& v_tokens) {
  char* pch = strtok (str, delims);
  while (pch != NULL) {
    v_tokens.push_back(atoi(pch));
    pch = strtok (NULL, delims);
  }
}

int select_points(const vnl_matrix<double>&pts,
    const std::vector<int>&index, vnl_matrix<double>& selected) {
  int n = index.size();
  int d = pts.cols();
  selected.set_size(n,d);
  for (int i = 0; i < n; ++i) {
    selected.update(pts.extract(1, d, index[i]), i);
  }
  return n;
}

void pick_indices(const vnl_matrix<double>&dist,
    std::vector<int>&row_index, std::vector<int>&col_index,
    const double threshold) {
  int m = dist.rows();
  int n = dist.cols();
  vnl_vector<int> row_flag, col_flag;
  col_flag.set_size(n);  col_flag.fill(0);
  row_flag.set_size(n);  row_flag.fill(0);
  for (int i = 0; i < m; ++i) {
    double min_dist = dist.get_row(i).min_value();
    if (min_dist < threshold) {
      for (int j = 0; j < n; ++j){
        if (dist(i,j) == min_dist && col_flag[j] == 0){
          row_index.push_back(i);
          row_flag[i] = 1;
          col_index.push_back(j);
          col_flag[j] = 1;
        }
      }
    }
  }
}

void pick_indices(const vnl_matrix<double>&dist,
    vnl_matrix<int>&pairs, const double threshold) {
  int m = dist.rows();
  int n = dist.cols();
  vnl_vector<int> row_flag, col_flag;
  col_flag.set_size(n);  col_flag.fill(0);
  row_flag.set_size(n);  row_flag.fill(0);
  std::vector<int> row_index,col_index;
  for (int i = 0; i < m; ++i) {
    double min_dist = dist.get_row(i).min_value();
    if (min_dist < threshold) {
      for (int j = 0; j < n; ++j) {
        if (dist(i,j)==min_dist && col_flag[j] == 0){
          row_index.push_back(i);
          row_flag[i] = 1;
          col_index.push_back(j);
          col_flag[j] = 1;
        }
      }
    }
  }
  pairs.set_size(2, row_index.size());
  for (int i = 0; i<pairs.cols(); ++i){
    pairs(0,i) = row_index[i];
    pairs(1,i) = col_index[i];
  }
}

int find_working_pair(const vnl_matrix<double>&M, const vnl_matrix<double>&S,
    const vnl_matrix<double>&Transformed_M, const double threshold,
    vnl_matrix<double>&working_M, vnl_matrix<double>&working_S) {
  vnl_matrix<double> dist;
  ComputeSquaredDistanceMatrix(Transformed_M, S, dist);
  std::vector<int> row_index, col_index;
  pick_indices(dist, row_index, col_index, threshold);
  select_points(M, row_index, working_M);
  select_points(S, col_index, working_S);
  return row_index.size();
}

int get_config_fullpath(const char* input_config,char* f_config) {
#ifdef WIN32
  char* lpPart[BUFSIZE]={NULL};
  int retval = GetFullPathName(input_config,
      BUFSIZE, f_config,
      lpPart);

  if (retval == 0) {
    // Handle an error condition.
    printf ("GetFullPathName failed (%d)\n", GetLastError());
    return -1;
  }
  else {
    //printf("The full path name is:  %s\n", f_config);
  }
#else
  strcpy(f_config,input_config);
#endif
  return 0;
}


void save_matrix( const char * filename, const vnl_matrix<double>& x) {
  if (strlen(filename)>0) {
    std::ofstream outfile(filename,std::ios_base::out);
    x.print(outfile);
  }
}

void save_vector( const char * filename, const vnl_vector<double>& x) {
  if (strlen(filename)>0) {
    std::ofstream outfile(filename,std::ios_base::out);
    outfile << x;
  }
}

void normalize(vnl_matrix<double>& x,
    vnl_vector<double>& centroid, double& scale) {
  int n = x.rows();
  if (n==0) return;
  int d = x.cols();
  centroid.set_size(d);

  vnl_vector<double> col;
  for (int i = 0; i < d; ++i) {
    col = x.get_column(i);
    centroid(i) = col.mean();
  }
  for (int i = 0; i < n; ++i) {
    x.set_row(i, x.get_row(i) - centroid);
  }
  scale = x.frobenius_norm() / sqrt(double(n));
  x = x / scale;
}

void denormalize(vnl_matrix<double>& x, const vnl_vector<double>& centroid, const double scale) {
  int n = x.rows();
  if (n==0) return;
  int d = x.cols();
  for (int i = 0; i < n; ++i) {
    x.set_row(i, x.get_row(i) * scale + centroid);
  }
}

void compute_P(const vnl_matrix<double>& x,const vnl_matrix<double>& y, vnl_matrix<double>& P, double &E, double sigma, int outliers) {
  double k;
  k = -2*sigma*sigma;

  //P.set_size(m,n); P.fill(0);
  //vnl_vector<double> v_ij;

  vnl_vector<double> column_sum;
  int m = x.rows();
  int s = y.rows();
  int d = x.cols();
  column_sum.set_size(s);
  column_sum.fill(0);
  double outlier_term = outliers*pow((2*sigma*sigma*3.1415926),0.5*d);
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < s; ++j) {
      double r = 0;
      for (int t = 0; t < d; ++t) {
        r += (x(i,t) - y(j,t))*(x(i,t) - y(j,t));
      }
      P(i,j) = exp(r/k);
      column_sum[j]+=P(i,j);
    }
  }


  if (outliers!=0) {
    for (int i = 0; i < s; ++i)
      column_sum[i] += outlier_term;
  }
  if (column_sum.min_value()>(1e-12)) {
    E = 0;
    for (int i = 0; i < s; ++i) {
      for (int j = 0; j < m; ++j){
        P(j,i) = P(j,i)/column_sum[i];
      }
      E -= log(column_sum[i]);
    }
    //vcl_cerr< < s;
    //vcl_cerr<<P.get_column(10);
  }
  else {
    P.empty();
  }
}


//template<class T>
void quaternion2rotation(vnl_vector<double> q, vnl_matrix<double>& R, vnl_matrix<double>& g1, vnl_matrix<double>& g2, vnl_matrix<double>& g3, vnl_matrix<double>& g4){
  double x,y,z,r;
  double x2,y2,z2,r2;
  x = q[0];  y = q[1];  z=q[2];  r = q[3];
  x2 = q[0] * q[0];
  y2 = q[1] * q[1];
  z2 = q[2] * q[2];
  r2 = q[3] * q[3];
  // fill diagonal terms
  R(0,0) = r2 + x2 - y2 - z2;
  R(1,1) = r2 - x2 + y2 - z2;
  R(2,2) = r2 - x2 - y2 + z2;
  // fill off diagonal terms
  R(0,1) = 2 * (x*y + r*z);
  R(0,2) = 2 * (z*x - r*y);
  R(1,2) = 2 * (y*z + r*x);
  R(1,0) = 2 * (x*y - r*z);
  R(2,0) = 2 * (z*x + r*y);
  R(2,1) = 2 * (y*z - r*x);
  double ss = (x2+y2+z2+r2);
  R = R/ss;
  double ssss = ss*ss;

  // derivative of R(0,0) = r2 + x2 - y2 - z2;
  g1(0,0) = 4*x*(y2+z2)/ssss;  g2(0,0) = -4*y*(x2+r2)/ssss;
  g3(0,0) = -4*z*(x2+r2)/ssss;  g4(0,0) = 4*r*(y2+z2)/ssss;
  // derivative of R(1,1) = r2 - x2 + y2 - z2;
  g1(1,1) = -4*x*(y2+r2)/ssss;  g2(1,1) = 4*y*(x2+z2)/ssss;
  g3(1,1) = -4*z*(y2+r2)/ssss;  g4(1,1) = 4*r*(x2+z2)/ssss;
  // derivative of R(2,2) = r2 - x2 - y2 + z2;
  g1(2,2) = -4*x*(z2+r2)/ssss;  g2(2,2) = -4*y*(r2+z2)/ssss;
  g3(2,2) = 4*z*(x2+y2)/ssss;  g4(2,2) = 4*r*(x2+y2)/ssss;

  // fill off diagonal terms
  // derivative of R(0,1) = 2 * (xy + rz);
  g1(0,1) = 2*y/ss - 2*x*R(0,1)/ssss;
  g2(0,1) = 2*x/ss - 2*y*R(0,1)/ssss;
  g3(0,1) = 2*r/ss - 2*z*R(0,1)/ssss;
  g4(0,1) = 2*z/ss - 2*r*R(0,1)/ssss;
  // derivative of R(0,2) = 2 * (zx - ry);
  g1(0,2) = 2*z/ss - 2*x*R(0,2)/ssss;
  g2(0,2) = -2*r/ss - 2*y*R(0,2)/ssss;
  g3(0,2) = 2*x/ss - 2*z*R(0,2)/ssss;
  g4(0,2) = -2*y/ss - 2*r*R(0,2)/ssss;
  // derivative of R(1,2) = 2 * (yz + rx);
  g1(1,2) = 2*r/ss - 2*x*R(1,2)/ssss;
  g2(1,2) = 2*z/ss - 2*y*R(1,2)/ssss;
  g3(1,2) = 2*y/ss - 2*z*R(1,2)/ssss;
  g4(1,2) = 2*x/ss - 2*r*R(1,2)/ssss;
  // derivative of R(1,0) = 2 * (xy - rz);
  g1(1,0) = 2*y/ss - 2*x*R(1,0)/ssss;
  g2(1,0) = 2*x/ss - 2*y*R(1,0)/ssss;
  g3(1,0) = -2*r/ss - 2*z*R(1,0)/ssss;
  g4(1,0) = -2*z/ss - 2*r*R(1,0)/ssss;
  // derivative of R(2,0) = 2 * (zx + ry);
  g1(2,0) = 2*z/ss - 2*x*R(2,0)/ssss;
  g2(2,0) = 2*r/ss - 2*y*R(2,0)/ssss;
  g3(2,0) = 2*x/ss - 2*z*R(2,0)/ssss;
  g4(2,0) = 2*y/ss - 2*r*R(2,0)/ssss;
  // derivative of R(2,1) = 2 * (yz - rx);
  g1(2,1) = -2*r/ss - 2*x*R(2,1)/ssss;
  g2(2,1) = 2*z/ss - 2*y*R(2,1)/ssss;
  g3(2,1) = 2*y/ss - 2*z*R(2,1)/ssss;
  g4(2,1) = -2*x/ss - 2*r*R(2,1)/ssss;
}


//template<class T>
void quaternion2rotation(vnl_vector<double> q, vnl_matrix<double>& R){
  double x,y,z,r;
  double x2,y2,z2,r2;
  x = q[0];  y = q[1];  z=q[2];  r = q[3];
  x2 = q[0] * q[0];
  y2 = q[1] * q[1];
  z2 = q[2] * q[2];
  r2 = q[3] * q[3];
  // fill diagonal terms
  R(0,0) = r2 + x2 - y2 - z2;
  R(1,1) = r2 - x2 + y2 - z2;
  R(2,2) = r2 - x2 - y2 + z2;
  // fill off diagonal terms
  R(0,1) = 2 * (x*y + r*z);
  R(0,2) = 2 * (z*x - r*y);
  R(1,2) = 2 * (y*z + r*x);
  R(1,0) = 2 * (x*y - r*z);
  R(2,0) = 2 * (z*x + r*y);
  R(2,1) = 2 * (y*z - r*x);
  double ss = (x2+y2+z2+r2);
  R = R/ss;
}

#ifndef WIN32
char *strupr(char *string) {
  char *s;
  if (string) {
    for (s = string;  *s;  ++s) {
      *s = toupper(*s);
    }
  }
  return string;
}

char *strlwr(char *string) {
  char *s;
  if (string) {
    for (s = string;  *s;  ++s) {
      *s = tolower(*s);
    }
  }
  return string;
}
#endif
