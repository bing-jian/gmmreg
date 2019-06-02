#ifndef FGT_UTILS_H_
#define FGT_UTILS_H_

// Code for computing the direct Gauss transform using a kd-tree.

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

#include "nanoflann.hpp"

using namespace nanoflann;

namespace gmmreg {

/// The default fgt epsilon.
const double DEFAULT_EPSILON = 1e-4;
const double DEFAULT_EPSILON_COEFF = 3.034854258770293;
// std::sqrt(std::log(1.0 / DEFAULT_EPSILON))

template <typename T>
using PointCloud = vnl_matrix<T>;

// This is the "dataset to kd-tree" adaptor class:
template <typename Derived>
struct PointCloudAdaptor {
  typedef typename Derived::element_type coord_t;

  const Derived& obj;  //!< A const ref to the data set origin

  /// The constructor that sets the data set source
  PointCloudAdaptor(const Derived& obj_) : obj(obj_) {}

  /// CRTP helper method
  inline const Derived& derived() const { return obj; }

  // Must return the number of data points
  inline size_t kdtree_get_point_count() const {
    // return derived().pts.size();
    return derived().rows();
  }

  // Returns the dim'th component of the idx'th point in the class:
  // Since this is inlined and the "dim" argument is typically an immediate
  // value, the
  //  "if/else's" are actually solved at compile time.
  inline coord_t kdtree_get_pt(const size_t idx, const size_t dim) const {
    return derived()(idx, dim);
  }

  // Optional bounding-box computation: return false to default to a standard
  // bbox computation loop.
  //   Return true if the BBOX was already computed by the class and returned in
  //   "bb" so it can be avoided to redo it again.
  //   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3
  //   for point clouds)
  template <class BBOX>
  bool kdtree_get_bbox(BBOX& /*bb*/) const {
    return false;
  }

};  // end of PointCloudAdaptor

template <typename T>
using MatrixAdaptor = PointCloudAdaptor<PointCloud<T>>;

template <typename T>
using MyKdTree =
    KDTreeSingleIndexAdaptor<L2_Simple_Adaptor<T, MatrixAdaptor<T>>,
                             MatrixAdaptor<T>>;

template <typename T>
struct NanoflannTree {
  NanoflannTree(const PointCloud<T>& source)
      : matrix_adaptor({source}),
        tree(int(source.cols()), matrix_adaptor,
             KDTreeSingleIndexAdaptorParams(10 /* max leaf */)) {}

  NanoflannTree(const NanoflannTree&) = delete;
  NanoflannTree& operator=(const NanoflannTree&) = delete;
  NanoflannTree& operator=(NanoflannTree&&) = delete;

  MatrixAdaptor<T> matrix_adaptor;
  MyKdTree<T> tree;
};

template <typename T>
T FastGaussTransform(const NanoflannTree<T>& tree, const PointCloud<T>& moving,
                     T scale, vnl_matrix<T>& gradient) {
  T* grad = gradient.data_block();

  int m = moving.rows();
  int dim = moving.cols();

  T cross_term = 0;
  #pragma omp parallel for
  for (int i = 0; i < m * dim; ++i) {
    grad[i] = 0;
  }
  const T* A = moving.data_block();
  const T* B = tree.matrix_adaptor.derived().data_block();
  int n = tree.matrix_adaptor.derived().rows();

  T cutoff_radius = scale * DEFAULT_EPSILON_COEFF;
  T r2 = cutoff_radius * cutoff_radius;

  T h2 = scale * scale;
  T inv_h2 = 1.0 / h2;

  nanoflann::SearchParams params;
  params.sorted = false;

  #pragma omp parallel for reduction(+ : cross_term)
  for (int j = 0; j < m; ++j) {
    std::vector<std::pair<size_t, T>> indices_distances;
    indices_distances.reserve(unsigned(n));

    size_t nfound =
        tree.tree.radiusSearch(&A[j * dim], r2, indices_distances, params);
    T cost_ij = 0;
    for (size_t i = 0; i < nfound; ++i) {
      const auto entry = indices_distances[i];
      cost_ij = std::exp(-entry.second * inv_h2);
      for (int d = 0; d < dim; ++d) {
        grad[j * dim + d] -=
            cost_ij * 2.0 * (A[j * dim + d] - B[entry.first * dim + d]);
      }
      cross_term += cost_ij;
    }
  }
  h2 *= m * n;
  inv_h2 = 1.0 / h2;

  #pragma omp parallel for
  for (int i = 0; i < m * dim; ++i) {
    grad[i] *= inv_h2;
  }
  return cross_term / (m * n);
}

template <typename T>
void FastNeighborSearch(const NanoflannTree<T>& tree,
                        const vnl_matrix<T>& moving, T scale,
                        std::vector<std::pair<int, int>>* edges) {
  edges->clear();
  int m = moving.rows();
  edges->reserve(10 * m);
  int dim = moving.cols();

  const T* A = moving.data_block();
  const T* B = tree.matrix_adaptor.derived().data_block();
  int n = tree.matrix_adaptor.derived().rows();

  T cutoff_radius = scale * DEFAULT_EPSILON_COEFF;
  T r2 = cutoff_radius * cutoff_radius;

  nanoflann::SearchParams params;
  params.sorted = false;

  #pragma omp for
  for (int j = 0; j < m; ++j) {
    std::vector<std::pair<size_t, T>> indices_distances;
    indices_distances.reserve(unsigned(n));

    size_t nfound =
        tree.tree.radiusSearch(&A[j * dim], r2, indices_distances, params);
    for (size_t i = 0; i < nfound; ++i) {
      const auto entry = indices_distances[i];
      if (entry.first < m && entry.first >= 0) {
        // TODO(bing.jian):
        // investigate why sometimes entry.first is out of range.
        edges->emplace_back(j, entry.first);
      }
    }
  }
}

template <typename T>
T FastSelfGaussTransform(const vnl_matrix<T>& pts,
                         const std::vector<std::pair<int, int>>& edges, T scale,
                         vnl_matrix<T>& gradient) {
  T* grad = gradient.data_block();

  int m = pts.rows();
  int dim = pts.cols();

  int total_dim = m * dim;
  T cross_term = 0;
  for (int i = 0; i < total_dim; ++i) {
    grad[i] = 0;
  }

  T h2 = scale * scale;
  T inv_h2 = 1.0 / h2;

  const T* A = pts.data_block();

  int total_edges = edges.size();
  #pragma omp parallel for reduction(+ : cross_term)
  for (int k = 0; k < total_edges; ++k) {
    int i = edges[k].first;
    int j = edges[k].second;
    T dist_ij = 0;
    for (int d = 0; d < dim; ++d) {
      dist_ij +=
          (A[i * dim + d] - A[j * dim + d]) * (A[i * dim + d] - A[j * dim + d]);
    }
    T cost_ij = std::exp(-dist_ij * inv_h2);
    for (int d = 0; d < dim; ++d) {
      grad[i * dim + d] -= cost_ij * (A[i * dim + d] - A[j * dim + d]);
      grad[j * dim + d] -= cost_ij * (A[j * dim + d] - A[i * dim + d]);
    }
    cross_term += cost_ij;
  }
  h2 *= m * m;
  inv_h2 = 1.0 / h2;
  #pragma omp parallel for
  for (int i = 0; i < m * dim; ++i) {
    grad[i] *= inv_h2;
  }
  return cross_term / (m * m);
}

}  // namespace gmmreg

#endif  // FGT_UTILS_H_
