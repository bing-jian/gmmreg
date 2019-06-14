#include "downsample_utils.h"

#include <unordered_set>
#include <vector>

namespace gmmreg {

/* Python version

def create_sparse_nodes(pts, radius):
     N = pts.shape[0]
     tree = KDTree(pts)
     available = set(range(N))
     assigned = set([])
     removed = set([])
     while(available):
          rand_id = np.random.randint(0, len(available))
          idx = list(available)[rand_id]
          assigned.add(idx)
          available.remove(idx)
          near_pts = tree.query_ball_point(pts[idx], r=radius)
          for pt_idx in near_pts:
               if pt_idx in available:
                    removed.add(pt_idx)
                    available.remove(pt_idx)
     return assigned, available, removed

 *
 */
template <typename T>
vnl_matrix<T> CreateSparseNodes(const NanoflannTree<T>& tree, T dist_radius) {
  const PointCloud<T>& pts = tree.matrix_adaptor.derived();
  int n = pts.rows();
  int dim = pts.cols();
  std::unordered_set<int> available;
  std::vector<int> assigned;
  available.reserve(n);
  assigned.reserve(n);
  for (int i = 0; i < n ; ++i) {
    available.insert(i);
  }

  nanoflann::SearchParams params;
  params.sorted = false;
  const T* A = pts.data_block();
  T r2 = dist_radius * dist_radius;
  while(!available.empty()) {
    int idx = *available.begin();
    assigned.push_back(idx);
    available.erase(idx);

    std::vector<std::pair<size_t, T>> indices_distances;
    indices_distances.reserve(unsigned(n));

    size_t nfound =
        tree.tree.radiusSearch(&A[idx * dim], r2, indices_distances, params);
    for (size_t i = 0; i < nfound; ++i) {
      const auto entry = indices_distances[i];
      int pt_idx = entry.first;
      if (available.find(pt_idx) != available.end()) {
          available.erase(pt_idx);
      }
    }
  }
  vnl_matrix<T> pts_kept;
  int n_kept = assigned.size();
  pts_kept.set_size(n_kept, dim);
  for (int i = 0; i < n_kept; ++i) {
    pts_kept.set_row(i, pts.get_row(assigned[i]));
  }
  return pts_kept;
}

}  // namespace gmmreg
