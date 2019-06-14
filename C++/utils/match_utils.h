#ifndef GMMREG_UTILS_MATCH_UTILS_H_
#define GMMREG_UTILS_MATCH_UTILS_H_

#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>

namespace gmmreg {

template <typename T>
void ComputeSquaredDistanceMatrix(const vnl_matrix<T>& A,
                                  const vnl_matrix<T>& B, vnl_matrix<T>& D);


template<typename T>
int SelectPoints(const vnl_matrix<T>& pts,
    const std::vector<int>& index, vnl_matrix<T>& selected);

template<typename T>
void PickIndices(const vnl_matrix<T>& dist,
    std::vector<int>& row_index, std::vector<int>& col_index);

template<typename T>
void FindNearestPairs(const vnl_matrix<T>& dist,
    vnl_matrix<int>& pairs, const T& threshold);


template<typename T>
int FindWorkingPairs(const vnl_matrix<T>& M,
                     const vnl_matrix<T>& S,
                     const vnl_matrix<T>& Transformed_M,
                     const T& threshold,
                     vnl_matrix<T>& working_M,
                     vnl_matrix<T>& working_S);

}  // namespace gmmreg

#include "match_utils.cc"

#endif // GMMREG_UTILS_MATCH_UTILS_H_
