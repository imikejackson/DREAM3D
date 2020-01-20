/*
 This code was extracted from HierarchicalSmooth, https://github.com/siddharth-maddali/HierarchicalSmooth,
 from commit 31ef680011f4bbaef59c0944876e84222ea7c49f
*/

// Copyright (c) 2016-2018, Siddharth Maddali
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of Carnegie Mellon University nor the names of its
//    contributors may be used to endorse or promote products derived from
//    this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

/*
 * Types.h -- contains all user-defined types.
 */

#pragma once

#include <cstdint>
#include <functional>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "Eigen/Eigen"
#include "Eigen/Sparse"
#include "Eigen/IterativeLinearSolvers"

namespace HierarchicalSmooth
{
/*
 * TriMesh:
 * Eigen array of integer triplets; the prototype of Delaunay
 * triangulations in this library.
 */
using TriMesh = Eigen::Array<uint64_t, Eigen::Dynamic, 3, Eigen::RowMajor>;

/*
 * MeshNode:
 * Eigen array of float triplets, each row representing
 * a 3D cartesian mesh node.
 */
using MeshNode = Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor>;

/*
 * FaceLabel:
 * This data type is analogous to Dream.3D's FaceLabels property
 * and represents a grain boundary patch by specifying the
 * grain IDs on either side of the patch.
 */
using FaceLabel = Eigen::Array<int32_t, Eigen::Dynamic, 2, Eigen::RowMajor>;

/*
 * NodeType:
 * Dream.3D-specific dataset which indicates the type of node
 * in a surface mesh: interior, triple junction or quad junction
 * ( denoted by 2, 3, 4 respectively if on the interior and
 * 12, 13, 14 if on the volume surface.
 */
using NodeType = Eigen::Array<int8_t, Eigen::Dynamic, 1>;

/*
 * IsSmoothed:
 * Boolean array specifying whether each node has been smoothed
 * or not. At the beginning, only the NodeTypes 4 and 14 should
 * be initialized to true, the others should be false.
 */
using IsSmoothed = Eigen::Array<bool, Eigen::Dynamic, 1>;

/*
 * MatIndex:
 * Special typedef of an Eigen vector of integers to indicate
 * an array of indices, to be used in slicing with libigl.
 * Also, need to set data type int and not int because of
 * compatibility with libigl.
 */
using MatIndex = Eigen::Matrix<int, Eigen::Dynamic, 1>;

/*
 * EdgePair, EdgeList:
 * Bookkeeping devices for edges in a Delaunay mesh, each
 * edge being represented by an ordered pair of integers.
 */
using EdgePair = std::pair<int32_t, int32_t>;
using EdgeList = std::vector<EdgePair>;

/*
 * EdgePairEqual:
 * Comparison function for two different EdgePairs, used in custom hash function.
 */
struct EdgePairEqual
{
  bool operator()(const EdgePair& lhs, const EdgePair& rhs) const
  {
    return (lhs.first == rhs.first) && (lhs.second == rhs.second);
  }
};

/*
 * Customized dictionary that takes EdgePairs as keys.
 */
template <typename T>
struct DictBase
{
  using EdgeDict = std::unordered_map<EdgePair, T, std::hash<EdgePair>, EdgePairEqual>;
};

/*
 * The dictionary initialization happens like this:
 *
 * DictBase< your-type >::EdgeDict my_dict;
 *
 */

/*
 * SparseMatrixF:
 * Shorthand for Eigen's sparse matrix type.
 *
 * TripletF:
 * Triplet containing a position indices for
 * a single sparse matrix element (i, j ) and
 * the floating point value at that position.
 * Defined in Eigen/Sparse
 *
 * SparseMatrixB:
 * Boolean mask for type SparseMatrixF
 *
 * TripletB:
 * Equivalent of T for boolean variables.
 */
using SparseMatrixF = Eigen::SparseMatrix<float>;
using TripletF = Eigen::Triplet<float>;
using SparseMatrixB = Eigen::SparseMatrix<bool>;
using TripletB = Eigen::Triplet<bool>;

/*
 * Typedef for the conjugate gradient solver for sparse systems.
 */
using Smoother = Eigen::ConjugateGradient<SparseMatrixF, Eigen::Upper | Eigen::Lower>;

} // namespace HierarchicalSmooth

/*
 * The actual hash function
 */
namespace std
{
template <>
struct hash<HierarchicalSmooth::EdgePair>
{
  std::size_t operator()(const HierarchicalSmooth::EdgePair& EP) const
  {
    return std::hash<std::size_t>{}(static_cast<std::size_t>(EP.first) << (sizeof(std::size_t) * 4) | static_cast<std::size_t>(EP.second));
  }
};
} // namespace std
