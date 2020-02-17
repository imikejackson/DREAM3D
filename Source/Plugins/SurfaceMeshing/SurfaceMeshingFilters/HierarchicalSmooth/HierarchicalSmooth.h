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
 *
 * HierarchicalSmooth.h -- Contains main and auxiliary routine prototypes
 * for smoothing according to the hierarchical filter method.
 *
 */

#pragma once

#include "HierarchicalSmoothTypes.h"

namespace HierarchicalSmooth
{

enum class Type
{
  Serial,
  Cyclic
};

SparseMatrixF laplacian2D(int32_t N, Type type = Type::Serial);
std::tuple<SparseMatrixF, std::vector<int32_t>> graphLaplacian(const TriMesh& tri);

MeshNode smooth(const MeshNode& nodes, Type type = Type::Serial, float threshold = 0.001f, uint64_t iterations = 53);
MeshNode smooth(const MeshNode& nodes, const MatIndex& nFixed, const SparseMatrixF& GL, float threshold = 0.001f, uint64_t iterations = 53);

std::tuple<SparseMatrixF, SparseMatrixF> getDirichletBVP(const SparseMatrixF& GL, const SparseMatrixF& y, const MatIndex& nFixed, const MatIndex& nMobile);
std::tuple<SparseMatrixF, SparseMatrixF> analyzeLaplacian(const SparseMatrixF& GL);

float getObjFn(Smoother& smth, float feps, const SparseMatrixF& fSmallEye, const SparseMatrixF& LTL, const SparseMatrixF& LTK, const SparseMatrixF& Data, const MatIndex& nMobile,
                const SparseMatrixF& yMobile, const SparseMatrixF& D, const SparseMatrixF& AyIn, SparseMatrixF& yOut);

} // namespace HSmoothMain
