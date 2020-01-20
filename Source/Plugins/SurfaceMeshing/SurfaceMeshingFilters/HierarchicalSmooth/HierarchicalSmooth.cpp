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
 * HierarchicalSmooth.cpp -- implementation of the routines declared in HierarchicalSmooth.h
 *
 */

#include "HierarchicalSmooth.h"

#include <cmath>

#include <QtCore/QDebug>

#include "Base.h"
#include "Slice.h"

//============================================================================================

SparseMatrixF HSmoothMain::laplacian2D(int32_t N, Type type)
{
  std::vector<TripletF> tripletList;
  tripletList.reserve(3 * N); // approx. number of nonzero elements
  for(int32_t i = 0; i < N; i++)
  {
    tripletList.push_back(TripletF(i, i, -1.0f));
    if(i != N - 1)
    {
      tripletList.push_back(TripletF(i, i + 1, 1.0f));
    }
  }
  SparseMatrixF temp(N, N);
  temp.setFromTriplets(tripletList.cbegin(), tripletList.cend());
  SparseMatrixF L = SparseMatrixF(temp.transpose()) * temp;

  switch(type)
  {
  case Type::Serial:
  {
    L.coeffRef(N - 1, N - 1) = 1.0f;
  }
  break;
  case Type::Cyclic:
  {
    L.coeffRef(0, 0) = 2.0;
    L.coeffRef(0, N - 1) = -1.0f;
    L.coeffRef(N - 1, 0) = -1.0f;
  }
  break;
  default:
  {
    qDebug() << "HSmoothMain::Laplacian2D: Unrecognized type.";
  }
  break;
  }

  L.makeCompressed();
  return L;
}

//============================================================================================

std::tuple<SparseMatrixF, std::vector<int32_t>> HSmoothMain::graphLaplacian(const TriMesh& tri)
{
  std::vector<int32_t> nUnique;
  for(Eigen::Index i = 0; i < tri.rows(); i++)
  {
    for(Eigen::Index j = 0; j < tri.cols(); j++)
    {
      nUnique.push_back(tri(i, j));
    }
  }

  std::sort(nUnique.begin(), nUnique.end());
  nUnique.erase(std::unique(nUnique.begin(), nUnique.end()), nUnique.end());

  std::vector<TripletF> tripletList;
  tripletList.reserve(nUnique.size() + tri.rows() * tri.cols() * 2);

  std::vector<float> fDiagCount(nUnique.size(), 0.0f);

  TriMesh nSubTri = HSmoothBase::isMember(tri, nUnique);

  DictBase<int32_t>::EdgeDict MyDict;
  for(int32_t i = 0; i < nSubTri.rows(); i++)
  {
    for(int32_t j = 0; j < 3; j++)
    {
      int32_t l = (j + 3) % 3;
      int32_t m = (j + 4) % 3;
      int32_t this_row = nSubTri(i, l);
      int32_t this_col = nSubTri(i, m);
      EdgePair EP = std::make_pair(std::min(this_row, this_col), std::max(this_row, this_col));
      DictBase<int32_t>::EdgeDict::const_iterator got = MyDict.find(EP);
      if(got == MyDict.end())
      {
        // not found yet...
        MyDict.insert({EP, i}); // i.e. the edge, and one of the triangles it belongs to.
        tripletList.push_back(TripletF(this_row, this_col, -1.0));
        tripletList.push_back(TripletF(this_col, this_row, -1.0));
        fDiagCount[this_row] += 1.0f;
        fDiagCount[this_col] += 1.0f;
      }
    }
  }

  for(int32_t i = 0; i < fDiagCount.size(); i++)
  {
    tripletList.push_back(TripletF(i, i, fDiagCount[i]));
  }

  SparseMatrixF MLap(nUnique.size(), nUnique.size());
  MLap.setFromTriplets(tripletList.begin(), tripletList.end());
  MLap.makeCompressed();
  return std::make_tuple(MLap, nUnique);
}

//============================================================================================

MeshNode HSmoothMain::smooth(const MeshNode& nodes, Type type, float threshold, uint64_t iterations)
{
  SparseMatrixF L = laplacian2D(nodes.rows(), type);
  std::vector<int32_t> vidx;

  if(type == Type::Serial)
  {
    vidx = std::vector<int32_t>{0, static_cast<int32_t>(L.cols() - 1)};
  }
  else
  {
    vidx = std::vector<int32_t>{};
  }

  MatIndex nFixed = HSmoothBase::getIndex(vidx);
  return smooth(nodes, nFixed, L, threshold, iterations);
}

//============================================================================================

MeshNode HSmoothMain::smooth(const MeshNode& nodes, const MatIndex& nFixed, const SparseMatrixF& GL, float threshold, uint64_t iterations)
{
  MatIndex nMobile = HSmoothBase::getComplement(nFixed, GL.cols());
  if(nMobile.size() == 0)
  {
    return nodes;
  }

  SparseMatrixF Data = nodes.sparseView();

  std::tuple<SparseMatrixF, SparseMatrixF> dbvp = getDirichletBVP(GL, Data, nFixed, nMobile);
  SparseMatrixF GLRed = std::get<0>(dbvp);
  SparseMatrixF fConst = std::get<1>(dbvp);

  std::tuple<SparseMatrixF, SparseMatrixF> pieces = analyzeLaplacian(GL);
  SparseMatrixF D = std::get<0>(pieces);
  SparseMatrixF A = std::get<1>(pieces);

  SparseMatrixF AyIn = A * Data;
  Eigen::MatrixXf mtemp = Eigen::MatrixXf::Zero(nMobile.size(), nMobile.size());
  mtemp.setIdentity();
  SparseMatrixF fSmallEye = mtemp.sparseView();
  fSmallEye.makeCompressed();
  SparseMatrixF yMobile;
  slice::slice(Data, nMobile, 1, yMobile);
  SparseMatrixF LTL = SparseMatrixF(GLRed.transpose() * GLRed);
  SparseMatrixF LTK = SparseMatrixF(GLRed.transpose() * fConst); // casting as SparseMatrixF to make column-major
  SparseMatrixF yOut = Data;

  Smoother smth;

  float fEps = 0.5f;
  float fStep = fEps / 2.0f;
  uint64_t nCount = 1;

  float fobj1, fobj2, fslope;
  fobj1 = getObjFn(smth, fEps, fSmallEye, LTL, LTK, Data, nMobile, yMobile, D, AyIn, yOut); // only the last parameter is actually modified
  fobj2 = getObjFn(smth, fEps + threshold, fSmallEye, LTL, LTK, Data, nMobile, yMobile, D, AyIn, yOut);
  fslope = (fobj2 - fobj1) / threshold;

  while(fabs(fslope) < threshold && nCount < iterations)
  {
    if(fStep > 0.0f)
    {
      fEps -= fStep;
    }
    else
    {
      fEps += fStep;
    }

    fEps /= 2.0f;
    fobj1 = getObjFn(smth, fEps, fSmallEye, LTL, LTK, Data, nMobile, yMobile, D, AyIn, yOut); // only the last parameter is actually modified
    fobj2 = getObjFn(smth, fEps + threshold, fSmallEye, LTL, LTK, Data, nMobile, yMobile, D, AyIn, yOut);
    fslope = (fobj2 - fobj1) / threshold;
    nCount++;
  }

  return Eigen::MatrixXf(yOut);
}

//============================================================================================

std::tuple<SparseMatrixF, SparseMatrixF> HSmoothMain::getDirichletBVP(const SparseMatrixF& GL, const SparseMatrixF& y, const MatIndex& nFixed, const MatIndex& nMobile)
{
  MatIndex nAll = HSmoothBase::matUnion(nFixed, nMobile);
  std::vector<int32_t> v;
  for(int32_t i = 0; i < y.cols(); i++)
  {
    v.push_back(i);
  }
  MatIndex dims = HSmoothBase::getIndex(v);

  SparseMatrixF GLRed(nMobile.size(), nMobile.size());
  SparseMatrixF sm1(nAll.size(), nFixed.size());
  SparseMatrixF sm2(nFixed.size(), dims.size());
  SparseMatrixF sm3(nAll.size(), dims.size());
  SparseMatrixF fConst(nMobile.size(), dims.size());

  slice::slice(GL, nMobile, nMobile, GLRed); // thank goodness for igl::slice!
  slice::slice(GL, nAll, nFixed, sm1);
  slice::slice(y, nFixed, dims, sm2);
  sm3 = sm1 * sm2;
  slice::slice(sm3, nMobile, dims, fConst);

  return std::make_tuple(GLRed, fConst);
}

//============================================================================================

std::tuple<SparseMatrixF, SparseMatrixF> HSmoothMain::analyzeLaplacian(const SparseMatrixF& GL)
{
  // NOTE: GL should be column-major

  SparseMatrixF D(GL.rows(), GL.cols());
  SparseMatrixF A(GL.rows(), GL.cols());
  std::vector<TripletF> DT;
  std::vector<TripletF> AT;

  for(int32_t k = 0; k < GL.outerSize(); ++k)
  {
    for(SparseMatrixF::InnerIterator iter(GL, k); iter; ++iter)
    {
      if(iter.value() < -0.5)
      {
        AT.push_back(TripletF(iter.row(), iter.col(), iter.value()));
      }
      else if(iter.value() > 0.5)
      {
        DT.push_back(TripletF(iter.row(), iter.col(), iter.value()));
      }
    }
  }
  D.setFromTriplets(DT.cbegin(), DT.cend());
  A.setFromTriplets(AT.cbegin(), AT.cend());

  D.makeCompressed();
  A.makeCompressed();

  return std::make_tuple(D, A);
}

//============================================================================================

float HSmoothMain::getObjFn(Smoother& smth, float feps, const SparseMatrixF& fSmallEye, const SparseMatrixF& LTL, const SparseMatrixF& LTK, const SparseMatrixF& Data, const MatIndex& nMobile,
                            const SparseMatrixF& yMobile, const SparseMatrixF& D, const SparseMatrixF& AyIn, SparseMatrixF& yOut)
{
  SparseMatrixF A = (1.0 - feps) * fSmallEye + feps * LTL;
  SparseMatrixF b = (1.0 - feps) * yMobile - feps * LTK;

  smth.compute(A); // smth in general changes with each function call
  SparseMatrixF ySmooth = smth.solve(b);

  HSmoothBase::merge(ySmooth, yOut, nMobile);

  Eigen::ArrayXXf yDeltaD = Eigen::MatrixXf(D * yOut + AyIn).array();

  return (yDeltaD * yDeltaD).sum(); // more efficient to calculate trace this way
}
