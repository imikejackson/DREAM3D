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

/* Triangulation.cpp -- implementation of routines declared in Triangulation.h */

#include "Triangulation.h"

#include "Base.h"

//===================================================================================

HierarchicalSmooth::Triangulation::Triangulation() = default;

HierarchicalSmooth::Triangulation::Triangulation(const TriMesh& tri)
{
  Mesh = tri;
  std::tuple<EdgeList, EdgeList> GetAllEdges = getEdges(Mesh);
  edge_list = std::get<0>(GetAllEdges);
  free_boundary = std::get<1>(GetAllEdges);
  differentiateFaces();
}

//===================================================================================

HierarchicalSmooth::TriMesh HierarchicalSmooth::Triangulation::connectivityList() const
{
  return Mesh;
}

//===================================================================================

HierarchicalSmooth::EdgeList HierarchicalSmooth::Triangulation::allEdges() const
{
  return edge_list;
}

//===================================================================================

std::tuple<HierarchicalSmooth::EdgeList, HierarchicalSmooth::EdgeList> HierarchicalSmooth::Triangulation::freeBoundary() const
{
  return std::make_tuple(free_boundary, free_boundary_segments);
}

//===================================================================================
void HierarchicalSmooth::Triangulation::differentiateFaces()
{
  int start = std::get<0>(free_boundary[0]);
  std::vector<int32_t> thissec{0};

  for(int32_t n = 1; n < free_boundary.size(); n++)
  {
    if(std::get<1>(free_boundary[n]) == start)
    {
      thissec.push_back(n);
      free_boundary_segments.push_back(std::make_pair(thissec[0], thissec[1]));
      thissec.clear();
    }
    else if(thissec.size() == 0)
    {
      start = std::get<0>(free_boundary[n]);
      thissec.push_back(n);
    }
  }
}

//===================================================================================

/*
 * NOTE:
 * The GetEdges method defined here contains code that is extraneous to its
 * primary function, in order to lay the groundwork for a much faster method
 * that computes the graph Laplacian, than the standalone function defined in
 * HierarchicalSmooth.cpp. This method, internal to the Triangulation object
 * is implemented in the hierarchical smooth algorithm while the standalone
 * is provided for choice.
 */

std::tuple<HierarchicalSmooth::EdgeList, HierarchicalSmooth::EdgeList> HierarchicalSmooth::Triangulation::getEdges(const TriMesh& tri)
{
  for(Eigen::Index i = 0; i < tri.rows(); i++)
  {
    for(Eigen::Index j = 0; j < tri.cols(); j++)
    {
      nUnique.push_back(tri(i, j));
    }
  }

  std::sort(nUnique.begin(), nUnique.end());
  nUnique.erase(std::unique(nUnique.begin(), nUnique.end()), nUnique.end());

  fDiagCount = std::vector<float>(nUnique.size(), 0.0f);
  nSubTri = HierarchicalSmooth::isMember(tri, nUnique);

  for(Eigen::Index i = 0; i < nSubTri.rows(); i++)
  {
    for(int32_t j = 0; j < 3; j++)
    {
      int32_t l = (j + 3) % 3;
      int32_t m = (j + 4) % 3;
      int32_t this_row = nSubTri(i, l);
      int32_t this_col = nSubTri(i, m);
      EdgePair EP = std::make_pair(std::min(this_row, this_col), std::max(this_row, this_col));
      DictBase<EdgeCount>::EdgeDict::iterator got = MyDict.find(EP);
      if(got == MyDict.end())
      {
        // not found yet; this is a new edge.
        EdgeCount EC(nUnique[this_row], nUnique[this_col]);
        MyDict.insert({EP, EC});
        fDiagCount[this_row] += 1.0f;
        fDiagCount[this_col] += 1.0f;
      }
      else
      {
        // this is definitely an interior edge.
        EdgeCount& ec = got->second;
        (ec.ncount)++;
      }
    }
  }

  for(auto iter = MyDict.cbegin(); iter != MyDict.cend(); iter++)
  {
    edge_list.push_back((iter->second).orig_pair);
    if((iter->second).ncount == 1)
    {
      free_boundary.push_back((iter->second).orig_pair);
    }
  }
  return std::make_tuple(edge_list, fastChainLinkSort(free_boundary));
}

//===================================================================================

HierarchicalSmooth::EdgeList HierarchicalSmooth::Triangulation::fastChainLinkSort(const EdgeList& list)
{
  std::unordered_map<int, std::vector<int>> WindingDict;
  for(size_t i = 0; i < list.size(); i++)
  {
    int32_t ltemp = std::get<0>(list[i]);
    int32_t rtemp = std::get<1>(list[i]);
    std::unordered_map<int, std::vector<int>>::iterator got = WindingDict.find(ltemp);
    if(got == WindingDict.end())
    {
      std::vector<int32_t> v{rtemp}; // yet another way to initialize a vector!
      WindingDict.insert({ltemp, v});
    }
    else
    {
      std::vector<int32_t>& vtemp = got->second;
      vtemp.push_back(rtemp);
    }
  }
  EdgeList outList;
  std::unordered_map<int32_t, std::vector<int32_t>>::iterator iter = WindingDict.begin();
  while(WindingDict.size() > 0)
  {
    // decimate the dictionary as chain-linked list is generated.
    int32_t next = (iter->second).back(); // pop backwards; as good a way to fill as any!
    EdgePair ptemp = std::make_pair(iter->first, next);
    outList.push_back(ptemp);
    (iter->second).pop_back();
    if((iter->second).size() == 0)
    {
      WindingDict.erase(iter);
    }
    iter = WindingDict.find(next);
    if(iter == WindingDict.end())
    {
      iter = WindingDict.begin();
    }
  }
  return outList;
}

//===================================================================================

std::tuple<HierarchicalSmooth::SparseMatrixF, HierarchicalSmooth::MatIndex> HierarchicalSmooth::Triangulation::graphLaplacian() const
{
  // most of the work already done in method GetEdges
  std::vector<TripletF> tripletList;
  tripletList.reserve(nUnique.size() + 2 * Mesh.rows() * Mesh.cols());
  for(auto iter = MyDict.cbegin(); iter != MyDict.cend(); iter++)
  {
    int32_t l = std::get<0>(iter->first);
    int32_t m = std::get<1>(iter->first);
    tripletList.push_back(TripletF(l, m, -1.0f));
    tripletList.push_back(TripletF(m, l, -1.0f));
  }
  for(int32_t i = 0; i < fDiagCount.size(); i++)
  {
    tripletList.push_back(TripletF(i, i, fDiagCount[i]));
  }

  SparseMatrixF GL = SparseMatrixF(nUnique.size(), nUnique.size());
  GL.setFromTriplets(tripletList.begin(), tripletList.end());
  GL.makeCompressed();

  return std::make_tuple(GL, HierarchicalSmooth::getIndex(nUnique));
}

//===================================================================================
