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
 * VolumcSolver.cpp:
 * Defines classes and routines declared in HierarchicalSmooth.h
 */

#include "VolumeSolver.h"

#include <sstream>

#include "Base.h"
#include "HierarchicalSmooth.h"
#include "Slice.h"
#include "Triangulation.h"

namespace
{
const HierarchicalSmooth::MatIndex k_One = HierarchicalSmooth::getIndex({0});
const HierarchicalSmooth::MatIndex k_Three = HierarchicalSmooth::getIndex({0, 1, 2});

HierarchicalSmooth::TriMesh sliceMesh(const Eigen::Ref<const HierarchicalSmooth::TriMesh>& mesh, const std::vector<int>& patches)
{
  HierarchicalSmooth::TriMesh triSub;
  slice::slice(mesh, HierarchicalSmooth::getIndex(patches), k_Three, triSub);
  return triSub;
}

void markSectionAsComplete(HierarchicalSmooth::IsSmoothed& status, const HierarchicalSmooth::MatIndex& idx)
{
  for(Eigen::Index i = 0; i < idx.size(); i++)
  {
    status(idx(i)) = true;
  }
}

} // namespace

void HierarchicalSmooth::hierarchicalSmooth(Eigen::Ref<TriMesh> volumeMesh, const Eigen::Ref<const MeshNode>& surfaceNodes, const Eigen::Ref<const FaceLabel>& faceLabels,
                                      const Eigen::Ref<const NodeType>& nodeTypes, Eigen::Ref<MeshNode> smoothedNodes, float threshold, uint64_t iterations, LogCallback logFunction)
{
  DictBase<std::vector<int>>::EdgeDict boundaryDict;

  IsSmoothed status = IsSmoothed(nodeTypes.size());
  for(Eigen::Index i = 0; i < nodeTypes.size(); i++)
  {
    // quad jn points considered already smoothed.
    status(i) = (nodeTypes(i) % 10 == 4);
  }

  smoothedNodes = surfaceNodes;

  float error = (surfaceNodes.row(volumeMesh(0, 0)).array() - surfaceNodes.row(volumeMesh(0, 1)).array()).matrix().norm();
  error = std::min(error, (surfaceNodes.row(volumeMesh(0, 1)).array() - surfaceNodes.row(volumeMesh(0, 2)).array()).matrix().norm());
  error = std::sqrt(3.0f * error * error);
  constexpr float errorThreshold = 2.0f;

  // Filling in boundary dictionary. This also takes care of
  // flipping the direction of a mesh element as computed by
  // Dream.3d so that every element has the same hendedness.
  std::stringstream ss;
  ss.str("");
  ss << "Adjusting Winding of Triangles....";
  logFunction(ss.str());

  for(Eigen::Index i = 0; i < faceLabels.rows(); i++)
  {
    int32_t min = std::min(faceLabels(i, 0), faceLabels(i, 1));
    int32_t max = std::max(faceLabels(i, 0), faceLabels(i, 1));
    if(faceLabels(i, 0) == min)
    {
      // this line ensures consistent handedness for entire surface
      std::swap(volumeMesh(i, 0), volumeMesh(i, 1));
    }
    EdgePair edgePair = std::make_pair(min, max);
    DictBase<std::vector<int32_t>>::EdgeDict::iterator iter = boundaryDict.find(edgePair);
    if(iter == boundaryDict.end())
    {
      // new boundary, new dictionary entry
      std::vector<int32_t> v;
      v.push_back(static_cast<int32_t>(i));
      boundaryDict.insert({edgePair, v});
    }
    else
    {
      // this patch belongs to an already found boundary
      std::vector<int32_t>& vec = iter->second;
      vec.push_back(static_cast<int32_t>(i));
    }
  }

  ss.str("");
  ss << "Looping on the boundary dictionary of " << boundaryDict.size() << " elements....";
  logFunction(ss.str());
  for(auto iter = boundaryDict.cbegin(); iter != boundaryDict.cend(); iter++)
  {
    TriMesh triSub = sliceMesh(volumeMesh, iter->second);
    HierarchicalSmooth::Triangulation triangulation(triSub);

    std::tuple<SparseMatrixF, MatIndex> topology = triangulation.graphLaplacian();
    SparseMatrixF GL = std::get<0>(topology);
    MatIndex nUniq = std::get<1>(topology);

    std::tuple<EdgeList, EdgeList> freeBoundryData = triangulation.freeBoundary();
    EdgeList freeBoundry = std::get<0>(freeBoundryData);
    EdgeList freeBoundrySegments = std::get<1>(freeBoundryData);

    for(size_t i = 0; i < freeBoundrySegments.size(); i++)
    {
      // smooth each free boundary first
      int32_t count = 0;
      int32_t start = std::get<0>(freeBoundrySegments[i]);
      int32_t stop = std::get<1>(freeBoundrySegments[i]);
      for(count = start; count <= stop; count++)
      {
        if(nodeTypes(std::get<0>(freeBoundry[count])) % 10 == 4)
        {
          break;
        }
      }

      std::vector<int32_t> vtemp;
      if(count > stop)
      {
        // no quad jns in this free boundary. smooth without constraints and get out.
        for(count = start; count <= stop; count++)
        {
          vtemp.push_back(std::get<0>(freeBoundry[count]));
        }

        MatIndex thisFreeBoundaryIdx = HierarchicalSmooth::getIndex(vtemp);
        MeshNode thisFreeBoundary;
        slice::slice(smoothedNodes, thisFreeBoundaryIdx, k_Three, thisFreeBoundary);
        MeshNode thisFreeBoundarySmooth = HierarchicalSmooth::smooth(thisFreeBoundary, HierarchicalSmooth::Type::Cyclic, threshold, iterations);
        HierarchicalSmooth::merge(thisFreeBoundarySmooth, smoothedNodes, thisFreeBoundaryIdx);
        markSectionAsComplete(status, thisFreeBoundaryIdx);
      }
      else
      {
        // triple line sections found; smooth separately.
        vtemp.push_back(std::get<0>(freeBoundry[count]));
        int32_t size = 1 + (stop - start);
        for(int32_t j = count + 1; j < 1 + count + size; j++)
        {
          int32_t effective_j = j % size;
          vtemp.push_back(std::get<0>(freeBoundry[effective_j]));
          if(nodeTypes(std::get<0>(freeBoundry[effective_j])) % 10 == 4)
          {
            // reached terminal quad point
            MatIndex thisTripleLineIndex = HierarchicalSmooth::getIndex(vtemp);
            IsSmoothed thisStatus;
            slice::slice(status, thisTripleLineIndex, k_One, thisStatus);
            if(!thisStatus.all())
            {
              MeshNode thisTripleLine, thisTripleLineSmoothed;
              slice::slice(smoothedNodes, thisTripleLineIndex, k_Three, thisTripleLine);
              thisTripleLineSmoothed = HierarchicalSmooth::smooth(thisTripleLine, HierarchicalSmooth::Type::Serial, threshold, iterations);
              HierarchicalSmooth::merge(thisTripleLineSmoothed, smoothedNodes, thisTripleLineIndex); /* HAVEN'T CHECKED FOR BUGS */
              markSectionAsComplete(status, thisTripleLineIndex);
            }
            vtemp.clear();
            vtemp.push_back(std::get<0>(freeBoundry[effective_j]));
          }
        }
      }
    }
    // NOW, smooth entire boundary subject to fixed triple points.
    MeshNode boundaryNode;
    slice::slice(smoothedNodes, nUniq, k_Three, boundaryNode);
    std::vector<int32_t> fixed;
    for(size_t i = 0; i < freeBoundry.size(); i++)
    {
      fixed.push_back(std::get<0>(freeBoundry[i]));
    }
    MatIndex nFixed = HierarchicalSmooth::getIndex(fixed, nUniq);
    MeshNode BoundaryNodeSmooth = HierarchicalSmooth::smooth(boundaryNode, nFixed, GL, threshold, iterations);
    HierarchicalSmooth::merge(BoundaryNodeSmooth, smoothedNodes, nUniq);
    markSectionAsComplete(status, nUniq);
    // Done. Get out.
  }

  Eigen::ArrayX3f tempArray = smoothedNodes.array() - surfaceNodes.array();
  Eigen::ArrayXf normArray = (tempArray * tempArray).rowwise().sum().sqrt() / error;
  if((normArray > errorThreshold).any())
  {
    for(Eigen::Index i = 0; i < status.rows(); i++)
      if(normArray(i, 0) > errorThreshold)
      {
        status(i, 0) = false;
        smoothedNodes.row(i) << surfaceNodes.row(i); // reset to old values.
      }
  }

  if(logFunction)
  {
    std::stringstream ss;
    if(!status.all())
    {
      ss << "WARNING: " << (normArray > errorThreshold).count() << " of " << smoothedNodes.rows() << " nodes not smoothed.";
    }
    else
    {
      ss << "All nodes smoothed";
    }
    logFunction(ss.str());
  }
}
