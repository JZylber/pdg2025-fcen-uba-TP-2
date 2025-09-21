//------------------------------------------------------------------------
//  Copyright (C) Gabriel Taubin
//  Time-stamp: <2025-08-05 16:34:29 taubin>
//------------------------------------------------------------------------
//
// PolygonMesh.cpp
//
// Software developed for the course
// Digital Geometry Processing
// Copyright (c) 2025, Gabriel Taubin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
//     * Redistributions of source code must retain the above
//       copyright notice, this list of conditions and the following
//       disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials
//       provided with the distribution.
//     * Neither the name of the Brown University nor the names of its
//       contributors may be used to endorse or promote products
//       derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL GABRIEL
// TAUBIN BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
// OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
// USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.

#include <iostream>
#include "PolygonMesh.hpp"
#include "Partition.hpp"

PolygonMesh::PolygonMesh(const int nVertices, const vector<int> &coordIndex) : HalfEdges(nVertices, coordIndex),
                                                                               _nPartsVertex(),
                                                                               _isBoundaryVertex()
{
  int nV = getNumberOfVertices();
  int nE = getNumberOfEdges(); // Edges method
  // int nF = getNumberOfFaces();
  int nC = getNumberOfCorners();

  // 1) classify the vertices as boundary or internal
  _isBoundaryVertex.resize(nV, false);
  for (int iE = 0; iE < nE; iE++)
  {
    if (isBoundaryEdge(iE))
    {
      int iV0 = getVertex0(iE);
      int iV1 = getVertex1(iE);
      _isBoundaryVertex[iV0] = true;
      _isBoundaryVertex[iV1] = true;
    }
  }
  // - for edge boundary iE label its two end vertices as boundary

  // 2) create a partition of the corners in the stack
  Partition partition(nC);
  // 3) for each regular edge
  //    - get the two half edges incident to the edge
  //    - join the two pairs of corresponding corners accross the edge
  //    - you need to take into account the relative orientation of
  //      the two incident half edges
  for (int iE = 0; iE < nE; iE++)
  {
    if (!isRegularEdge(iE))
      continue;
    int nHE = getNumberOfEdgeHalfEdges(iE);
    if (nHE != 2)
      continue; // should not happen
    int iC0 = getEdgeHalfEdge(iE, 0);
    int iC1 = getEdgeHalfEdge(iE, 1);
    int iV00 = getSrc(iC0);
    int iV01 = getDst(iC0);
    int iV10 = getSrc(iC1);
    int iV11 = getDst(iC1);
    if (iV00 == iV11 && iV01 == iV10)
    {
      // consistently oriented
      partition.join(iC0, iC1);
      partition.join((iC0 + 1) % nC, (iC1 + 1) % nC);
    }
    else
    {
      // opposite orientation
      partition.join(iC0, (iC1 + 1) % nC);
      partition.join((iC0 + 1) % nC, iC1);
    }
  }
  // consistently oriented
  /* \                  / */
  /*  \ iC01 <-- iC00  /  */
  /*   X ---- iE ---- X   */
  /*  / iC10 --> iC11  \  */
  /* /                  \ */

  // oposite orientation
  /* \                  / */
  /*  \ iC01 --> iC00  /  */
  /*   X ---- iE ---- X   */
  /*  / iC10 --> iC11  \  */
  /* /                  \ */

  // a decision has to be made about inconsistently oriented faces
  // incident to the same edge, as well as how to deal with singular
  // edges; for the moment let's assume that the mesh does not have
  // singular edges, and that pairs of corners corresponding to the
  // same vertex across inconsistently oriented faces will be joined

  // note that the partition will end up with the corner separators as
  // singletons, but it doesn't matter for the last step, and
  // the partition will be deleted upon return

  // 4) count number of parts per vertex
  //    - initialize _nPartsVertex array to 0's
  //    - for each corner iC which is a representative of its subset,
  //    - get the corresponding vertex index iV and increment _nPartsVertex[iV]
  //    - note that all the corners in each subset share a common
  //      vertex index, but multiple subsets may correspond to the
  //      same vertex index, indicating that the vertex is singular
  _nPartsVertex.resize(nV, 0);
  for (int iC = 0; iC < nC; iC++)
  {
    if (partition.find(iC) == iC)
    {
      int iV = _coordIndex[iC];
      if (iV >= 0)
        _nPartsVertex[iV]++;
    }
  }
}

int PolygonMesh::getNumberOfFaces() const
{
  int nC = getNumberOfCorners();
  int nF = 0;
  for (int iC = 0; iC < nC; iC++)
  {
    if (_coordIndex[iC] < 0)
      nF++;
  }
  return nF;
}

int PolygonMesh::getNumberOfEdgeFaces(const int iE) const
{
  return getNumberOfEdgeHalfEdges(iE);
}

int PolygonMesh::getEdgeFace(const int iE, const int j) const
{
  if (!isValidEdgeHalfEdge(iE, j))
    return -1;
  int iC = getEdgeHalfEdge(iE, j);
  return getFace(iC);
}

bool PolygonMesh::isEdgeFace(const int iE, const int iF) const
{
  if (!isValidIe(iE))
    return false;
  int nF = getNumberOfFaces();
  if (!(0 <= iF && iF < nF))
    return false;
  int nHE = getNumberOfEdgeHalfEdges(iE);
  for (int j = 0; j < nHE; j++)
  {
    int iC = getEdgeHalfEdge(iE, j);
    if (getFace(iC) == iF)
      return true;
  }
}

// classification of edges

bool PolygonMesh::isBoundaryEdge(const int iE) const
{
  if (!isValidIe(iE))
    return false;
  return (getNumberOfEdgeHalfEdges(iE) == 1);
}

bool PolygonMesh::isRegularEdge(const int iE) const
{
  if (!isValidIe(iE))
    return false;
  return (getNumberOfEdgeHalfEdges(iE) == 2);
}

bool PolygonMesh::isSingularEdge(const int iE) const
{
  if (!isValidIe(iE))
    return false;
  return (getNumberOfEdgeHalfEdges(iE) > 2);
}

// classification of vertices

bool PolygonMesh::isBoundaryVertex(const int iV) const
{
  int nV = getNumberOfVertices();
  return (0 <= iV && iV < nV) ? _isBoundaryVertex[iV] : false;
}

bool PolygonMesh::isSingularVertex(const int iV) const
{
  int nV = getNumberOfVertices();
  return (0 <= iV && iV < nV && _nPartsVertex[iV] > 1);
}

// properties of the whole mesh

bool PolygonMesh::isRegular() const
{
  // TODO
  return false;
}

bool PolygonMesh::hasBoundary() const
{
  // TODO
  return false;
}
