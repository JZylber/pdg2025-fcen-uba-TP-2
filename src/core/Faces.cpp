//------------------------------------------------------------------------
//  Copyright (C) Gabriel Taubin
//  Time-stamp: <2025-08-04 22:09:56 gtaubin>
//------------------------------------------------------------------------
//
// Faces.cpp
//
// Written by: <Your Name>
//
// Software developed for the course
// Digital Geometry Processing
// Copyright (c) 2025, Gabriel Taubin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Brown University nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL GABRIEL TAUBIN BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <math.h>
#include "Faces.hpp"

Faces::Faces(const int nV, const vector<int> &coordIndex)
{
  _coordIndex = coordIndex;
  int maxVertex = nV - 1;
  size_t currentFace = 0;
  for (size_t i = 0; i < _coordIndex.size(); i++)
  {
    if (currentFace >= _faceStartIndex.size())
    {
      _faceStartIndex.push_back(i);
    }
    if (_coordIndex[i] == -1)
    {
      _coordIndex[i] = -1 * (currentFace + 1); // mark the -1 separator with a negative face index
      currentFace++;
    }
    if (_coordIndex[i] > maxVertex)
    {
      maxVertex = _coordIndex[i];
    }
  }
  _nV = maxVertex + 1;
}

bool Faces::isValidFaceIndex(int iF) const
{
  return (iF >= 0 && iF < getNumberOfFaces());
}

bool Faces::isNonSeparatorVertex(const int iC) const
{
  return iC >= 0 && iC < getNumberOfCorners() && _coordIndex[iC] != -1;
}

int Faces::getNumberOfVertices() const
{
  return _nV;
}

int Faces::getNumberOfFaces() const
{
  return (int)_faceStartIndex.size();
}

int Faces::getNumberOfCorners() const
{
  return (int)_coordIndex.size();
}

int Faces::getFaceSize(const int iF) const
{
  if (!isValidFaceIndex(iF))
    return 0;
  int faceStart = getFaceFirstCorner(iF);
  int faceEnd = (iF + 1 < getNumberOfFaces()) ? getFaceFirstCorner(iF + 1) : (int)_coordIndex.size();
  return faceEnd - faceStart - 1; // exclude the -1 separator
}
int Faces::getFaceFirstCorner(const int iF) const
{
  if (!isValidFaceIndex(iF))
    return -1;
  return _faceStartIndex[iF];
}

int Faces::getFaceVertex(const int iF, const int j) const
{
  if (!isValidFaceIndex(iF))
    return -1;
  int faceSize = getFaceSize(iF);
  if (j >= 0 && j < faceSize)
  {
    int faceStart = getFaceFirstCorner(iF);
    return _coordIndex[faceStart + j];
  }
  return -1;
}

int Faces::getCornerFace(const int iC) const
{
  if (!isNonSeparatorVertex(iC))
    return -1;
  size_t i = (size_t)iC;
  while (_coordIndex[i] >= 0)
    i++;
  int faceNumber = -1 * _coordIndex[i] - 1;
  return faceNumber;
}

int Faces::getNextCorner(const int iC) const
{
  if (!isNonSeparatorVertex(iC))
    return -1;
  int faceNumber = getCornerFace(iC);
  int faceSize = getFaceSize(faceNumber);
  int faceStart = getFaceFirstCorner(faceNumber);
  return (iC + 1 - faceStart) % faceSize + faceStart;
}
