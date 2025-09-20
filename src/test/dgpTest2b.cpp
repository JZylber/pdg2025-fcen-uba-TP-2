//------------------------------------------------------------------------
//  Copyright (C) Gabriel Taubin
//  Time-stamp: <2025-08-05 16:38:52 taubin>
//------------------------------------------------------------------------
//
// dgpTest2b.cpp
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

#include <string>
#include <iostream>

using namespace std;

#include <wrl/SceneGraphTraversal.hpp>
#include <io/AppLoader.hpp>
#include <io/AppSaver.hpp>
#include <io/LoaderPly.hpp>
#include <io/LoaderStl.hpp>
#include <io/LoaderWrl.hpp>
#include <io/SaverPly.hpp>
#include <io/SaverStl.hpp>
#include <io/SaverWrl.hpp>
#include "dgpPrt.hpp"

class Data
{
public:
  bool _debug;
  bool _binaryOutput;
  bool _removeNormal;
  bool _removeColor;
  bool _removeTexCoord;
  string _inFile;
  string _outFile;

public:
  Data() : _debug(false),
           _binaryOutput(false),
           _removeNormal(false),
           _removeColor(false),
           _removeTexCoord(false),
           _inFile(""),
           _outFile("")
  {
  }
};

void options(Data &D)
{
  cout << "   -d|-debug               [" << tv(D._debug) << "]" << endl;
  cout << "   -b|-binaryOutput        [" << tv(D._binaryOutput) << "]" << endl;
  cout << "   -r|-removeProperties    [" << tv(D._removeNormal) << "]" << endl;
  cout << "  -rn|-removeNormal        [" << tv(D._removeNormal) << "]" << endl;
  cout << "  -rc|-removeColor         [" << tv(D._removeColor) << "]" << endl;
  cout << "  -rt|-removeTexCoord      [" << tv(D._removeTexCoord) << "]" << endl;
}

void usage(Data &D)
{
  cout << "USAGE: dgpTest2b [options] inFile outFile" << endl;
  cout << "   -h|-help" << endl;
  options(D);
  cout << endl;
  exit(0);
}

void error(const char *msg)
{
  cout << "ERROR: dgpTest2b | " << ((msg) ? msg : "") << endl;
  exit(0);
}

//////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{

  Data D;

  if (argc == 1)
    usage(D);

  for (int i = 1; i < argc; i++)
  {
    if (string(argv[i]) == "-h" || string(argv[i]) == "-help")
    {
      usage(D);
    }
    else if (string(argv[i]) == "-d" || string(argv[i]) == "-debug")
    {
      D._debug = !D._debug;
    }
    else if (string(argv[i]) == "-b" || string(argv[i]) == "-binaryOutput")
    {
      D._binaryOutput = !D._binaryOutput;
    }
    else if (string(argv[i]) == "-r" || string(argv[i]) == "-removeProperties")
    {
      D._removeNormal = !D._removeNormal;
      D._removeColor = !D._removeColor;
      D._removeTexCoord = !D._removeTexCoord;
    }
    else if (string(argv[i]) == "-rn" || string(argv[i]) == "-removeNormal")
    {
      D._removeNormal = !D._removeNormal;
    }
    else if (string(argv[i]) == "-rc" || string(argv[i]) == "-removeColor")
    {
      D._removeColor = !D._removeColor;
    }
    else if (string(argv[i]) == "-rt" || string(argv[i]) == "-removeTexCoord")
    {
      D._removeTexCoord = !D._removeTexCoord;
    }
    else if (string(argv[i])[0] == '-')
    {
      error("unknown option");
    }
    else if (D._inFile == "")
    {
      D._inFile = string(argv[i]);
    }
    else if (D._outFile == "")
    {
      D._outFile = string(argv[i]);
    }
    else
    {
      error("too many command line arguments");
    }
  }

  if (D._inFile == "")
    error("no inFile");
  if (D._outFile == "")
    error("no outFile");

  if (D._debug)
  {
    cout << "dgpTest2b {" << endl;
    cout << endl;
    options(D);
    cout << endl;
    cout << "  inFile  = " << D._inFile << endl;
    cout << "  outFile = " << D._outFile << endl;
    cout << endl;
    fflush(stderr);
  }

  bool success;

  //////////////////////////////////////////////////////////////////////
  // create loader and saver factories
  AppLoader loaderFactory;
  AppSaver saverFactory;

  // register input file loaders
  LoaderPly *plyLoader = new LoaderPly();
  loaderFactory.registerLoader(plyLoader);
  LoaderStl *stlLoader = new LoaderStl();
  loaderFactory.registerLoader(stlLoader);
  LoaderWrl *wrlLoader = new LoaderWrl();
  loaderFactory.registerLoader(wrlLoader);

  // register output file savers
  SaverPly *plySaver = new SaverPly();
  saverFactory.registerSaver(plySaver);
  SaverStl *stlSaver = new SaverStl();
  saverFactory.registerSaver(stlSaver);
  SaverWrl *wrlSaver = new SaverWrl();
  saverFactory.registerSaver(wrlSaver);

  SaverStl::FileType stlFt =
      (D._binaryOutput) ? SaverStl::FileType::BINARY : SaverStl::FileType::ASCII;
  stlSaver->setFileType(stlFt);

  Ply::DataType plyDt =
      (D._binaryOutput) ? Ply::DataType::BINARY_LITTLE_ENDIAN : Ply::DataType::ASCII;
  SaverPly::setDefaultDataType(plyDt);

  if (D._debug)
  {
    SaverPly::setOstream(&cout);
    SaverPly::setIndent("    ");
  }

  //////////////////////////////////////////////////////////////////////
  // read ScheneGraph

  SceneGraph wrl; // create empty scene graph

  if (D._debug)
  {
    cout << "  loading inFile {" << endl;
  }

  success = loaderFactory.load(D._inFile.c_str(), wrl);

  if (D._debug)
  {
    cout << "    success        = " << tv(success) << endl;
    cout << "  } loading inFile" << endl;
    cout << endl;
  }

  if (success == false)
    return -1;

  //////////////////////////////////////////////////////////////////////
  // process

  if (D._removeNormal || D._removeColor || D._removeTexCoord)
  {
    if (D._debug)
      cout << "  processing {" << endl;

    // traverse scene graph looking for IndexedFaceSet nodes
    Node *node;
    SceneGraphTraversal sgt(wrl);
    for (int iIfs = 0; (node = sgt.next()) != (Node *)0; iIfs++)
    {
      Shape *shape = dynamic_cast<Shape *>(node);
      if (shape == (Shape *)0)
        continue;
      const string &shapeName = shape->getName();

      IndexedFaceSet *ifs =
          dynamic_cast<IndexedFaceSet *>(shape->getGeometry());
      if (ifs == (IndexedFaceSet *)0)
        continue;

      if (D._debug)
      {
        cout << "  before processing" << endl;
        printIndexedFaceSetInfo(cout, shapeName, iIfs, *ifs, "    ");
      }

      // TODO ...
      if (D._removeNormal)
      {
        if (D._debug)
          cout << "    removing normals" << endl;
        ifs->getNormal().clear();
        ifs->getNormalIndex().clear();
        ifs->setNormalPerVertex(false);
      }
      if (D._removeColor)
      {
        if (D._debug)
          cout << "    removing colors" << endl;
        ifs->getColor().clear();
        ifs->getColorIndex().clear();
        ifs->setColorPerVertex(false);
      }
      if (D._removeTexCoord)
      {
        if (D._debug)
          cout << "    removing texCoords" << endl;
        ifs->getTexCoord().clear();
        ifs->getTexCoordIndex().clear();
      }

      if (D._debug)
      {
        cout << "  after processing" << endl;
        printIndexedFaceSetInfo(cout, shapeName, iIfs, *ifs, "    ");
      }

      if (D._debug)
        cout << "" << endl;
    }

    if (D._debug)
      cout << "  }" << endl;
  }

  //////////////////////////////////////////////////////////////////////
  // write

  if (D._debug)
  {
    cout << "  saving outFile {" << endl;
  }

  success = saverFactory.save(D._outFile.c_str(), wrl);

  if (D._debug)
  {
    cout << "    success        = " << tv(success) << endl;
    cout << "  } saving outFile" << endl;
    cout << endl;
  }

  //////////////////////////////////////////////////////////////////////

  if (D._debug)
  {
    cout << "} dgpTest2b" << endl;
    fflush(stderr);
  }

  return 0;
}
