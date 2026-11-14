# License notices

The notices below apply to their respective components. The combined package
uses [GPL-3.0-or-later](../LICENSE). Full copies of
[GPL version 2](GPL-2.0-or-later.txt) and
[LGPL version 2.1](LGPL-2.1-only.txt) accompany the component notices.

## Tessella independent material

This grant applies to independently authored Tessella material used separately;
it does not relicense third-party adaptations.

```text
MIT License

Copyright (c) 2026 Jake Liu

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

## Gmsh

```text
Gmsh is provided under the terms of the GNU General Public License
(GPL), Version 2 or later, with the following exception:

  The copyright holders of Gmsh give you permission to combine Gmsh
  with code included in the standard release of Netgen (from Joachim
  Sch"oberl), METIS (from George Karypis at the University of
  Minnesota), OpenCASCADE (from Open CASCADE S.A.S) and ParaView
  (from Kitware, Inc.) under their respective licenses. You may copy
  and distribute such a system following the terms of the GNU GPL for
  Gmsh and the licenses of the other code concerned, provided that
  you include the source code of that other code when and as the GNU
  GPL requires distribution of source code.

  Note that people who make modified versions of Gmsh are not
  obligated to grant this special exception for their modified
  versions; it is their choice whether to do so. The GNU General
  Public License gives permission to release a modified version
  without this exception; this exception also makes it possible to
  release a modified version which carries forward this exception.

End of exception.
```

The full GPL version 2 text is in [GPL-2.0-or-later.txt](GPL-2.0-or-later.txt).

## HXT

```text
Hxt is provided under the terms of the GNU General Public License
(GPL), Version 2 or later, with the following exception:

  The copyright holders of Hxt give you permission to combine Hxt
  with code included in the standard release of Gmsh (from C. Geuzaine
  and J.-F Remacle) under the Gmsh license. You may copy and distribute
  such a system following the terms of the Gmsh license.

  Note that people who make modified versions of Hxt are not
  obligated to grant this special exception for their modified
  versions; it is their choice whether to do so. The GNU General
  Public License gives permission to release a modified version
  without this exception; this exception also makes it possible to
  release a modified version which carries forward this exception.

End of exception.
```

The full GPL version 2 text is in [GPL-2.0-or-later.txt](GPL-2.0-or-later.txt).

## Tetgen/BR

```text
From Gmsh 4.15.2 src/mesh/tetgenBR.cxx and tetgenBR.h (upstream preamble):

// This file is part of Tetgen/BR, the Boundary Recovery code of TetGen
//
// Copyright 2015-2016 Weierstrass Institute for Applied Analysis and
// Stochastics
//
// This file is relicensed under the terms of LICENSE.txt for use in Gmsh thanks
// to a Software License Agreement between Weierstrass Institute for Applied
// Analysis and Stochastics and GMESH SPRL.

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TetGen                                                                    //
//                                                                           //
// A Quality Tetrahedral Mesh Generator and A 3D Delaunay Triangulator       //
//                                                                           //
// Version 1.5                                                               //
// May 31, 2014                                                              //
//                                                                           //
// Copyright (C) 2002--2016                                                  //
//                                                                           //
// TetGen is freely available through the website: http://www.tetgen.org.    //
//   It may be copied, modified, and redistributed for non-commercial use.   //
//   Please consult the file LICENSE for the detailed copyright notices.     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

```

Tessella's adaptation uses the Gmsh distribution of Tetgen/BR under the
relicensing notice above. The earlier TetGen banner is retained as upstream
attribution. See the [Gmsh grant and exception](#gmsh) and
[full GPL version 2](GPL-2.0-or-later.txt). Julia adaptation by Tessella
contributors, 2026.

## NetworkX

```text
NetworkX is distributed with the 3-clause BSD license.

::

   Copyright (c) 2004-2025, NetworkX Developers
   Aric Hagberg <hagberg@lanl.gov>
   Dan Schult <dschult@colgate.edu>
   Pieter Swart <swart@lanl.gov>
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are
   met:

     * Redistributions of source code must retain the above copyright
       notice, this list of conditions and the following disclaimer.

     * Redistributions in binary form must reproduce the above
       copyright notice, this list of conditions and the following
       disclaimer in the documentation and/or other materials provided
       with the distribution.

     * Neither the name of the NetworkX Developers nor the names of its
       contributors may be used to endorse or promote products derived
       from this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
   OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

## HP/SGI STL

```text
# Adapted from the original HP/SGI STL stl_algo.h and stl_heap.h,
# https://github.com/karottc/sgi-stl/blob/master/stl_algo.h
# https://github.com/karottc/sgi-stl/blob/master/stl_heap.h
# Full notices from both permissive source headers are retained below.

#
# Copyright (c) 1994
# Hewlett-Packard Company
#
# Permission to use, copy, modify, distribute and sell this software
# and its documentation for any purpose is hereby granted without fee,
# provided that the above copyright notice appear in all copies and
# that both that copyright notice and this permission notice appear
# in supporting documentation.  Hewlett-Packard Company makes no
# representations about the suitability of this software for any
# purpose.  It is provided "as is" without express or implied warranty.
#
#
# Copyright (c) 1996
# Silicon Graphics Computer Systems, Inc.
#
# Permission to use, copy, modify, distribute and sell this software
# and its documentation for any purpose is hereby granted without fee,
# provided that the above copyright notice appear in all copies and
# that both that copyright notice and this permission notice appear
# in supporting documentation.  Silicon Graphics makes no
# representations about the suitability of this software for any
# purpose.  It is provided "as is" without express or implied warranty.

#
# Copyright (c) 1994
# Hewlett-Packard Company
#
# Permission to use, copy, modify, distribute and sell this software
# and its documentation for any purpose is hereby granted without fee,
# provided that the above copyright notice appear in all copies and
# that both that copyright notice and this permission notice appear
# in supporting documentation.  Hewlett-Packard Company makes no
# representations about the suitability of this software for any
# purpose.  It is provided "as is" without express or implied warranty.
#
# Copyright (c) 1997
# Silicon Graphics Computer Systems, Inc.
#
# Permission to use, copy, modify, distribute and sell this software
# and its documentation for any purpose is hereby granted without fee,
# provided that the above copyright notice appear in all copies and
# that both that copyright notice and this permission notice appear
# in supporting documentation.  Silicon Graphics makes no
# representations about the suitability of this software for any
# purpose.  It is provided "as is" without express or implied warranty.
```

## MinGW-w64

```text
Windows native exp and trig adaptations, modified 2026-10-09.
Upstream: https://github.com/mingw-w64/mingw-w64/tree/2bcf06abc9ec5f8f73c743896b815a7a7654f9ae

exp.def.h (ZPL2.1):
/*
 This Software is provided under the Zope Public License (ZPL) Version 2.1.

 Copyright (c) 2009, 2010 by the mingw-w64 project

 See the AUTHORS file for the list of contributors to the mingw-w64 project.

 This license has been certified as open source. It has also been designated
 as GPL compatible by the Free Software Foundation (FSF).

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

   1. Redistributions in source code must retain the accompanying copyright
      notice, this list of conditions, and the following disclaimer.
   2. Redistributions in binary form must reproduce the accompanying
      copyright notice, this list of conditions, and the following disclaimer
      in the documentation and/or other materials provided with the
      distribution.
   3. Names of the copyright holders must not be used to endorse or promote
      products derived from this software without prior written permission
      from the copyright holders.
   4. The right to distribute this software or to use it for any purpose does
      not give you the right to use Servicemarks (sm) or Trademarks (tm) of
      the copyright holders.  Use of them is covered by separate agreement
      with the copyright holders.
   5. If any files are modified, you must cause the modified files to carry
      prominent notices stating that you changed the files and the date of
      any change.

 Disclaimer

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY EXPRESSED
 OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
 EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT,
 INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

sinl_internal.S/cosl_internal.S (public domain):
/**
 * This file has no copyright assigned and is placed in the Public Domain.
 * This file is part of the mingw-w64 runtime package.
 * No warranty is given; refer to the file DISCLAIMER.PD within this package.
 */
```

## GCC Runtime Library Exception

```text
GCC RUNTIME LIBRARY EXCEPTION

Version 3.1, 31 March 2009

Copyright (C) 2009 Free Software Foundation, Inc. <http://fsf.org/>

Everyone is permitted to copy and distribute verbatim copies of this
license document, but changing it is not allowed.

This GCC Runtime Library Exception ("Exception") is an additional
permission under section 7 of the GNU General Public License, version
3 ("GPLv3"). It applies to a given file (the "Runtime Library") that
bears a notice placed by the copyright holder of the file stating that
the file is governed by GPLv3 along with this Exception.

When you use GCC to compile a program, GCC may combine portions of
certain GCC header files and runtime libraries with the compiled
program. The purpose of this Exception is to allow compilation of
non-GPL (including proprietary) programs to use, in this way, the
header files and runtime libraries covered by this Exception.

0. Definitions.

A file is an "Independent Module" if it either requires the Runtime
Library for execution after a Compilation Process, or makes use of an
interface provided by the Runtime Library, but is not otherwise based
on the Runtime Library.

"GCC" means a version of the GNU Compiler Collection, with or without
modifications, governed by version 3 (or a specified later version) of
the GNU General Public License (GPL) with the option of using any
subsequent versions published by the FSF.

"GPL-compatible Software" is software whose conditions of propagation,
modification and use would permit combination with GCC in accord with
the license of GCC.

"Target Code" refers to output from any compiler for a real or virtual
target processor architecture, in executable form or suitable for
input to an assembler, loader, linker and/or execution
phase. Notwithstanding that, Target Code does not include data in any
format that is used as a compiler intermediate representation, or used
for producing a compiler intermediate representation.

The "Compilation Process" transforms code entirely represented in
non-intermediate languages designed for human-written code, and/or in
Java Virtual Machine byte code, into Target Code. Thus, for example,
use of source code generators and preprocessors need not be considered
part of the Compilation Process, since the Compilation Process can be
understood as starting with the output of the generators or
preprocessors.

A Compilation Process is "Eligible" if it is done using GCC, alone or
with other GPL-compatible software, or if it is done without using any
work based on GCC. For example, using non-GPL-compatible Software to
optimize any GCC intermediate representations would not qualify as an
Eligible Compilation Process.

1. Grant of Additional Permission.

You have permission to propagate a work of Target Code formed by
combining the Runtime Library with Independent Modules, even if such
propagation would otherwise violate the terms of GPLv3, provided that
all Target Code was generated by Eligible Compilation Processes. You
may then convey such a combination under terms of your choice,
consistent with the licensing of the Independent Modules.

2. No Weakening of GCC Copyleft.

The availability of this Exception does not imply any general
presumption that third-party software is unaffected by the copyleft
requirements of the license of GCC.
```

## OCCT exception

```text
Open CASCADE exception (version 1.0) to GNU LGPL version 2.1.

The object code (i.e. not a source) form of a "work that uses the Library"
can incorporate material from a header file that is part of the Library.
As a special exception to the GNU Lesser General Public License version 2.1,
you may distribute such object code incorporating material from header files
provided with the Open CASCADE Technology libraries (including code of CDL
generic classes) under terms of your choice, provided that you give
prominent notice in supporting documentation to this code that it makes use
of or is based on facilities provided by the Open CASCADE Technology software.
```
