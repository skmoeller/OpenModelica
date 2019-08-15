/*
 * This file is part of OpenModelica.
 *
 * Copyright (c) 1998-2014, Open Source Modelica Consortium (OSMC),
 * c/o Linköpings universitet, Department of Computer and Information Science,
 * SE-58183 Linköping, Sweden.
 *
 * All rights reserved.
 *
 * THIS PROGRAM IS PROVIDED UNDER THE TERMS OF GPL VERSION 3 LICENSE OR
 * THIS OSMC PUBLIC LICENSE (OSMC-PL) VERSION 1.2.
 * ANY USE, REPRODUCTION OR DISTRIBUTION OF THIS PROGRAM CONSTITUTES
 * RECIPIENT'S ACCEPTANCE OF THE OSMC PUBLIC LICENSE OR THE GPL VERSION 3,
 * ACCORDING TO RECIPIENTS CHOICE.
 *
 * The OpenModelica software and the Open Source Modelica
 * Consortium (OSMC) Public License (OSMC-PL) are obtained
 * from OSMC, either from the above address,
 * from the URLs: http://www.ida.liu.se/projects/OpenModelica or
 * http://www.openmodelica.org, and in the OpenModelica distribution.
 * GNU version 3 is obtained from: http://www.gnu.org/copyleft/gpl.html.
 *
 * This program is distributed WITHOUT ANY WARRANTY; without
 * even the implied warranty of  MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE, EXCEPT AS EXPRESSLY SET FORTH
 * IN THE BY RECIPIENT SELECTED SUBSIDIARY LICENSE CONDITIONS OF OSMC-PL.
 *
 * See the full OSMC Public License conditions for more details.
 *
 */

 encapsulated package SymbolicHessian
 " file:        SymbolicHessian.mo
   package:     SymbolicHessian
   description: This package contains stuff to calculate the symbolic Hessian Matrix with the help of the symbolic Jacobian."


 public import Absyn;
 public import BackendDAE;
 public import DAE;
 public import FCore;
 public import FGraph;

 protected
 import Array;
 import BackendDAEOptimize;
 import BackendDAETransform;
 import BackendDAEUtil;
 import BackendDump;
 import BackendEquation;
 import BackendVariable;
 import BackendVarTransform;
 import BaseHashSet;
 import Ceval;
 import ClockIndexes;
 import Config;
 import ComponentReference;
 import Debug;
 import Differentiate;
 import DynamicOptimization;
 import SymbolicJacobian;
 import ElementSource;
 import ExecStat.execStat;
 import ExpandableArray;
 import Expression;
 import ExpressionDump;
 import ExpressionSimplify;
 import Error;
 import Flags;
 import GC;
 import Global;
 import Graph;
 import HashSet;
 import IndexReduction;
 import List;
 import System;
 import Util;
 import Values;
 import ValuesUtil;

 // =============================================================================
 // section for postOptModule >>symbolicJacobian<<
 //
 // Detects the sparse pattern of the ODE system and calculates also the symbolic
 // Jacobian if flag "--generateSymbolicJacobian" is enabled.
 // =============================================================================
