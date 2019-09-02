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
 //
 //
 //
 //
 // =============================================================================

public function generateSymbolicHessian
  "Function to generate the symbolic hessian with respect to the stats of an dynamic optimization problem."
  input BackendDAE.BackendDAE inBackendDAE "Input BackendDAE";
  output BackendDAE.BackendDAE outHessian "second derivates-> this is the hessian"; //Improve it by using special hessian struct
protected
  BackendDAE.SymbolicJacobians linearModelMatrixes;
  BackendDAE.SymbolicJacobian jacA;
  BackendDAE.SymbolicJacobian jacB;
  BackendDAE.SymbolicJacobian jacC;
  BackendDAE.SymbolicJacobian jacD;
  Option<list<DAE.ComponentRef>> lambdas;
algorithm
  outHessian:=SymbolicJacobian.generateSymbolicLinearizationPast(inBackendDAE);
  linearModelMatrixes:=BackendDAEUtil.getSharedSymJacs(outHessian.shared);
  (SOME(jacA),_,_)::(SOME(jacB),_,_)::(SOME(jacC),_,_)::(SOME(jacD),_,_)::linearModelMatrixes:=linearModelMatrixes;
  outHessian:=generateSymbolicHessianA(jacA);
  print("\n\nOutput after linearization\n\n");
  BackendDump.printBackendDAE(outHessian);
end generateSymbolicHessian;

protected function generateSymbolicHessianA
  "Function sets the lagrange factors and multiplies the vector to the jacobian.
   Then it runs the jacobian routine again!"
  input  BackendDAE.SymbolicJacobian A "Symbolic Jacobian Matrix A";
  output BackendDAE.BackendDAE HessA "Symbolic Hessian Matrix for A";
protected
  list<BackendDAE.Var> stats;
  Option<list<DAE.ComponentRef>> lambdas;
algorithm
  (HessA,_,_,_,_,_):=A;
  /*
  // somehow get states
  (_,_,stats,_,_,_):=A;
   lambdas:=if Flags.getConfigBool(Flags.GENERATE_SYMBOLIC_HESSIAN) then SOME(getLambdaList(listLength(stats))) else NONE();

  if isSome(lambdas) then
    A:=multiplyLambdas(lambdas,A);
  end if;

// add up equations to one equation
// differentiate equation wrt all states
*/
end generateSymbolicHessianA;

protected function getLambdaList
  "Function sets the lambdas to the system."
  input Integer lambdaCount "Number of lambdas";
  output list<DAE.ComponentRef> lambdas = {} "List of componentrefs for lambdas";
algorithm
  for i in lambdaCount:-1:1 loop
    lambdas := DAE.CREF_IDENT("$lambda", DAE.T_REAL_DEFAULT, {DAE.INDEX(DAE.ICONST(i))}) ::lambdas;
  end for;
end getLambdaList;

protected function multiplyLambdas
  "Function sets the "
  input Option<list<DAE.ComponentRef>> lambdas;
  input output Option<BackendDAE.SymbolicJacobian> jac;
algorithm
  /*
  get ordered equations from jac
  traverse and multiply each lambda on rhs
  e1.rhs -> e1.rhs * lambda1
  BackendDAEUtil.traverseArrayNoCopyWithUpdate
  */
end multiplyLambdas;

protected function addEquations
  input output BackendDAE.EquationArray eqs;
algorithm

end addEquations;
annotation(__OpenModelica_Interface="backend");
end SymbolicHessian;
