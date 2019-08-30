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
  output BackendDAE.BackendDAE outHessian "second derivates-> this is the hessian";
protected
  BackendDAE.SymbolicJacobians linearModelMatrixes;
  BackendDAE.SymbolicJacobian JacA;
  Option<list<DAE.ComponentRef>> lambdas;
algorithm
  outHessian:=SymbolicJacobian.generateSymbolicLinearizationPast(inBackendDAE);
  linearModelMatrixes:=BackendDAEUtil.getSharedSymJacs(outHessian.shared);
  (SOME(JacA),_,_)::linearModelMatrixes:=linearModelMatrixes;
  (outHessian,_,_,_,_,_):=JacA;
  //A::linearModelMatrixes = linearModelMatrixes; get first matrix
  //HessA = generateSymbolicHessianA(A);
  print("\n\nOutput after linearization\n\n");
  BackendDump.printBackendDAE(outHessian);
  //outHessian:=transformJacobian(outHessian,lambda); //Umbauen des Gleichungssystems der jacobi matrix damit dann JAcobi erneut verwendet werden kann
  //outHessian:=SymbolicJacobian.symbolicJacobian(outHessian); //Erzeuge jetzt Hessematrix dies wird dann genauso zurueck gegeben!
end generateSymbolicHessian;

protected function transformJacobian
  "Function sets the lagrange factors to the jacobian matrix
  and reduces the left side of the equation by setting it to zero."
  input BackendDAE.BackendDAE inJacobian "Jacobian after normal calculation.";
  input list<Real> lambda "Lagrange factors.";
  output BackendDAE.BackendDAE outJacobian "Jacobian with leftside zero and multiplication of the lagrange factors.";
protected
  BackendDAE.EqSystem eqSys;
  BackendDAE.EquationArray eqns;
algorithm
  outJacobian:=reduceJacobian(inJacobian);
  outJacobian:=lagrangeJacobian(outJacobian,lambda);
  eqSys:=listGet(outJacobian.eqs,1);
  eqns:=eqSys.orderedEqs;
end transformJacobian;

protected function generateSymbolicHessianA
  input  BackendDAE.SymbolicJacobian A;
  output BackendDAE.BackendDAE HessA;
algorithm
  // somehow get states
  /*
   lambdas = if true then SOME(getLambdaList(listLength(states))) else NONE();

  if isSome(lambdas) then
    linearModelMatrix = multiplyLambdas(lambdas, linearModelMatrix);
  end if;
*/

// add up equations to one equation
// differentiate equation wrt all states
end generateSymbolicHessianA;

protected function reduceJacobian
input BackendDAE.BackendDAE inJacobian;
output BackendDAE.BackendDAE outReducedJacobian;
algorithm
  outReducedJacobian:=inJacobian;
end reduceJacobian;

protected function lagrangeJacobian
 input BackendDAE.BackendDAE reducedJacobian;
 input list<Real> lambda;
 output BackendDAE.BackendDAE lagrangeGradient;
algorithm
  lagrangeGradient:=reducedJacobian;
end lagrangeJacobian;

protected function getLambdaList
  input Integer lambdaCount;
  output list<DAE.ComponentRef> lambdas = {};
algorithm
  for i in lambdaCount:-1:1 loop
    lambdas := DAE.CREF_IDENT("$lambda", DAE.T_REAL_DEFAULT, {DAE.INDEX(DAE.ICONST(i))}) ::lambdas;
  end for;
end getLambdaList;

protected function multiplyLambdas
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
annotation(__OpenModelica_Interface="backend");
end SymbolicHessian;
