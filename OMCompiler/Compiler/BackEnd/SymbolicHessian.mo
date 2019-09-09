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
  BackendDAE.SymbolicJacobians linearModelMatrixes; //All Matrices A,B,C,D
  BackendDAE.SymbolicJacobian jacA; //Matrix A
  BackendDAE.SymbolicJacobian jacB; //Matrix B
  BackendDAE.SymbolicJacobian jacC; //Matrix C
  BackendDAE.SymbolicJacobian jacD; //Matrix D
  Option<list<DAE.ComponentRef>> lambdas; //Lagrange factors
algorithm
  outHessian := SymbolicJacobian.generateSymbolicLinearizationPast(inBackendDAE); //Generates Matrices A,B,C,D
  linearModelMatrixes := BackendDAEUtil.getSharedSymJacs(outHessian.shared); //Get the Matrices
  (SOME(jacA),_,_)::(SOME(jacB),_,_)::(SOME(jacC),_,_)::(SOME(jacD),_,_)::linearModelMatrixes:=linearModelMatrixes; //Isolate A,B,C,D
  //Generate the Hessian for Matrix A
  outHessian:=generateSymbolicHessianA(jacA);
end generateSymbolicHessian;

protected function generateSymbolicHessianA
  "Function sets the lagrange factors and multiplies the vector to the jacobian.
   Then it runs the jacobian routine again!"
  input  BackendDAE.SymbolicJacobian A "Symbolic Jacobian Matrix A";
  output BackendDAE.BackendDAE HessA "Symbolic Hessian Matrix for A";
protected
  list<BackendDAE.Var> stats;
  Option<list<DAE.ComponentRef>> lambdas;
  BackendDAE.SymbolicJacobians linearModelMatrixes;
  BackendDAE.SymbolicJacobian jacA;
algorithm
  // somehow get states -> wich stats are important for what matrix???
  (_,_,stats,_,_,_):=A;
   lambdas:=if Flags.getConfigBool(Flags.GENERATE_SYMBOLIC_HESSIAN) then SOME(getLambdaList(listLength(stats))) else NONE();
  if isSome(lambdas) then
    HessA := multiplyLambdas(lambdas,A);
  end if;
// add up equations to one equation
// differentiate equation wrt all states
end generateSymbolicHessianA;

protected function getLambdaList
  "Function sets the lambdas to the system."
  input Integer lambdaCount "Number of lambdas";
  output list<DAE.ComponentRef> lambdas = {} "List of componentrefs for lambdas";
algorithm
  for i in lambdaCount:-1:1 loop //Iteration: num_lambda==num_stats
    lambdas := DAE.CREF_IDENT("$lambda", DAE.T_ARRAY_REAL_NODIM, {DAE.INDEX(DAE.ICONST(i))}) ::lambdas;
  end for;
end getLambdaList;

protected function multiplyLambdas
  "Function sets the Lagrangenfactors to the jacobian"
  input Option<list<DAE.ComponentRef>> lambdas_option;
  input BackendDAE.SymbolicJacobian jac;
  output BackendDAE.BackendDAE lambdaJac;
protected
  list<DAE.ComponentRef> lambdas;
  BackendDAE.EquationArray eqns;
  BackendDAE.EqSystem eqs;
  BackendDAE.Equation eq;
  Integer indexEq;
  DAE.Exp eqExpr;
algorithm
  SOME(lambdas):=lambdas_option;
  (lambdaJac,_,_,_,_,_):=jac;
  {eqs} := lambdaJac.eqs;
  BackendDAE.EQSYSTEM(orderedEqs=eqns) := eqs;
  /*get ordered equations from jac
  traverse and multiply each lambda on rhs
  e1.rhs -> e1.rhs * lambda1*/
  indexEq:=1;
  for lambdaList in lambdas loop
    eq := ExpandableArray.get(indexEq,eqns);
    eqExpr := getExpression(eq);
    eqExpr := multiplyLambda2Expression(eqExpr,lambdaList);
    eq := setExpression(eq,eqExpr);
    eqns := ExpandableArray.update(indexEq,eq,eqns);
    indexEq:=indexEq+1;
  end for;
  /*Updating the DAE*/
  eqs.orderedEqs:=eqns;
  lambdaJac.eqs := {eqs};
end multiplyLambdas;

protected function getExpression
  input BackendDAE.Equation inEq;
  output DAE.Exp rhs;
algorithm
  rhs := match (inEq)
    local
      DAE.Exp res;
    case (BackendDAE.EQUATION(scalar = res)) then res;
    case (BackendDAE.COMPLEX_EQUATION(right = res)) then res;
    case (BackendDAE.ARRAY_EQUATION(right = res)) then res;
    case (BackendDAE.SOLVED_EQUATION(exp = res)) then res;
    case (BackendDAE.RESIDUAL_EQUATION(exp = res)) then res;
    else
    algorithm
      print("\n\n***Error, used a unknown Equationcase!***\n\n");
    then fail();
  end match;
end getExpression;

protected function setExpression
  input BackendDAE.Equation inEq;
  input DAE.Exp rhsWithLambda;
  output BackendDAE.Equation outEq;
algorithm
  outEq := match (inEq)
    local
      BackendDAE.Equation localEq;
    case localEq as (BackendDAE.EQUATION())
      equation
        localEq.scalar = rhsWithLambda;
      then localEq;
    case localEq as (BackendDAE.COMPLEX_EQUATION())
      equation
        localEq.right = rhsWithLambda;
      then localEq;
    case localEq as (BackendDAE.ARRAY_EQUATION())
      equation
        localEq.right = rhsWithLambda;
      then localEq;
    case localEq as (BackendDAE.SOLVED_EQUATION())
      equation
        localEq.exp = rhsWithLambda;
      then localEq;
    case localEq as (BackendDAE.RESIDUAL_EQUATION())
      equation
        localEq.exp = rhsWithLambda;
      then localEq;
    else
      equation
        localEq = inEq;
      then localEq;
  end match;
end setExpression;

protected function multiplyLambda2Expression
  input DAE.Exp inExp;
  input DAE.ComponentRef crefLambda;
  output DAE.Exp outExp;
protected
  DAE.Exp expLambda;
algorithm
  /*make cref to an expression*/
  expLambda := Expression.crefToExp(crefLambda);
  outExp := Expression.expMul(inExp,expLambda);
end multiplyLambda2Expression;
annotation(__OpenModelica_Interface="backend");
end SymbolicHessian;
