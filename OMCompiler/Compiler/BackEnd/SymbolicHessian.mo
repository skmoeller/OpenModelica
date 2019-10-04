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


 public import BackendDAE;
 public import DAE;

 protected
 import Array;
 import BackendDAEUtil;
 import BackendDump;
 import BackendEquation;
 import BackendVariable;
 import ComponentReference;
 import SymbolicJacobian;
 import ExpandableArray;
 import Expression;
 import Flags;
 import List;
 import Util;

 // =============================================================================
 //
 //
 //
 //
 // =============================================================================

public function generateSymbolicHessian
  "Function to generate the symbolic hessian with respect to the stats of an dynamic optimization problem."
  input BackendDAE.BackendDAE inBackendDAE "Input BackendDAE";
  output BackendDAE.BackendDAE outHessian "second derivates-> this is the hessian"; //Improve it by using special hessian struct -> need to added!!!
protected
  BackendDAE.SymbolicJacobians linearModelMatrixes; //All Matrices A,B,C,D
  list< Option< BackendDAE.BackendDAE > > SymbolicHessians = {};
  Option< list< DAE.ComponentRef > > lambdas; //Lagrange factors
algorithm
  outHessian := SymbolicJacobian.generateSymbolicLinearizationPast(inBackendDAE); //Generates Matrices A,B,C,D and calculates the second derivative
  linearModelMatrixes := BackendDAEUtil.getSharedSymJacs(outHessian.shared); //Get the Matrices from shared
  for jacobian in linearModelMatrixes loop
    _ := match jacobian
    local
      BackendDAE.SymbolicJacobian symJac;
    case ((SOME(symJac),_,_))
    equation
      SymbolicHessians = wrapperCreateSymbolicHessian(symJac)::SymbolicHessians;//Generate the Hessians by adding the lambdas and add up all equations -> write in list!!!
    then "";
    else then "";
    end match;
  end for;
  SymbolicHessians := listReverse(SymbolicHessians);
end generateSymbolicHessian;

protected function wrapperCreateSymbolicHessian
  input BackendDAE.SymbolicJacobian InSymJac;
  output Option< BackendDAE.BackendDAE > OutSymJac;
algorithm
  OutSymJac := match InSymJac
    local String nameMatrix;
          BackendDAE.BackendDAE dae;
          list<BackendDAE.Var> states;
    case (dae,nameMatrix,states,_,_,_) guard nameMatrix == "A" then SOME(createSymbolicHessian(dae, nameMatrix,states));
  else then NONE();
  end match;
end wrapperCreateSymbolicHessian;

protected function createSymbolicHessian
  "Function sets the lagrange factors and multiplies the vector to the jacobian.
   Then it runs the jacobian routine again!"
  input BackendDAE.BackendDAE JacDAE "Symbolic Jacobian Matrix";
  input String nameMatrix "Name of the matrix";
  input list<BackendDAE.Var> states "States form Matrix A";
  output BackendDAE.BackendDAE Hessian "Symbolic Hessian Matrix";
protected
  Option< list< DAE.ComponentRef > > lambdas;
algorithm
  Hessian := JacDAE;
  lambdas := if Flags.getConfigBool(Flags.GENERATE_SYMBOLIC_HESSIAN) then SOME(getLambdaList(listLength(states))) else NONE();
  /*Section for setting the lambdas & add up equations to one equation*/
  if isSome(lambdas) then
    Hessian := multiplyLambdas(lambdas, Hessian, nameMatrix); //multiple the lagrange factors and add the equations
  end if;
  print("\n\n Hessian for "+nameMatrix+"\n\n");
  BackendDump.dumpDAE(Hessian);
end createSymbolicHessian;

protected function getLambdaList
  "Function sets the lambdas to the system."
  input Integer lambdaCount "Number of lambdas";
  output list< DAE.ComponentRef > lambdas = {} "List of componentrefs for lambdas";
algorithm
  for i in 1:lambdaCount loop //Iteration: num_lambda==num_stats -> for matrix A!!!
    lambdas := DAE.CREF_IDENT("$lambda", DAE.T_ARRAY_REAL_NODIM, {DAE.INDEX(DAE.ICONST(i))}) ::lambdas;
  end for;
end getLambdaList;

protected function multiplyLambdas
  "Function calculates '(lambda)^T*Jacobian' -> First multiply the lambdas, then add the equation up to one."
  input Option< list< DAE.ComponentRef > > lambdasOption; //List of lambdas
  input BackendDAE.BackendDAE jac; //Jacobi Matrix
  input String matrixName; //Name of the jacobian (A,B,C,D)
  output BackendDAE.BackendDAE lambdaJac; //Equationsystem with only one equation -> that is the hessian!
protected
  list<DAE.ComponentRef> lambdas; //Lambdas without option typ!
  BackendDAE.EquationArray eqns, jacEqns; //Array for the hessian!
  BackendDAE.EqSystem eqs;
  BackendDAE.Variables vars;
  BackendDAE.Equation eq;
  Integer indexEq, NumOfEqs;
  DAE.Exp eqExpr;
  list<DAE.Exp> hessExpr = {};
algorithm

  SOME(lambdas) := lambdasOption;
  lambdaJac := jac;
  {eqs} := lambdaJac.eqs;
  BackendDAE.EQSYSTEM(orderedVars = vars, orderedEqs = eqns) := eqs;
  jacEqns := eqns;
  NumOfEqs := ExpandableArray.getNumberOfElements(eqns);

  /*get ordered equations from jac
  traverse and multiply each lambda on rhs
  e1.rhs -> e1.rhs * lambda[1]*/
  indexEq := 1;
  for lambdaList in lambdas loop
    eq := ExpandableArray.get(indexEq, eqns);
    eqExpr := BackendEquation.getEquationRHS(eq);
    eqExpr := multiplyLambda2Expression(eqExpr, lambdaList);
    hessExpr := eqExpr::hessExpr;
    indexEq := indexEq+1;
  end for;
  (vars, eqns) := addEquations(hessExpr, eq, matrixName, vars);
  eqns := addConstraintEquation(jacEqns, eqns, NumOfEqs, indexEq);
  /*Updating the DAE*/
  eqs.orderedEqs := eqns;
  eqs.orderedVars := vars;
  lambdaJac.eqs := {eqs};
end multiplyLambdas;

protected function multiplyLambda2Expression
  "Function takes RHS of an Equation and multiplies the lagrange factor."
  input DAE.Exp inExp;
  input DAE.ComponentRef crefLambda;
  output DAE.Exp outExp;
protected
  DAE.Exp expLambda;
algorithm
  /*make cref to an expression*/
  expLambda := Expression.crefToExp(crefLambda);
  outExp := Expression.expMul(inExp, expLambda);
end multiplyLambda2Expression;

protected function addEquations
  "Helper Function for adding up all the Equations."
  input list<DAE.Exp> hessExpr;
  input BackendDAE.Equation inEq;
  input String matrixName;
  input output BackendDAE.Variables vars;
  output BackendDAE.EquationArray eqns;
protected
  DAE.Exp EqSum;
  BackendDAE.Equation eq;
algorithm
  EqSum := Expression.makeSum(hessExpr);
  eqns := ExpandableArray.new(1, inEq);
  (vars, eq) := setHessian(inEq, EqSum, matrixName, vars);
  eqns := ExpandableArray.set(1, eq, eqns);
end addEquations;

protected function setHessian
  "Function sets the new Equation for the Hessian on the RHS and creates a new Variable for the LHS -> '$HessA = lambda[i]*eq[i]+...'."
  input BackendDAE.Equation inEq;
  input DAE.Exp rhsWithLambda;
  input String matrixName;
  input output BackendDAE.Variables vars;
  output BackendDAE.Equation outEq;
protected
  DAE.ComponentRef hessCref = ComponentReference.makeCrefIdent(Util.hessianIdent + matrixName, DAE.T_ARRAY_REAL_NODIM, {});
algorithm
  outEq := match (inEq)
    local
      BackendDAE.Equation localEq;
    case localEq as (BackendDAE.EQUATION())
      equation
        // maybe wrong check with compiled code
        localEq.exp = Expression.crefToExp(hessCref);
        localEq.scalar = rhsWithLambda;
      then localEq;
    case localEq as (BackendDAE.COMPLEX_EQUATION())
      equation
        localEq.left = Expression.crefToExp(hessCref);
        localEq.right = rhsWithLambda;
      then localEq;
    case localEq as (BackendDAE.ARRAY_EQUATION())
      equation
        localEq.left = Expression.crefToExp(hessCref);
        localEq.right = rhsWithLambda;
      then localEq;
    case localEq as (BackendDAE.SOLVED_EQUATION())
      equation
        localEq.componentRef = hessCref;
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
  vars := BackendVariable.addVar(BackendVariable.makeVar(hessCref), vars);
end setHessian;

protected function addConstraintEquation
  input BackendDAE.EquationArray jacEqs;
  input output BackendDAE.EquationArray hessEqs;
  input Integer numOfEqs;
  input Integer idx;
algorithm
  for i in idx:numOfEqs loop
    hessEqs := ExpandableArray.set(i, ExpandableArray.get(i, jacEqs), hessEqs);
  end for;
end addConstraintEquation;

annotation(__OpenModelica_Interface="backend");
end SymbolicHessian;
