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
  outHessian := SymbolicJacobian.generateSymbolicLinearizationPast(inBackendDAE); //Generates Matrices A,B,C,D and calculates the second derivativ
  linearModelMatrixes := BackendDAEUtil.getSharedSymJacs(outHessian.shared); //Get the Matrices from shared
  for jacobian in linearModelMatrixes loop
    _ := match jacobian
    local
      BackendDAE.SymbolicJacobian symJac;
    case ((SOME(symJac),_,_))
    equation
      SymbolicHessians = SOME(createSymbolicHessian(symJac))::SymbolicHessians;//Generate the Hessians by adding the lambdas and add up all equations -> write in list!!!
    then "";
    else then "";
    end match;
  end for;
  SymbolicHessians := listReverse(SymbolicHessians);
end generateSymbolicHessian;

protected function createSymbolicHessian
  "Function sets the lagrange factors and multiplies the vector to the jacobian.
   Then it runs the jacobian routine again!"
  input  BackendDAE.SymbolicJacobian Jac "Symbolic Jacobian Matrix";
  output BackendDAE.BackendDAE Hessian "Symbolic Hessian Matrix";
protected
  Option< list< DAE.ComponentRef > > lambdas;
  String nameMatrix;
  BackendDAE.SymbolicJacobians linearModelMatrixes;
  list< BackendDAE.Equation > eqns, derivedEquations;
  DAE.FunctionTree functions;
algorithm
  (Hessian, nameMatrix,_,_,_,_) := Jac;
  _ := match Jac
      local list<BackendDAE.Var> states;
    case (_,_,states,_,_,_)
    equation
      lambdas = if Flags.getConfigBool(Flags.GENERATE_SYMBOLIC_HESSIAN) then SOME(getLambdaList(listLength(states))) else NONE();
    then "";
    else then "";
    end match;
  /*Section for setting the lambdas & add up equations to one equation -> New vars!!!*/
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
  for i in lambdaCount:-1:1 loop //Iteration: num_lambda==num_stats
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
  BackendDAE.EquationArray eqns; //Array for the hessian!
  BackendDAE.EqSystem eqs; //
  BackendDAE.Equation eq;
  Integer indexEq;
  DAE.Exp eqExpr;
  list<DAE.Exp> hessExpr = {};
algorithm
  SOME(lambdas) := lambdasOption;
  lambdaJac := jac;
  {eqs} := lambdaJac.eqs;
  BackendDAE.EQSYSTEM(orderedEqs = eqns) := eqs;
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
  eqns := addEquations(hessExpr, eq, matrixName);
  /*Updating the DAE*/
  eqs.orderedEqs := eqns;
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
  output BackendDAE.EquationArray eqns;
protected
  DAE.Exp EqSum;
  BackendDAE.Equation eq;
algorithm
  EqSum := Expression.makeSum(hessExpr);
  eqns := ExpandableArray.new(1, inEq);
  eq := setHessian(inEq, EqSum, matrixName);
  eqns := ExpandableArray.set(1, eq, eqns);
end addEquations;

protected function setHessian
  "Function sets the new Equation for the Hessian on the RHS and creates a new Variable for the LHS -> '$HessA = lambda[i]*eq[i]+...'."
  input BackendDAE.Equation inEq;
  input DAE.Exp rhsWithLambda;
  input String matrixName;
  output BackendDAE.Equation outEq;
algorithm
  outEq := match (inEq)
    local
      BackendDAE.Equation localEq;
    case localEq as (BackendDAE.EQUATION())
      equation
        // maybe wrong check with compiled code
        localEq.exp = Expression.crefToExp(ComponentReference.makeCrefIdent(Util.hessianIdent + matrixName, DAE.T_ARRAY_REAL_NODIM, {}));
        localEq.scalar = rhsWithLambda;
      then localEq;
    case localEq as (BackendDAE.COMPLEX_EQUATION())
      equation
        localEq.left = Expression.crefToExp(ComponentReference.makeCrefIdent(Util.hessianIdent + matrixName, DAE.T_ARRAY_REAL_NODIM, {}));
        localEq.right = rhsWithLambda;
      then localEq;
    case localEq as (BackendDAE.ARRAY_EQUATION())
      equation
        localEq.left = Expression.crefToExp(ComponentReference.makeCrefIdent(Util.hessianIdent + matrixName, DAE.T_ARRAY_REAL_NODIM, {}));
        localEq.right = rhsWithLambda;
      then localEq;
    case localEq as (BackendDAE.SOLVED_EQUATION())
      equation
        localEq.componentRef = ComponentReference.makeCrefIdent(Util.hessianIdent + matrixName, DAE.T_ARRAY_REAL_NODIM, {});
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
end setHessian;

annotation(__OpenModelica_Interface="backend");
end SymbolicHessian;
