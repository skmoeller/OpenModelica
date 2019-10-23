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
 import DynamicOptimization;
 import SymbolicJacobian;
 import ExpandableArray;
 import Expression;
 import Flags;
 import List;
 import Util;

 // =============================================================================
 //                         !!! Place your ad here !!!
 // =============================================================================

public function generateSymbolicHessian
  "Function to generate the symbolic hessian with respect to the stats of an dynamic optimization problem."
  input BackendDAE.BackendDAE inBackendDAE "Input BackendDAE";
  output BackendDAE.BackendDAE outHessian "second derivates-> this is the hessian";
protected
  BackendDAE.SymbolicJacobians linearModelMatrixes; //All Matrices A,B,C,D
  BackendDAE.SymbolicHessians symHesss = {};
  Option< list< DAE.ComponentRef > > lambdas; //Lagrange factors
  /*Stuff for the Jacobian*/
  BackendDAE.EqSystems eqs;
  BackendDAE.Shared shared;
  DAE.FunctionTree funcs, functionTree;
algorithm
  BackendDAE.DAE(eqs=eqs,shared=shared) := inBackendDAE;
  (linearModelMatrixes, funcs) := SymbolicJacobian.createLinearModelMatrixes(inBackendDAE, true ,true);

  for jacobian in linearModelMatrixes loop
    _ := match jacobian
    local
      BackendDAE.SymbolicJacobian symJac;
    case ((SOME(symJac),_,_))
    equation
      symHesss = createSymbolicHessian(symJac)::symHesss;//Generate the Hessians by adding the lambdas and add up all equations -> write in list!!!
    then "";
    else then "";
    end match;
  end for;
  symHesss := listReverse(symHesss);

  if Flags.isSet(Flags.DUMP_HESSIAN) then
    printHessian(symHesss);
  end if;

  shared := BackendDAEUtil.setSharedSymHesss(shared, symHesss);
  functionTree := BackendDAEUtil.getFunctions(shared);
  functionTree := DAE.AvlTreePathFunction.join(functionTree, funcs);
  shared := BackendDAEUtil.setSharedFunctionTree(shared, functionTree);
  outHessian := BackendDAE.DAE(eqs,shared);
end generateSymbolicHessian;

protected function createSymbolicHessian
  "Function creates the symbolic Hessians for the Jacobians A, B and C"
  input BackendDAE.SymbolicJacobian InSymJac "Symbolic Jacobian Matrix";
  output Option< BackendDAE.SymbolicHessian > Hessian "Symbolic Hessian Matrix";
algorithm
  Hessian := match InSymJac
    local
    BackendDAE.BackendDAE dae;
    String nameMatrix;
    case (dae,nameMatrix,_,_,_,_) guard nameMatrix == "A" or nameMatrix == "B" or nameMatrix == "C" then SOME((multiplyLambdas(dae, nameMatrix), nameMatrix)); //multiple the lagrange factors and add the equations
    else then NONE();
  end match;
end createSymbolicHessian;

protected function multiplyLambdas
  "Function calculates '(lambda)^T*Jacobian' -> First multiply the lambdas, then add the equation up to one."
  input BackendDAE.BackendDAE jac; //Jacobi Matrix
  input String matrixName; //Name of the jacobian (A,B,C,D)
  output BackendDAE.BackendDAE lambdaJac; //Equationsystem with only one equation -> that is the hessian!
protected
  Option< list< DAE.ComponentRef > > lambdasOption;//Lambdas with option
  list<DAE.ComponentRef> lambdas; //Lambdas without option typ!
  BackendDAE.EquationArray eqns, jacEqns; //Array for the hessian!
  BackendDAE.EqSystem eqs;
  BackendDAE.Variables vars;
  BackendDAE.Equation eq;
  DAE.Exp eqExpr;
  list<DAE.Exp> hessExpr = {};
  list<BackendDAE.Equation> innerEqns = {}, residualEqns = {};
algorithm
  lambdaJac := jac;
  {eqs} := lambdaJac.eqs;
  BackendDAE.EQSYSTEM(orderedVars = vars, orderedEqs = eqns) := eqs;
  jacEqns := eqns;
  (innerEqns,residualEqns) := BackendEquation.traverseEquationArray(eqns, assignEqnToInnerOrResidual, (innerEqns, residualEqns));
  //BackendDump.dumpEquationArray(eqns, "filtered eqns");
  lambdasOption := if Flags.getConfigBool(Flags.GENERATE_SYMBOLIC_HESSIAN) and (listLength(residualEqns)<>0) then SOME(getLambdaList(listLength(residualEqns))) else NONE();

  if isSome(lambdasOption) then
    SOME(lambdas) := lambdasOption;
    //BackendDump.dumpEquationList(innerEqns, "inner Equations");
    //BackendDump.dumpEquationList(residualEqns, "residualEqns");

    /*get ordered equations from jac
    traverse and multiply each lambda on rhs
    e1.rhs -> e1.rhs * lambda[1]*/
    for lambdaList in lambdas loop
      eq :: residualEqns := residualEqns;
      eqExpr := BackendEquation.getEquationRHS(eq);
      eqExpr := multiplyLambda2Expression(eqExpr, lambdaList);
      hessExpr := eqExpr::hessExpr;
    end for;
    (vars, eqns) := addEquations(hessExpr, eq, matrixName, vars);
    vars := removeStateVars(vars);
    eqns := BackendEquation.addList(innerEqns, eqns);
    /*Updating the DAE*/
    eqs.orderedEqs := eqns;
    eqs.orderedVars := vars;
    eqs.matching := BackendDAE.NO_MATCHING();
    lambdaJac.eqs := {eqs};
  else
    print("\n***Error***\n NO RESIDUAL EQUATIONS GIVEN!\n\n");
    fail();
  end if;

  lambdaJac :=  BackendDAEUtil.transformBackendDAE(lambdaJac,SOME((BackendDAE.NO_INDEX_REDUCTION(),BackendDAE.EXACT())),NONE(),NONE());
end multiplyLambdas;

protected function getLambdaList
  "Function sets the lambdas to the system."
  input Integer lambdaCount "Number of lambdas";
  output list< DAE.ComponentRef > lambdas = {} "List of componentrefs for lambdas";
algorithm
  for i in 1:lambdaCount loop //Iteration: num_lambda==(num_stats+opt. vars) -> for all the matrixes!!!
    lambdas := DAE.CREF_IDENT("$lambda", DAE.T_ARRAY_REAL_NODIM, {DAE.INDEX(DAE.ICONST(i))}) ::lambdas;
  end for;
end getLambdaList;

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
  DAE.ComponentRef hessCref = ComponentReference.makeCrefIdent(Util.hessianIdent + matrixName, DAE.T_REAL_DEFAULT, {});
algorithm
  outEq := match (inEq)
    local
      BackendDAE.Equation localEq;
      Boolean isDer;
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

protected function assignEqnToInnerOrResidual
  "Functions splits the equation in a list of residual and inner equations (needed for the second derivatives)"
  input output BackendDAE.Equation eqn;
  input output tuple<list<BackendDAE.Equation>, list<BackendDAE.Equation>> eqnTpl;
protected
  list<BackendDAE.Equation> innerEqns, residualEqns;
algorithm
  (innerEqns, residualEqns) := eqnTpl;
  if isResidualEqn(eqn) then
    residualEqns := eqn :: residualEqns;
  elseif isInnerEqn(eqn) then
    innerEqns := eqn :: innerEqns;
  end if;
  eqnTpl := (innerEqns, residualEqns);
end assignEqnToInnerOrResidual;

protected function isResidualEqn
  "Function checks if a given Equation is a normal derivativ"
  input BackendDAE.Equation eqn;
  output Boolean b;
algorithm
  b := match eqn
    local
      DAE.ComponentRef cr;
    /* other eqn types relevant? */
    case BackendDAE.EQUATION(exp = DAE.CREF(componentRef = cr))
      guard(ComponentReference.isSecondPartialDerivativeHessian(cr) or ComponentReference.isMayerOrLagrange(cr) or ComponentReference.isConstraint(cr))
    then true;
    else false;
  end match;
end isResidualEqn;

protected function isInnerEqn
  "Functions checks if a given Equation is an inner derivativ"
  input BackendDAE.Equation eqn;
  output Boolean b;
algorithm
  b := match eqn
    local
      DAE.ComponentRef cr;
    /* other eqn types relevant? */
    case BackendDAE.EQUATION(exp = DAE.CREF(componentRef = cr))
      guard(Util.stringStartsWith("$DER",ComponentReference.crefFirstIdent(cr)))
    then false;
    else true;
  end match;
end isInnerEqn;

protected function removeStateVars
  "Functions removes the first derivatives from the variables"
  input output BackendDAE.Variables variables;
protected
  BackendDAE.VariableArray varriableArray;
  Integer NumOfVars;
  array<Option<BackendDAE.Var>> varArr;
algorithm
  NumOfVars := variables.numberOfVars;
  varriableArray := variables.varArr;
  varArr := varriableArray.varOptArr;
  for i in 1:NumOfVars loop
    _ := match varArr[i]
      local
        DAE.ComponentRef cr;
      case SOME(BackendDAE.VAR(varName = cr)) guard(Util.stringStartsWith("$DER",ComponentReference.crefFirstIdent(cr)) or ComponentReference.isConstraint(cr) or ComponentReference.isMayerOrLagrange(cr))
        equation
          variables = BackendVariable.deleteVar(cr, variables);
        then "";
      else then "";
    end match;
  end for;
end removeStateVars;

protected function printHessian
  "Functions prints the Hessians A,B,C"
  input BackendDAE.SymbolicHessians symHesss;
algorithm
  for hessian in symHesss loop
    _:= match hessian
      local BackendDAE.SymbolicHessian symHe;
            BackendDAE.BackendDAE dae;
            String matrixName;
      case SOME(symHe)
      equation
        (dae,matrixName) = symHe;
        print("\n\n########################################\nHessian for "+matrixName+"\n########################################\n\n");
        BackendDump.dumpDAE(dae);
      then "";
      case NONE() then "";
    end match;
  end for;
end printHessian;

annotation(__OpenModelica_Interface="backend");
end SymbolicHessian;
