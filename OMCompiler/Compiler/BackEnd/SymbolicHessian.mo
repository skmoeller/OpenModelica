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
 import ExpandableArray;
 import Expression;
 import ExpressionDump;
 import Flags;
 import SymbolicJacobian;
 import List;
 import Util;

 // =============================================================================
 //                         !!! Place your ad here !!!
 // =============================================================================

public function generateSymbolicHessian
  "Function to generate the symbolic hessian with respect to the stats of an dynamic optimization problem."
  input BackendDAE.BackendDAE inBackendDAE;
  output BackendDAE.BackendDAE outBackendDAE;
protected
  BackendDAE.EqSystems eqs;
  BackendDAE.Shared shared;
  BackendDAE.SymbolicJacobians symJacs;
  BackendDAE.SymbolicHessians symHesss = {};
algorithm
  BackendDAE.DAE(eqs=eqs,shared=shared) := inBackendDAE;
  symJacs := shared.symjacs;
  for jacobian in symJacs loop
    _ := match jacobian
    local
      BackendDAE.SymbolicJacobian symJac;
    case ((SOME(symJac),_,_))
    algorithm
      symHesss := createSymbolicHessian(symJac)::symHesss; //multiply lambdas then derive second time
    then "";
  else then "";
    end match;
  end for;
  /*Set list to correct order*/
  symHesss := listReverse(symHesss);
  /*Dump the dae of the hessian*/
  if Flags.isSet(Flags.DUMP_HESSIAN) then
    printHessian(symHesss);
  end if;
  /*Update shared*/
  shared := BackendDAEUtil.setSharedSymHesss(shared, symHesss);
  outBackendDAE := BackendDAE.DAE(eqs,shared);
end generateSymbolicHessian;

protected function createSymbolicHessian
  "Function creates the symbolic Hessians for the Jacobians A, B and C (D is not needed)
   Matrix A :  Differentiate the eqns system  w.r.t. states
   Matrix B :  Differentiate the eqns system including constraints  w.r.t. states & inputs
   Matrix C :  Differentiate the eqns system including constraints & mayer & lagrange term w.r.t. states & inputs
   "
  input BackendDAE.SymbolicJacobian InSymJac "Symbolic Jacobian Matrix";
  output Option<BackendDAE.SymbolicHessian> Hessian "Symbolic Hessian Matrix";
algorithm
  Hessian := match InSymJac
    local
      BackendDAE.BackendDAE backendDAE, hessDae;
      BackendDAE.SymbolicJacobian symjac;
      list< BackendDAE.Var > diffVars ,diffedVars, allDiffedVars;
      BackendDAE.Variables v,globalKnownVars,statesarr,inputvarsarr,paramvarsarr, optimizer_vars;
      list<BackendDAE.Var> lambdaVars, varlst, knvarlst,  states, inputvars, paramvars, states_inputs;
      DAE.FunctionTree funcs, functionTree;
      String nameMatrix, nameHessian;

    case (backendDAE,nameMatrix,diffVars,diffedVars,allDiffedVars,_) guard stringEqual(nameMatrix,"A") or stringEqual(nameMatrix,"B") or stringEqual(nameMatrix,"C")
    algorithm

      hessDae := BackendDAEUtil.copyBackendDAE(backendDAE);
      (hessDae,lambdaVars) := multiplyLambdas(hessDae, nameMatrix); //multiple the lagrange factors
      hessDae := BackendDAEOptimize.collapseIndependentBlocks(hessDae);
      hessDae := BackendDAEUtil.transformBackendDAE(hessDae,SOME((BackendDAE.NO_INDEX_REDUCTION(),BackendDAE.EXACT())),NONE(),NONE());
      BackendDAE.DAE({BackendDAE.EQSYSTEM(orderedVars = v)}, BackendDAE.SHARED(globalKnownVars = globalKnownVars)) := hessDae;

      // Prepare all needed variables
      nameHessian := (nameMatrix+"1");

      varlst := BackendVariable.varList(v);
      states := List.select(diffVars, BackendVariable.isStateVar);
      //BackendDump.dumpVarList(states, "states");
      knvarlst := BackendVariable.varList(globalKnownVars);
      //BackendDump.dumpVarList(knvarlst,"knvarlst");
      inputvars := List.select(allDiffedVars,BackendVariable.isInput);
      //BackendDump.dumpVarList(inputvars,"inputvars");
      paramvars := List.select(knvarlst, BackendVariable.isParam);
      //BackendDump.dumpVarList(paramvars,"paramvars");

      states_inputs := diffVars;
      //BackendDump.dumpVarList(states_inputs, "states_inputs");

      statesarr := BackendVariable.listVar1(diffedVars);
      //print("statesarr hes\n");
      //BackendDump.printVariables(statesarr);
      //print("\n\n");
      inputvarsarr := BackendVariable.listVar1(inputvars);
      //print("statesarr hes\n");
      //BackendDump.printVariables(statesarr);
      //print("\n\n")
      paramvarsarr := BackendVariable.listVar1(paramvars);
      //print("statesarr hes\n");
      //BackendDump.printVariables(statesarr);
      //print("\n\n")
      optimizer_vars := BackendVariable.listVar1(diffedVars);
      //print("statesarr hes\n");
      //BackendDump.printVariables(statesarr);
      //print("\n\n")

      //print("System for HESSIAN\n\n");
      //print("DAE Hes\n");
      //BackendDump.dumpDAE(hessDae);
      //print("\n");
      //BackendDump.dumpVarList(states_inputs, "states_inputs hes");
      //print("statesarr hes\n\n");
      //BackendDump.printVariables(statesarr);
      //print("inputvarsarr hes\n\n");
      //BackendDump.printVariables(inputvarsarr);
      //print("paramvarsarr hes\n\n");
      //BackendDump.printVariables(paramvarsarr);
      //print("optimizer_vars hes\n\n");
      //BackendDump.printVariables(optimizer_vars);
      //BackendDump.dumpVarList(varlst, "varlst hes");

      (SOME(symjac), functionTree,_,_) := SymbolicJacobian.generateGenericJacobian(hessDae, states_inputs, statesarr, inputvarsarr, paramvarsarr, optimizer_vars, varlst ,nameHessian, false, true); //generate second derivates
      (hessDae,_,_,diffedVars,allDiffedVars,_) := symjac;

      hessDae := BackendDAEUtil.setFunctionTree(hessDae, functionTree);
      hessDae := setHessianMatrix(hessDae,nameMatrix); //add up the equations
    then SOME((hessDae,nameMatrix,InSymJac,diffVars,diffedVars,allDiffedVars,lambdaVars));

    else NONE();
  end match;
end createSymbolicHessian;

protected function multiplyLambdas
  "Function calculates '(lambda)^T*Jacobian'-> multiply lambda[i] to the equation i."
  input BackendDAE.BackendDAE jac; //Jacobi Matrix
  input String matrixName; //Name of the jacobian (A,B,C)
  output BackendDAE.BackendDAE lambdaJac; //Equationsystem with only one equation -> that is the hessian!
  output list<BackendDAE.Var> lambdaVars = {}; //All the lambda variables
protected
  Option< list< DAE.ComponentRef > > lambdasOption;//Lambdas with option
  list<DAE.ComponentRef> lambdas; //Lambdas without option typ!
  BackendDAE.EquationArray eqns, jacEqns; //Array for the hessian!
  BackendDAE.EqSystem eqs;
  BackendDAE.Shared shared;
  BackendDAE.Equation eq;
  DAE.Exp eqExpr;
  list<BackendDAE.Equation> innerEqns = {}, residualEqns = {}, lambdaEqns ={};
algorithm
  /*more or less the dae's should be the same, just multiple a lambda*/
  lambdaJac := jac;
  //BackendDump.dumpDAE(jac);
  /*get the ordered equations*/
  {eqs} := lambdaJac.eqs;
  BackendDAE.EQSYSTEM(orderedEqs = eqns) := eqs;
  jacEqns := eqns;

  shared := lambdaJac.shared;

  (innerEqns,residualEqns) := BackendEquation.traverseEquationArray(eqns, assignEqnToInnerOrResidualFirstDerivatives, (innerEqns, residualEqns));
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
      eq::residualEqns := residualEqns;
      eqExpr := BackendEquation.getEquationRHS(eq);
      eqExpr := multiplyLambda2Expression(eqExpr, lambdaList);
      eq := BackendEquation.setEquationRHS(eq, eqExpr);
      lambdaEqns := eq::lambdaEqns;
      /*Create the lambda Vars*/
      lambdaVars := createLambdaVar(lambdaList)::lambdaVars;
    end for;
    //BackendDump.dumpEquationList(lambdaEqns,"lambda Equations1");
    lambdaEqns := listAppend(lambdaEqns,innerEqns);
    //BackendDump.dumpEquationList(lambdaEqns,"lambda Equations2");
    eqns := BackendEquation.listEquation(lambdaEqns);
    /*Updating the DAE*/
    eqs.orderedEqs := eqns;
    lambdaJac.eqs := {eqs};
    shared.globalKnownVars := BackendVariable.addVars(lambdaVars,shared.globalKnownVars);
    lambdaJac.shared := shared;
  else
    print("\n***Error***\n NO RESIDUAL EQUATIONS GIVEN!\n\n");
    fail();
  end if;
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

protected function createLambdaVar
  "Creates for given cref 'lambda[i]' a BackendDAE.Var, with kind: PARAM()"
  input DAE.ComponentRef lambdaCref;
  output BackendDAE.Var lambdaVar;
algorithm
  lambdaVar := BackendDAE.VAR(lambdaCref, BackendDAE.PARAM(), DAE.INPUT(), DAE.NON_PARALLEL(), ComponentReference.crefLastType(lambdaCref), NONE(), NONE(), {}, DAE.emptyElementSource, NONE(), NONE(), DAE.BCONST(false), NONE(),DAE.NON_CONNECTOR(), DAE.NOT_INNER_OUTER(), true);
end createLambdaVar;

protected function setHessianMatrix
  "Function add all equations for the the new Variable for the LHS -> '$HessA = lambda[i]*eq[i]+...'."
  input BackendDAE.BackendDAE inHessDae;
  input String matrixName;
  output BackendDAE.BackendDAE outHessDae;
protected
  BackendDAE.EquationArray eqns, hessEqns; //Array for the hessian!
  BackendDAE.EqSystem eqs;
  BackendDAE.Shared shared;
  BackendDAE.Variables vars;
  BackendDAE.Equation eq;
  DAE.Exp eqExpr;
  list<DAE.Exp> hessExpr = {};
  list<BackendDAE.Equation> mayerLagrange = {}, residualEqns = {};
  DAE.ComponentRef hessCref = ComponentReference.makeCrefIdent(Util.hessianIdent + matrixName, DAE.T_REAL_DEFAULT, {});
algorithm
  outHessDae := inHessDae;

  {eqs} := outHessDae.eqs;
  BackendDAE.EQSYSTEM(orderedVars = vars, orderedEqs = eqns) := eqs;
  hessEqns := eqns;

  shared := outHessDae.shared;

  (mayerLagrange,residualEqns) := BackendEquation.traverseEquationArray(eqns, assignEqnToInnerOrResidualSecondDerivatives, (mayerLagrange, residualEqns));
  //BackendDump.dumpEquationArray(eqns, "filtered eqns");
  //BackendDump.dumpEquationList(residualEqns, "residualEqns");
  //BackendDump.dumpEquationList(mayerLagrange, "mayerLagrange");

  /*get ordered equations from the given dae
  traverse and sum up in one variable*/
  if listLength(residualEqns)<>0 then
    for equ in residualEqns loop
      eqExpr := BackendEquation.getEquationRHS(equ);
      hessExpr := eqExpr::hessExpr;
      eq := equ;
    end for;
  else
    hessExpr := DAE.RCONST(real= 0.0)::hessExpr;
    eq := BackendDAE.EQUATION(exp= DAE.RCONST(real= 0.0), scalar= DAE.RCONST(real= 0.0) ,source= DAE.emptyElementSource, attr= BackendDAE.EQ_ATTR_DEFAULT_UNKNOWN);
  end if;
  (vars, eqns) := addEquations(hessExpr, eq, matrixName, vars);
  vars := removeStateVars(vars, matrixName);
  if matrixName=="C" then
    BackendEquation.addList(mayerLagrange, eqns);
  end if;
  /*Updating the DAE*/
  eqs.orderedEqs := eqns;
  eqs.orderedVars := vars;
  eqs.matching := BackendDAE.NO_MATCHING();
  outHessDae.eqs := {eqs};
  /*Reset of the dae*/
  outHessDae :=  BackendDAEUtil.transformBackendDAE(outHessDae,SOME((BackendDAE.NO_INDEX_REDUCTION(),BackendDAE.EXACT())),NONE(),NONE());
end setHessianMatrix;

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
  //print("Hessian: "+ExpressionDump.printExpStr(EqSum)+"\n");
  //BackendDump.dumpVariables(vars,"Given Vars");
  eqns := ExpandableArray.new(1, inEq);
  (vars, eq) := setHessianElement(inEq, EqSum, matrixName, vars);
  eqns := ExpandableArray.set(1, eq, eqns);
end addEquations;

protected function setHessianElement
  "Function sets the new Equation for the Hessian on the RHS and creates a new Variable for the LHS -> '$HessA = lambda[i]*eq[i]+...'."
  input BackendDAE.Equation inEq;
  input DAE.Exp rhs;
  input String matrixName;
  input output BackendDAE.Variables vars;
  output BackendDAE.Equation outEq;
protected
  DAE.ComponentRef hessCref = ComponentReference.makeCrefIdent(Util.hessianIdent + matrixName, DAE.T_REAL_DEFAULT, {});
  BackendDAE.Var hessVar;
algorithm
  outEq := match (inEq)
    local
      BackendDAE.Equation localEq;
      Boolean isDer;
    case localEq as (BackendDAE.EQUATION())
      equation
        // maybe wrong check with compiled code
        localEq.exp = Expression.crefToExp(hessCref);
        localEq.scalar = rhs;
      then localEq;
    case localEq as (BackendDAE.COMPLEX_EQUATION())
      equation
        localEq.left = Expression.crefToExp(hessCref);
        localEq.right = rhs;
      then localEq;
    case localEq as (BackendDAE.ARRAY_EQUATION())
      equation
        localEq.left = Expression.crefToExp(hessCref);
        localEq.right = rhs;
      then localEq;
    case localEq as (BackendDAE.SOLVED_EQUATION())
      equation
        localEq.componentRef = hessCref;
        localEq.exp = rhs;
      then localEq;
    case localEq as (BackendDAE.RESIDUAL_EQUATION())
      equation
        localEq.exp = rhs;
      then localEq;
    else
      equation
        localEq = inEq;
      then localEq;
  end match;
  hessVar := BackendVariable.makeVar(hessCref);
  hessVar := BackendVariable.setVarKind(hessVar,BackendDAE.HESS_VAR());
  vars := BackendVariable.addVar(hessVar, vars);
end setHessianElement;

protected function assignEqnToInnerOrResidualFirstDerivatives
  "Functions splits the equation in a list of residual and inner equations (For multiply the lambdas to it -> just need $DER..)"
  input output BackendDAE.Equation eqn;
  input output tuple<list<BackendDAE.Equation>, list<BackendDAE.Equation>> eqnTpl;
protected
  list<BackendDAE.Equation> innerEqns, residualEqns;
algorithm
  (innerEqns, residualEqns) := eqnTpl;
  if not isInnerEqn(eqn) then
    residualEqns := eqn :: residualEqns;
  else
    innerEqns := eqn :: innerEqns;
  end if;
  eqnTpl := (innerEqns, residualEqns);
end assignEqnToInnerOrResidualFirstDerivatives;

protected function assignEqnToInnerOrResidualSecondDerivatives
  "Functions splits the equation in a list of residual and inner equations (needed for the second derivatives)"
  input output BackendDAE.Equation eqn;
  input output tuple<list<BackendDAE.Equation>, list<BackendDAE.Equation>> eqnTpl;
protected
  list<BackendDAE.Equation> mayerLagrange, residualEqns;
algorithm
  (mayerLagrange, residualEqns) := eqnTpl;
  if isResidualEqn(eqn) and not isObjectiveFunction(eqn) then
    residualEqns := eqn :: residualEqns;
  elseif isObjectiveFunction(eqn) then
    mayerLagrange := eqn :: mayerLagrange;
  end if;
  eqnTpl := (mayerLagrange, residualEqns);
end assignEqnToInnerOrResidualSecondDerivatives;

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
      guard(ComponentReference.isSecondPartialDerivativeHessian(cr) or ComponentReference.isMayerOrLagrange(cr))
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

protected function isObjectiveFunction
  "Functions checks if a given Equation is Mayer or Lagrange Term "
  input BackendDAE.Equation eqn;
  output Boolean b;
algorithm
  b := match eqn
    local
      DAE.ComponentRef cr;
    /* other eqn types relevant? */
    case BackendDAE.EQUATION(exp = DAE.CREF(componentRef = cr))
      guard(Util.stringStartsWith("$OMC$",ComponentReference.crefFirstIdent(cr)))
    then true;
  else false;
  end match;
end isObjectiveFunction;

protected function removeStateVars
  "Functions removes the first derivatives from the variables"
  input output BackendDAE.Variables variables;
  input String matrixName;
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
          if not (matrixName=="C") then
            variables = BackendVariable.deleteVar(cr, variables);
          elseif not (ComponentReference.isMayerOrLagrange(cr)) then
            variables = BackendVariable.deleteVar(cr, variables);
          end if;
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
        (dae,matrixName,_,_,_,_,_) = symHe;
        print("\n\n########################################\nHessian for "+matrixName+"\n########################################\n\n");
        BackendDump.dumpDAE(dae);
      then "";
      case NONE() then "";
    end match;
  end for;
end printHessian;

annotation(__OpenModelica_Interface="backend");
end SymbolicHessian;
