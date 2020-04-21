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
   description: This package contains methodes to calculate the symbolic Hessian Matrix. Therefore the module SymbolicJacobian.mo is mainly used."


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
 //         Section to create Hessians with the Data from a Jacobian
 // =============================================================================

public function generateSymbolicHessian
  "Function to generate the symbolic hessian with respect to the stats of an dynamic optimization problem."
  /*Input and Output DAE*/
  input BackendDAE.BackendDAE inBackendDAE;
  output BackendDAE.BackendDAE outBackendDAE;
protected
  /*Variables for DAE and Jacobian, Hessian*/
  BackendDAE.EqSystems eqs;
  BackendDAE.Shared shared;
  BackendDAE.SymbolicJacobians symJacs;
  BackendDAE.SymbolicHessians symHesss = {};
algorithm
  BackendDAE.DAE(eqs=eqs,shared=shared) := inBackendDAE;
  symJacs := shared.symjacs;
  /*Get the Jacobians from the system*/
  for jacobian in symJacs loop
    _ := match jacobian
    local
      BackendDAE.SymbolicJacobian symJac;
    case ((SOME(symJac),_,_))
    algorithm
      /*Create one Hessian and add it to the list of Hessians*/
      symHesss := createSymbolicHessian(symJac)::symHesss;
    then "";
  else then "";
    end match;
  end for;
  /*Set list to correct order*/
  symHesss := listReverse(symHesss);
  /*Dump the Hessian DAE System*/
  if Flags.isSet(Flags.DUMP_HESSIAN) then
    printHessian(symHesss);
  end if;
  /*Update shared object*/
  shared := BackendDAEUtil.setSharedSymHesss(shared, symHesss);
  outBackendDAE := BackendDAE.DAE(eqs,shared);
end generateSymbolicHessian;

protected function createSymbolicHessian
  "Function creates the symbolic Hessians for the Jacobians A, B and C (D is not needed yet)
   Matrix A :  Differentiate the eqns system  w.r.t. states
   Matrix B :  Differentiate the eqns system including constraints & lagrange term w.r.t. states & inputs
   Matrix C :  Differentiate the eqns system including constraints lagrange & mayer term w.r.t. states & inputs
   "
  /*Symbolic Jacobian Matrix*/
  input BackendDAE.SymbolicJacobian InSymJac;
  /*Symbolic Hessian Matrix*/
  output Option<BackendDAE.SymbolicHessian> Hessian;
algorithm
  Hessian := match InSymJac
    local
      /*Local Variables for the Variables of the Jacobians*/
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
      /*multiple the lagrange factors to system -> result is an Vector!*/
      (hessDae,lambdaVars) := multiplyLambdas(hessDae, nameMatrix);
      hessDae := BackendDAEOptimize.collapseIndependentBlocks(hessDae);
      hessDae := BackendDAEUtil.transformBackendDAE(hessDae,SOME((BackendDAE.NO_INDEX_REDUCTION(),BackendDAE.EXACT())),NONE(),NONE());
      BackendDAE.DAE({BackendDAE.EQSYSTEM(orderedVars = v)}, BackendDAE.SHARED(globalKnownVars = globalKnownVars)) := hessDae;

      /*Prepare all needed variables*/
      nameHessian := (nameMatrix+"1");
      varlst := BackendVariable.varList(v);
      states := List.select(diffVars, BackendVariable.isStateVar);
      knvarlst := BackendVariable.varList(globalKnownVars);
      inputvars := List.select(diffVars,BackendVariable.isInput);
      paramvars := List.select(knvarlst, BackendVariable.isParam);

      /*Set the transfer parameters for next function call*/
      states_inputs := diffVars;
      statesarr := BackendVariable.listVar1(states);
      inputvarsarr := BackendVariable.listVar1(inputvars);
      paramvarsarr := BackendVariable.listVar1(paramvars);
      optimizer_vars := BackendVariable.listVar1(diffedVars);
      
      /*Init the shared of DAE*/
      hessDae.shared := BackendDAEUtil.setSharedGlobalKnownVars(hessDae.shared,paramvarsarr);

      /*Verbose Dump for all touched Variables*/
      if Flags.isSet(Flags.DUMP_HESSIAN_VERBOSE) then
        /*Variables directly from the Jacobian*/
        BackendDump.dumpVarList(states, "States");
        BackendDump.dumpVarList(knvarlst,"Global known Variables");
        BackendDump.dumpVarList(inputvars,"Input Variables");
        BackendDump.dumpVarList(paramvars,"Parameter Variables");

        /*Structurs that will be used to derive a second time*/
        print("\n\n########################################\nHESSIAN DAE System for Matrix "+nameMatrix+"\n########################################\n\n");
        print("Full DAE System\n");
        BackendDump.dumpDAE(hessDae);
        print("\n");
        BackendDump.dumpVarList(states_inputs, "States and Input Variables");
        print("Array with the states for the Hessian:\n\n");
        BackendDump.printVariables(statesarr);
        print("Array with the Input Variables for the Hessian:\n\n");
        BackendDump.printVariables(inputvarsarr);
        print("Array with the Parameter Variables for the Hessian:\n\n");
        BackendDump.printVariables(paramvarsarr);
        print("Variables for the Optimization process\n\n");
        BackendDump.printVariables(optimizer_vars);
        BackendDump.dumpVarList(varlst, "All Variables from Equationsystem");
      end if;

      /*Derive Vector second time; thus results in the Hessian of the system*/
      (SOME(symjac), functionTree,_,_) := SymbolicJacobian.generateGenericJacobian(hessDae, states_inputs, statesarr, inputvarsarr, paramvarsarr, optimizer_vars, varlst ,nameHessian, false, true);
      (hessDae,_,_,diffedVars,allDiffedVars,_) := symjac;
      hessDae := BackendDAEUtil.setFunctionTree(hessDae, functionTree);
      /*add up the  second derivatives -> Element of Hessian*/
      hessDae := setHessianMatrix(hessDae,nameMatrix);

    then SOME((hessDae,nameMatrix,InSymJac,diffVars,diffedVars,allDiffedVars,lambdaVars));

    else NONE();
  end match;
end createSymbolicHessian;

 // =============================================================================
 //      Section to determine the Vector that will be derived a second time
 // =============================================================================

protected function multiplyLambdas
  "Function calculates Vector-Matrix product:'(lambda)^T*Jacobian' (multiplies lambda[i] to equation i)."
  /*Given Jacobian*/
  input BackendDAE.BackendDAE jac;
  input String matrixName;
  /*Result: Equationsystem with only one equation (Thats the Hessian)*/
  output BackendDAE.BackendDAE lambdaJac;
  /*The Lambda Variables -> Identifies Tensor-Level and used in IPOPT*/
  output list<BackendDAE.Var> lambdaVars = {};
protected
  Option< list< DAE.ComponentRef > > lambdasOption;
  list<DAE.ComponentRef> lambdas;
  BackendDAE.EquationArray eqns, jacEqns;
  BackendDAE.EqSystem eqs;
  BackendDAE.Shared shared;
  BackendDAE.Equation eq;
  DAE.Exp eqExpr;
  /*Different Equationtypes for the optimization*/
  list<BackendDAE.Equation> innerEqns = {}, stateEqs = {}, constraints = {}, costFunction = {}, optimizerEqns = {}, lambdaEqns = {};
  /*Adjacency Matrix and the Strong Components of the System -> get order of the Equation*/
  array<Integer> assE2V, assV2E;
  BackendDAE.StrongComponents comps;
algorithm
  lambdaJac := jac;

  if Flags.isSet(Flags.DUMP_HESSIAN_VERBOSE) then
    print("DAE-System before multiply the lambda Vector:\n\n");
    BackendDump.dumpDAE(jac);
  end if;

  {eqs} := lambdaJac.eqs;
  BackendDAE.EQSYSTEM(orderedEqs = eqns, matching = BackendDAE.MATCHING(ass1 = assE2V, ass2 = assV2E, comps = comps)) := eqs;
  jacEqns := eqns;
  eqns := BackendEquation.listEquation(BackendDAEUtil.getStrongComponentsEquations(comps, eqns));
  shared := lambdaJac.shared;

  /*Dump for the order of the system*/
  if Flags.isSet(Flags.DUMP_HESSIAN_VERBOSE) then
    print("Matched Jacobian System\n\n");
    BackendDump.dumpMatchingVars(assV2E);
    BackendDump.dumpComponents(comps);
    BackendDump.dumpEquationArray(eqns, "Ordered Equations");
  end if;

  /*Filter different types of equations: inner Equations do not get an lambda but needed for second derivatives*/
  (innerEqns, stateEqs, constraints, costFunction) := BackendEquation.traverseEquationArray(eqns, assignEqnToInnerOrResidualFirstDerivatives, (innerEqns, stateEqs, constraints, costFunction));

  /*Dump the different Equationtypes of the DAE system*/
  if Flags.isSet(Flags.DUMP_HESSIAN_VERBOSE) then
    BackendDump.dumpEquationArray(eqns, "Equationsystem from the Jacobian");
    BackendDump.dumpEquationList(stateEqs, "Equations of the states");
    BackendDump.dumpEquationList(constraints, "Constraint Equations");
    BackendDump.dumpEquationList(costFunction, "Objective function");
  end if;

  stateEqs := listReverse(stateEqs);
  constraints := listReverse(constraints);
  costFunction := listReverse(costFunction);

  optimizerEqns := listAppend(stateEqs,constraints);
  optimizerEqns := listAppend(optimizerEqns,costFunction);

  lambdasOption := if Flags.getConfigBool(Flags.GENERATE_SYMBOLIC_HESSIAN) and (listLength(optimizerEqns)<>0) then SOME(getLambdaList(listLength(optimizerEqns))) else NONE();

  if isSome(lambdasOption) then
    SOME(lambdas) := lambdasOption;
    lambdas := listReverse(lambdas);

    if Flags.isSet(Flags.DUMP_HESSIAN_VERBOSE) then
      BackendDump.dumpEquationList(innerEqns, "Inner Equations");
      BackendDump.dumpEquationList(optimizerEqns, "Equations used for the Hessian");
    end if;

    /*Traverse and multiply each lambda on RHS of single Equation (ei.rhs -> ei.rhs * lambda[i])*/
    for lambdaList in lambdas loop
      eq::optimizerEqns := optimizerEqns;

      if Flags.isSet(Flags.DUMP_HESSIAN_VERBOSE) then
        print("Lambda:= "+ComponentReference.printComponentRefStr(lambdaList)+"\n");
      end if;

      /*Take the RHS of the current Equation*/
      eqExpr := BackendEquation.getEquationRHS(eq);
      /*Multiply the lambda Expression to the RHS*/
      eqExpr := multiplyLambda2Expression(eqExpr, lambdaList);
      /*Update the RHS of current Equation*/
      eq := BackendEquation.setEquationRHS(eq, eqExpr);
      /*Put updated Equation in a list of Equation -> Representation of an Vector*/
      lambdaEqns := eq::lambdaEqns;

      /*Create Lambdas as Variable*/
      lambdaVars := createLambdaVar(lambdaList)::lambdaVars;
    end for;

    /*Dump Equations with multiplied lambda vector*/
    if Flags.isSet(Flags.DUMP_HESSIAN_VERBOSE) then
      BackendDump.dumpEquationList(lambdaEqns,"Equations with multiplied lambda");
    end if;

    /*Add the inner Equations (without lambda)*/
    lambdaEqns := listAppend(lambdaEqns,innerEqns);

    /*Dump updated Equationsystem*/
    if Flags.isSet(Flags.DUMP_HESSIAN_VERBOSE) then
      BackendDump.dumpEquationList(lambdaEqns,"Equations with multiplied lambda and inner Equations");
    end if;

    /*Update the DAE*/
    eqns := BackendEquation.listEquation(lambdaEqns);
    eqs.orderedEqs := eqns;
    lambdaJac.eqs := {eqs};
    shared.globalKnownVars := BackendVariable.addVars(lambdaVars,shared.globalKnownVars);
    lambdaJac.shared := shared;
  else
    /*No Equations that could be derive -> Error!!!*/
    print("\n************Error************\n NO EQUATIONS FOUND!\n\n");
    fail();
  end if;
end multiplyLambdas;

 // =============================================================================
 //                 Subsection to handle the lambda Variables
 // =============================================================================

protected function getLambdaList
  "Initial Function for the lambda expressions."
  /*Number of lambdas: number of Equations for optimizer (State-Equations, Constraints, Objective Function)*/
  input Integer lambdaCount;
  /*List of lambdas as Component References*/
  output list< DAE.ComponentRef > lambdas = {};
algorithm
  for i in 1:lambdaCount loop
    lambdas := DAE.CREF_IDENT("$lambda", DAE.T_ARRAY_REAL_NODIM, {DAE.INDEX(DAE.ICONST(i))}) ::lambdas;
  end for;
end getLambdaList;

protected function multiplyLambda2Expression
  "Function takes RHS of Equation 'i' and multiplies '$lambda[i]' as a factor to it."
  input DAE.Exp inExp;
  input DAE.ComponentRef crefLambda;
  output DAE.Exp outExp;
protected
  DAE.Exp expLambda;
algorithm
  expLambda := Expression.crefToExp(crefLambda);
  outExp := Expression.expMul(inExp, expLambda);
end multiplyLambda2Expression;

protected function createLambdaVar
  "Creates for the Component Reference 'lambda[i]' an BackendDAE.Var Type (kind: PARAM())"
  input DAE.ComponentRef lambdaCref;
  output BackendDAE.Var lambdaVar;
algorithm
  lambdaVar := BackendDAE.VAR(lambdaCref, BackendDAE.PARAM(), DAE.INPUT(), DAE.NON_PARALLEL(), ComponentReference.crefLastType(lambdaCref), NONE(), NONE(), {}, DAE.emptyElementSource, NONE(), NONE(), DAE.BCONST(false), NONE(),DAE.NON_CONNECTOR(), DAE.NOT_INNER_OUTER(), true);
end createLambdaVar;

 // =============================================================================
 //             Section to set the Hessian 
 // =============================================================================

protected function setHessianMatrix
  "Function add up all the equations and set Variable '$Hess' as LHS ('$HessA = lambda[i]*eq[i]+...')."
  input BackendDAE.BackendDAE inHessDae;
  input String matrixName;
  output BackendDAE.BackendDAE outHessDae;
protected
  BackendDAE.EquationArray eqns, hessEqns;
  BackendDAE.EqSystem eqs;
  BackendDAE.Shared shared;
  BackendDAE.Variables vars;
  BackendDAE.Equation eq;
  DAE.Exp eqExpr;
  list<DAE.Exp> hessExpr = {};
  list<BackendDAE.Equation> innerEqns = {}, resEqns = {};
  DAE.ComponentRef hessCref = ComponentReference.makeCrefIdent(Util.hessianIdent + matrixName, DAE.T_REAL_DEFAULT, {});
algorithm
  outHessDae := inHessDae;
  {eqs} := outHessDae.eqs;
  BackendDAE.EQSYSTEM(orderedVars = vars, orderedEqs = eqns) := eqs;
  hessEqns := eqns;
  shared := outHessDae.shared;

  /*Split the equations in inner derivatives and derivates of residual equations*/
  (innerEqns,resEqns) := BackendEquation.traverseEquationArray(eqns, assignEqnToInnerOrResidualSecondDerivatives, (innerEqns, resEqns));

  /*Dump the different types of derived equations*/
  if Flags.isSet(Flags.DUMP_HESSIAN_VERBOSE) then
    BackendDump.dumpEquationArray(eqns, "Whole Equationsystem");
    BackendDump.dumpEquationList(resEqns, "Residuale Equations");
    BackendDump.dumpEquationList(innerEqns, "Inner Derivatives");
  end if;

  /*get residuale equations from the DAE*/
  if listLength(resEqns)<>0 then
    for equ in resEqns loop
      eqExpr := BackendEquation.getEquationRHS(equ);
      hessExpr := eqExpr::hessExpr;
      eq := equ;
    end for;
  else
    hessExpr := DAE.RCONST(real= 0.0)::hessExpr;
    eq := BackendDAE.EQUATION(exp= DAE.RCONST(real= 0.0), scalar= DAE.RCONST(real= 0.0) ,source= DAE.emptyElementSource, attr= BackendDAE.EQ_ATTR_DEFAULT_UNKNOWN);
  end if;

  /*Add up the Equation and store the result in a variable*/
  (vars, eqns) := addEquations(hessExpr, eq, matrixName, vars);
  /*remove first derivatives from system*/
  vars := removeDerVars(vars, matrixName);
  eqns := BackendEquation.addList(innerEqns, eqns);

  /*Dump the updated equations and variables*/
  if Flags.isSet(Flags.DUMP_HESSIAN_VERBOSE) then
    BackendDump.dumpEquationArray(eqns, "Hessian Equations");
    BackendDump.printVariables(vars);
   end if;

  /*Update DAE*/
  eqs.orderedEqs := eqns;
  eqs.orderedVars := vars;
  eqs.matching := BackendDAE.NO_MATCHING();
  outHessDae.eqs := {eqs};
  /*Reset the causalization because of removing variables from DAE*/
  outHessDae :=  BackendDAEUtil.transformBackendDAE(outHessDae,SOME((BackendDAE.NO_INDEX_REDUCTION(),BackendDAE.EXACT())),NONE(),NONE());
end setHessianMatrix;

 // =============================================================================
 //                  Section for the Helper functions
 // =============================================================================
 
 protected function addEquations
  "Helper Function to add up the Equations for the Hessian."
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

  /*Dump the the sum with the variables of the system*/
  if Flags.isSet(Flags.DUMP_HESSIAN_VERBOSE) then
    print("Hessian for Matrix "+matrixName+" :"+ExpressionDump.printExpStr(EqSum)+"\n");
    BackendDump.dumpVariables(vars,"Variables of the DAE");
  end if;

  /*Store the Equation in an Array*/
  eqns := ExpandableArray.new(1, inEq);
  (vars, eq) := setHessianElement(inEq, EqSum, matrixName, vars);
  eqns := ExpandableArray.set(1, eq, eqns);
end addEquations;
 
protected function assignEqnToInnerOrResidualFirstDerivatives
  "Functions filters the equation that contain inner derivatives and split the other equations in resiudual, constraint and cost equations."
  input output BackendDAE.Equation eqn;
  input output tuple<list<BackendDAE.Equation>, list<BackendDAE.Equation>, list<BackendDAE.Equation>, list<BackendDAE.Equation>> eqnTpl;
protected
  list<BackendDAE.Equation> innerEqns, stateEqs, constraints, costFunction;
algorithm
  (innerEqns, stateEqs, constraints, costFunction) := eqnTpl;
  if not isInnerEqn(eqn) then
    stateEqs := eqn :: stateEqs;
  elseif isConstraint(eqn) then
    constraints := eqn :: constraints;
  elseif isObjective(eqn) then
    costFunction := eqn :: costFunction;
  else
    innerEqns := eqn :: innerEqns;
  end if;
  eqnTpl := (innerEqns, stateEqs, constraints, costFunction);
end assignEqnToInnerOrResidualFirstDerivatives;

protected function assignEqnToInnerOrResidualSecondDerivatives
  "Functions splits the equation in a list of residual equations and inner derivatives."
  input output BackendDAE.Equation eqn;
  input output tuple<list<BackendDAE.Equation>, list<BackendDAE.Equation>> eqnTpl;
protected
  list<BackendDAE.Equation> mayerLagrange, residualEqns;
algorithm
  (mayerLagrange, residualEqns) := eqnTpl;
  if isResidualEqn(eqn) or isConstraint(eqn) then
    residualEqns := eqn :: residualEqns;
  else
    mayerLagrange := eqn :: mayerLagrange;
  end if;
  eqnTpl := (mayerLagrange, residualEqns);
end assignEqnToInnerOrResidualSecondDerivatives;

protected function isConstraint
  "Functions checks if an Equation has the type constraint."
  input BackendDAE.Equation eqn;
  output Boolean b;
algorithm
  b := match eqn
    local
      DAE.ComponentRef cr;
    case BackendDAE.EQUATION(exp = DAE.CREF(componentRef = cr))
      guard(Util.stringStartsWith("$con",ComponentReference.crefFirstIdent(cr)))
    then true;
  else false;
  end match;
end isConstraint;

protected function isInnerEqn
  "Functions checks if an Equation is inner derivativ."
  input BackendDAE.Equation eqn;
  output Boolean b;
algorithm
  b := match eqn
    local
      DAE.ComponentRef cr;
    case BackendDAE.EQUATION(exp = DAE.CREF(componentRef = cr))
      guard(Util.stringStartsWith("$DER",ComponentReference.crefFirstIdent(cr)))
    then false;
    else true;
  end match;
end isInnerEqn;

protected function isObjective
  "Functions checks if an Equation has annotation Mayer or Lagrange."
  input BackendDAE.Equation eqn;
  output Boolean b;
algorithm
  b := match eqn
    local
      DAE.ComponentRef cr;
    case BackendDAE.EQUATION(exp = DAE.CREF(componentRef = cr))
      guard(Util.stringStartsWith("$OMC$object",ComponentReference.crefFirstIdent(cr)))
    then true;
  else false;
  end match;
end isObjective;

protected function isResidualEqn
  "Function checks if an equation is an residual equation."
  input BackendDAE.Equation eqn;
  output Boolean b;
algorithm
  b := match eqn
    local
      DAE.ComponentRef cr;
    case BackendDAE.EQUATION(exp = DAE.CREF(componentRef = cr))
      guard(ComponentReference.isSecondPartialDerivativeHessian(cr) or ComponentReference.isMayerOrLagrange(cr))
    then true;
    else false;
  end match;
end isResidualEqn;

protected function printHessian
  "Print the Hessian Matrices A,B and C"
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

protected function removeDerVars
  "Functions removes Variables of the type '$DER' from the current variable array"
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
          variables = BackendVariable.deleteVar(cr, variables);
        then "";
      else then "";
    end match;
  end for;
end removeDerVars;

protected function setHessianElement
  "Helper Function to set the the Hessian. The sum of Equations is the RHS. For the LHS a new Variable is created ('$Hess = lambda[i]*eq[i]+...')."
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

annotation(__OpenModelica_Interface="backend");
end SymbolicHessian;
