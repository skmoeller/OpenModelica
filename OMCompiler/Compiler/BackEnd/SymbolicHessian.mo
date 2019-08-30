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
  input list<Real> lambda                  "Lagrange factors, at the moment used as known.";
  output BackendDAE.BackendDAE outHessian "second derivates-> this is the hessian";
algorithm
  outHessian := matchcontinue(inBackendDAE)
  local
    BackendDAE.EqSystems eqs;
    BackendDAE.Shared shared;
    BackendDAE.SymbolicJacobians linearModelMatrixes;
    BackendDAE.SymbolicJacobian A;
    BackendDAE.BackendDAE HessA;
    DAE.FunctionTree funcs, functionTree;
    list< .DAE.Constraint> constraints;
    Option<list<DAE.ComponentRef>> lambdas;

  case(_) equation
    true = Flags.getConfigBool(Flags.GENERATE_SYMBOLIC_HESSIAN);
    BackendDAE.DAE(eqs=eqs,shared=shared) = inBackendDAE;
    (linearModelMatrixes, funcs) = SymbolicJacobian.createLinearModelMatrixes(inBackendDAE, Config.acceptOptimicaGrammar());

    //A::linearModelMatrixes = linearModelMatrixes; get first matrix
    //HessA = generateSymbolicHessianA(A);

    shared = BackendDAEUtil.setSharedSymJacs(shared, linearModelMatrixes);
    functionTree = BackendDAEUtil.getFunctions(shared);
    functionTree = DAE.AvlTreePathFunction.join(functionTree, funcs);
    shared = BackendDAEUtil.setSharedFunctionTree(shared, functionTree);
    outHessian = BackendDAE.DAE(eqs,shared);
  then outHessian;

  else inBackendDAE;
  end matchcontinue;
  print("\n\nOutput after linearization\n\n");
  BackendDump.printBackendDAE(outHessian);
  //outHessian:=transformJacobian(outHessian,lambda); //Umbauen des Gleichungssystems der jacobi matrix damit dann JAcobi erneut verwendet werden kann
  //outHessian:=SymbolicJacobian.symbolicJacobian(outHessian); //Erzeuge jetzt Hessematrix dies wird dann genauso zurueck gegeben!
end generateSymbolicHessian;

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
/*
protected function calculateSymbolicHessian
  "Function calculates the symbolic Hessian by using the algorithm for the jacobians."
 input BackendDAE.BackendDAE inBackendDAE "reducedDAE (variables and equations needed to calculate resVars)";
 input list<BackendDAE.Var> inVars        "independent vars";
 input BackendDAE.Variables inDiffedVars  "resVars";
 input BackendDAE.Variables inSeedVars    "Seed Variables";
 input BackendDAE.Variables inStateVars   "Stats";
 input BackendDAE.Variables inInputVars   "Input Variables";
 input BackendDAE.Variables inParamVars   "globalKnownVars";
 input String inMatrixName                "Name of the Matrix";
 output BackendDAE.BackendDAE outHessian  "The symbolic Hessian";
 output DAE.FunctionTree outFunctions;
algorithm
 (outHessian,outFunctions) := matchcontinue(inBackendDAE, inVars, inDiffedVars, inSeedVars, inStateVars, inInputVars, inParamVars, inMatrixName)
   local
     BackendDAE.BackendDAE bDAE;
     DAE.FunctionTree functions;
     list<DAE.ComponentRef> vars, comref_diffvars, comref_diffedvars;
     DAE.ComponentRef x;
     String dummyVarName;

     BackendDAE.Variables diffVarsArr;
     BackendDAE.Variables stateVars;
     BackendDAE.Variables inputVars;
     BackendDAE.Variables paramVars;
     BackendDAE.Variables diffedVars "resVars";
     BackendDAE.BackendDAE jacobian;

     // BackendDAE
     BackendDAE.Variables orderedVars, jacOrderedVars; // ordered Variables, only states and alg. vars
     BackendDAE.Variables globalKnownVars, jacKnownVars; // Known variables, i.e. constants and parameters
     BackendDAE.EquationArray orderedEqs, jacOrderedEqs; // ordered Equations
     BackendDAE.EquationArray removedEqs, jacRemovedEqs; // Removed equations a=b
     // end BackendDAE

     list<BackendDAE.Var> diffVars "independent vars", derivedVariables, diffedVarLst;
     list<BackendDAE.Equation> eqns, derivedEquations;

     list<list<BackendDAE.Equation>> derivedEquationslst;


     FCore.Cache cache;
     FCore.Graph graph;
     BackendDAE.Shared shared;

     String matrixName;
     array<Integer> ass2;
     list<Integer> assLst;

     BackendDAE.DifferentiateInputData diffData;

     BackendDAE.ExtraInfo ei;
     Integer size;

   case(BackendDAE.DAE(shared=BackendDAE.SHARED(cache=cache, graph=graph, info=ei, functionTree=functions)), {}, _, _, _, _, _, _) equation
     jacobian = BackendDAE.DAE( {BackendDAEUtil.createEqSystem(BackendVariable.emptyVars(), BackendEquation.emptyEqns())},
                                BackendDAEUtil.createEmptyShared(BackendDAE.JACOBIAN(), ei, cache, graph));
   then (jacobian, functions);

   case( BackendDAE.DAE( BackendDAE.EQSYSTEM(orderedVars=orderedVars, orderedEqs=orderedEqs, matching=BackendDAE.MATCHING(ass2=ass2))::{},
                        BackendDAE.SHARED(globalKnownVars=globalKnownVars, cache=cache,graph=graph, functionTree=functions, info=ei) ),
         diffVars, diffedVars, _, _, _, _, matrixName ) equation
     // Generate tmp variables
     dummyVarName = ("dummyVar" + matrixName);
     x = DAE.CREF_IDENT(dummyVarName,DAE.T_REAL_DEFAULT,{});

     // differentiate the equation system
     if Flags.isSet(Flags.JAC_DUMP2) then
       print("*** analytical Jacobians -> derived all algorithms time: " + realString(clock()) + "\n");
     end if;
     //Add prints here to see how the vars looking
     diffVarsArr = BackendVariable.listVar1(diffVars);
     comref_diffvars = List.map(diffVars, BackendVariable.varCref);
     diffData = BackendDAE.emptyInputData;
     diffData.independenentVars = SOME(diffVarsArr);
     diffData.dependenentVars = SOME(diffedVars);
     diffData.knownVars = SOME(globalKnownVars);
     diffData.allVars = SOME(orderedVars);
     diffData.diffCrefs = comref_diffvars;
     diffData.matrixName = SOME(matrixName);
     eqns = BackendEquation.equationList(orderedEqs);
     if Flags.isSet(Flags.JAC_DUMP2) then
       print("*** analytical Jacobians -> before derive all equation: " + realString(clock()) + "\n");
     end if;
     (derivedEquations, functions) = deriveAll(eqns, arrayList(ass2), x, diffData, functions);
     if Flags.isSet(Flags.JAC_DUMP2) then
       print("*** analytical Jacobians -> after derive all equation: " + realString(clock()) + "\n");
     end if;
     // replace all der(x), since ExpressionSolve can't handle der(x) proper
     derivedEquations = BackendEquation.replaceDerOpInEquationList(derivedEquations); //Hier wird etwas verändert, checken wie das gemacht wird, vielleicht ergibt sich daraus noch Anwendung; bei minimal beispiel passiert hier nichts!!!
     if Flags.isSet(Flags.JAC_DUMP2) then
       print("*** analytical Jacobians -> created all derived equation time: " + realString(clock()) + "\n");
     end if;

     // create BackendDAE.DAE with differentiated vars and equations

     // all variables for new equation system
     // d(ordered vars)/d(dummyVar)
     diffVars = BackendVariable.varList(orderedVars);
     derivedVariables = createAllDiffedVars(diffVars, x, diffedVars, matrixName);

     jacOrderedVars = BackendVariable.listVar1(derivedVariables);
     // known vars: all variable from original system + seed
     size = BackendVariable.varsSize(orderedVars) +
            BackendVariable.varsSize(globalKnownVars) +
            BackendVariable.varsSize(inSeedVars);
     jacKnownVars = BackendVariable.emptyVarsSized(size);
     jacKnownVars = BackendVariable.addVariables(orderedVars, jacKnownVars);
     jacKnownVars = BackendVariable.addVariables(globalKnownVars, jacKnownVars);
     jacKnownVars = BackendVariable.addVariables(inSeedVars, jacKnownVars);
     (jacKnownVars,_) = BackendVariable.traverseBackendDAEVarsWithUpdate(jacKnownVars, BackendVariable.setVarDirectionTpl, (DAE.INPUT()));
     jacOrderedEqs = BackendEquation.listEquation(derivedEquations);
    //todo print whats in the vars of jacobian!!!!!

     shared = BackendDAEUtil.createEmptyShared(BackendDAE.JACOBIAN(), ei, cache, graph);

     jacobian = BackendDAE.DAE( BackendDAEUtil.createEqSystem(jacOrderedVars, jacOrderedEqs)::{},
                                BackendDAEUtil.setSharedGlobalKnownVars(shared, jacKnownVars) );
   then (jacobian, functions);

   else
    equation
     Error.addInternalError("function generateSymbolicJacobian failed", sourceInfo());
   then fail();
 end matchcontinue;
end generateSymbolicHessian;
*/

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
