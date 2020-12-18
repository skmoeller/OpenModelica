model MyBatchReactor
  Real x2(start = 0, fixed = true);
  Real x1(start = 1, fixed = true);
  input Real u(min=0, max = 5, nominal = 1.0, start = 1.0);
equation
  der(x1) = -(u+u^2/2)*x1;
  der(x2) = u*x1;
end MyBatchReactor;

optimization nmpcMyBatchReactorConstraints(objectiveIntegrand = -x2)
  extends MyBatchReactor;
  Real c1 = min(exp(2*time),c3)-0.5;
  Real c2 = 1+ exp(time) ;
  Real c3 = 2+ exp(time) ;
constraint
 u >= c1;
 u >= c2;
 u <= c3;
end nmpcMyBatchReactorConstraints;