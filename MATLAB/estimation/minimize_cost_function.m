function [solution,fval,exitflag,output]=minimize_cost_function(theta0, data)
 
% summary of the optimization problem-general form
% fmincon attempts to solve problems of the form:
%     min J(x)  subject to:  Aineq*theta  <= Bineq, Aeq*theta  = Beq (linear constraints)
%      x                    C(x) <= 0, Ceq(x) = 0   (nonlinear constraints)
%                               lb <= theta <= ub        (bounds)
 
 
Aineq=[];
Bineq=[];
Aeq=[];
Beq=[];
lb = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1];  % Lower bounds
ub = [2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1];   % Upper bounds
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% keyboard
% algorithm can be 'interior-point', 'SQP','active set', and 'trust region reflective'
% set the options
options = optimoptions(@fmincon,'Algorithm','interior-point','Display','iter','UseParallel',true);
 
[solution,fval,exitflag,output] = fmincon(@(theta)cost_function(theta0, data),theta0,Aineq,Bineq,Aeq,Beq,lb,ub,[],options);
 
 
function J=cost_function(theta, data)
    J = neg_ll(theta, data);    
end
 
end

