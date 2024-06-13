function [c, ceq] = constraint_phi(theta, data)
   % No equality constraints
   ceq = []; 
  
   % Define bounds for |phi|
    lower_bound = 0;  % Replace with your desired lower bound
    upper_bound = 1.4;  % Replace with your desired upper bound
    
    
   % calculate phi
   [~, phi] = neg_ll_diag(theta, data);
   
   % Inequality constraints to ensure lower_bound <= |phi| <= upper_bound
    c = [
        lower_bound - abs(phi);  % abs(phi) >= lower_bound
        abs(phi) - upper_bound    % abs(phi) <= upper_bound 
    ]; 
end