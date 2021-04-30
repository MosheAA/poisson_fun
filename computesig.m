function [newsigsq] = computesig(qnew, signewsq, a);

% Mstep for sigma^2_e assuming a fixed ic at x1
% doesn't use qnew(1) 


N = size(qnew, 2);

term1 = 2*(sum(qnew(3:end).^2) + sum(signewsq(3:end))); 

  a       = a(3:end);
  covcalc = signewsq(4:end).*a;   
  
term2  = -2*(sum(covcalc) + sum(qnew(3:N-1).*qnew(4:N)) ); 

term3 = - qnew(end)^2 - signewsq(end); 

term4  = qnew(2)^2 - 2*qnew(2)*qnew(3);    


newsigsq = (term1 + term2 + term3 + term4)/(N-2);

