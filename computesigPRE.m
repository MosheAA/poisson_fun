function [newsigsq] = computesigPRE(qnew, signewsq, a);

% This function computes sigma^2_e assuming
% an unknown start distribution

N           = size(qnew, 2);

qnewt       = qnew(3:N);
qnewtm1     = qnew(2:N-1);
signewsqt   = signewsq(3:N);

term1 = 2*(sum(qnewt.^2) + sum(signewsqt));  %sum 2 to K

a       = a(2:end);
covcalc = signewsqt.*a;    

term2  = -2*(sum(covcalc) + sum(qnew(2:N-1).*qnew(3:N)) );  %sum 2 to K

term3  = qnew(2)^2 + 2*signewsq(2);
term4  = -qnew(end)^2 - signewsq(end); 
    
newsigsq = (term1 + term2 + term3 + term4)/(N-1);


