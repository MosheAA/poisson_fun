function [q, timefail] = newtonsolve(obs, qold, sigoldsq);
    							
% This function solves the posterior mode equation using Newtons method
% Variable  it(i) is the estimate of posterior mode
%
% (J. D. Scalon 07/27/2004)

it(1) = qold + sigoldsq*(-exp(qold) + obs);

for i = 1:200
   g         = it(i) - qold - sigoldsq*(-exp(it(i)) + obs);
   gprime    = 1 + sigoldsq*exp(it(i));
   it(i+1)   = it(i) - g/gprime;
   q         = it(i+1);

   if abs(q - it(i)) < 1e-12
      %fprintf(2,'cvrged (1e-12) in %d iters \n', i);
      timefail = 0; 
      return
   end
end

if(i >= 200) 
   %fprintf(2,'failed to converge')
   %abs(q-it(i))
   timefail = 1;
        q=it(i);
   return 
end
