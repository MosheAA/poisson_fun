
function [e1, ppm, e2] = getcls(qg, sg, conv_factor);
%
% This function gets the estimated mean firing rate and
% the 95% simulation bands using Monte Carlo simulation.
%
%------------------------------------------------
for i = 1:length(qg)

% pp = exp(  qg(i)  +  sqrt(sg(i))*randn(10000,1) );
 pp = exp(qg(i)  +  normrnd(0, sqrt(sg(i)), 10000, 1));

 pp = conv_factor*pp;
 
 psort = sort(pp);

 e1(i)  = psort(.05*10000);
 ppm(i) = psort(.50*10000);
 e2(i)  = psort(.95*10000);

end
%----------------------------------------------------
%----------------------------------------------------
