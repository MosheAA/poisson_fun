function [qnew, signewsq, a] = backpoisson(q, qold, sigsq, sigsqold, startiter);

% Given the sequence of posterior mean (q) and variance (sigsq) estimates,
% this function computes smooth estimates of mean and variance
% using the fixed-interval smoothing algorithm.
%
% (j.d.scalon 07/19/2004) (a.c.smith 11/10/04)

T = size(q,2);

qnew(T)     = q(T);
signewsq(T) = sigsq(T);

for i = T-1 :-1: startiter-1
   a(i)        = sigsq(i)/sigsqold(i+1);
   qnew(i)     = q(i) + a(i)*(qnew(i+1) - qold(i+1));
   signewsq(i) = sigsq(i) + a(i)*a(i)*(signewsq(i+1)-sigsqold(i+1));
end
