function [ ppm, e1, e2, stats ] = runpoisson( I, qRev, sigNewRev, time)

% POISSON COUNT DATA 
% This function runs subroutines to analyze count data
% of spikes obtained from single neurons
% using a state-space paradigm.
% The state process, modulating the neural spiking activity, is a gaussian AR(1).
% Given the latent process, the space models is defined by the Poisson
% distribution.
% The spike times should be mesured in ms with a 1 ms resolution.
%
% Variables to be reset by user:
%
%        I:                       vector of count data (number of spikes/trial)
%        sige:                    SD of the state process
%        ccrt:                    Convergence criteria of the algorithm
%        qguess:                  Initial condition of the state process
%        time:                    Time in msec of experiment
%
% Other variables
%        q, s:                    Vectors of state process and its variance (forward estimate)
%        qnew, signewsq:          Vectors of state process and its variance (backward estimate)
%        newsigsq:                Estimate of state process variance from EM
%
%        (J. D. Scalon, 09/20/04)
%        (A. C. Smith,  01/03/05)
%
%-----------------------------------------------------------------
qguess              = qRev; %start mean
sigsqguess          = sigNewRev; %start variance
%--------------------------------------------------------------- 
  number_of_steps  = 2000; %maximum time steps for EM
  ccrt  = 1e-5;            %convergence criteria for the EM 
  startiter = 3;           %start recursive filter at trial 2 
  sige         = 0.50;   %EM estimates this value - this is a starting guess  
%-----------------------------------------------------------------------
% Loop for EM
for jk = 1:number_of_steps
%--------------------------------------------------------------
% Do forward filter 
   [q, s, qold, sold] = recpoisson(I, sige, qguess, sigsqguess, startiter);
%----------------------------------------------------------------
% Do backward filter
   [qnew, signewsq, a] = backpoisson(q, qold, s, sold, startiter);
%------------------------------------------------------------------
% Mstep      
   [newsigsq(jk)]      = computesig(qnew, signewsq, a);   
%-------------------------------------------------------------------------------
 % Set starting values for the next step jk 
   sige  = sqrt(newsigsq(jk)); %from EM      
%    if(sige < 1e-5)
%        fprintf(2, 'STOPPED because estimate got too small')
%        break
%    end      
%    fprintf(2, 'Final run iteration %d sige %f \n',jk, sige);
%------------------------------------------------------------------------------ 
%check for convergence
%    if(jk>1)
%       a1 = abs(newsigsq(jk) - newsigsq(jk-1));      
%      if(a1 < ccrt)
%          fprintf(2,'EM estimates converged in %d iterations \n',jk)
%          break
%      end
%    end
end
%-------------------------------------------------------------------------

% if(jk == number_of_steps)
%  fprintf(2,'EM failed to converge after %d steps; convergence criterion was %f \n', jk, ccrt)
% end
%-------------------------------------------------------------------
%confidence intervals for estimates
%take off the zeroth value as this is not estimated 
qnew(1) =[]; signewsq(1) = []; a(1) =[];
qnew(1) = qguess; 
[e1, ppm, e2] = getcls( qnew, signewsq, 1 );

stats.qNew = qnew;
stats.sigqnew = signewsq;
stats.a = a;
%-----------------------------------------------------------------------
%------------------------------------------------------------------------------
% Get firing rates (spikes/s or Hz)
%  Ir           = (I*1000)/time;
%  ppmr         = (ppm*1000)/time;
%  e1r          = (e1*1000)/time;
%  e2r          = (e2*1000)/time;
%---------------------------------------------------------------------------
% %plot up the figures
%  figure(1); subplot(211);
%   t = 1:1:size(I,2);
%   hold on;
%   %plot(t, Ir,'ko');
%   hold on;
%   plot(t, ppmr,'b-');
%   hold on;
%   plot(t, e1r,'r-');
%   hold on;
%   plot(t, e2r,'r-');
%   hold on;
%   xlabel('Trial Number');
%   ylabel('Spike Count or Rate (change this)'); box on
% figure(2)
%  bar(t, I, 0.01,'k');
% xlabel('Trial');
% ylabel('Counts (spikes/trial)');
% figure(3),
%  plot(newsigsq, 'b');
%  title('Convergence of sigma e');  
%---------------------------------------------------------
% save ppmr ppmr
% save e1r e1r
% save e2r e2r
% save qnew qnew
% save a a
%trialtotrial
   
   