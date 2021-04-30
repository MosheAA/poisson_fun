function [qstart, sigstart] = runpoissonPRE(I, time)

%This runs the EM algorithm on the reversed data and 
%estimates the start point for the next EM code

I = fliplr(I); %reverse the data

%---------------------------------------------------------------

  number_of_steps  = 1500;
  ccrt             = 1e-5;

%--------------------------------------------------------------------------
 % Provide initial condition about the mean of the state process and
 % parameter MU
  
  qguess     = 3; 
  
  startiter  = 2; %this is where the kalman filter starts
  		  
%-------------------------------------------------------------------------
% Make the first guess of the random walk variance		  

  sige         = 0.10; 
  sigsqguess   = sige^2;

%-----------------------------------------------------------------------
% Loop for EM

for jk = 1:number_of_steps

%forwards filter 

    [q, s, qold, sold] = recpoisson(I, sige, qguess, sigsqguess, startiter);

%backwards filter

   [qnew, signewsq, a] = backpoisson(q, qold, s, sold, startiter);

   
% Mstep   
    
   [newsigsq(jk)]      = computesigPRE(qnew, signewsq, a);
   qnewsave(jk)        = qnew(2); %store start estimates to check convergence
   
%-------------------------------------------------------------------------------
 % Set starting values for the next step jk
 
   sige       = sqrt(newsigsq(jk)); %from EM
   
   sigsqguess = sige^2;  %start for kalman filter at next iteration
   qguess     = qnew(2); %start for kalman filter at next iteration 
 
   if(sige < 1e-5)
       fprintf(2, 'STOPPED because sigme got to small')
       break
   end

  % fprintf(2, 'Pre-run iteration %d x0 %f sige %f \n',jk, qguess, sige);

%------------------------------------------------------------------------------ 
%check for convergence

  
   if(jk>1)
      a1 = abs(newsigsq(jk) - newsigsq(jk-1));
      a2 = abs(qnewsave(jk) - qnewsave(jk-1));
      
     if(a1 < ccrt & a2 < ccrt )
   %      fprintf(2,'EM estimates converged in %d iterations \n',jk)
         break
     end
   end

%-------------------------------------------------------------------------
end

qnew(1) =[]; signewsq(1) = []; a(1) =[];

[e1, ppm, e2] = getcls(qnew, signewsq, 1);


%-----------------------------------------------------------------------
%------------------------------------------------------------------------------
% Get firing rates (spikes/s or Hz)

 Ir           = (I*1000)/time;
 ppmr         = (ppm*1000)/time;
 e1r          = (e1*1000)/time;
 e2r          = (e2*1000)/time;

%------------------------------------------------------------------------------

%qstart is the new start value for the next code 
qstart   = qnew(end);
sigstart = signewsq(end);


if(jk == number_of_steps)
    qstart = NaN;
    sigstart = NaN;
% fprintf(2,'EM failed to converge after %d steps; convergence criterion was %f \n', jk, ccrt)
end