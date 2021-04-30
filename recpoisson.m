function  [qhat, sigsq, qhatold, sigsqold] = recpoisson(I, sige, qguess, ...
                                                        sigsqguess, startiter)

% Implements the forward recursive filtering algorithm
% on the count data of spikes (I)
% Variables:
%        qhatold      one-step prediction
%        sigsqold     one-step prediction variance
%        qhat         posterior mode
%        sigsq        posterior variance
%
%        (J. D. Scalon - 08/13/2004)
%
%---------------------------------------------------
%set up some initial values

T                     = size(I,2); % Time

qhat(startiter-1)     = qguess;
sigsq(startiter-1)    = sigsqguess;  

number_fail = []; %number_fail saves the time steps if Newton method fails

%------------------------------------------------------
%-------------------------------------------------------

% Loop through all time

for t = startiter:T+1

   qhatold(t)  = qhat(t-1);
   sigsqold(t) = sigsq(t-1) + sige^2;
%----------------------------------------------------------------

% Calls newton solve to find solution to nonlinear posterior prediction estimate

   [qreturn, flagfail] = newtonsolve(I(t-1), qhatold(t), sigsqold(t));

   if flagfail > 0
      number_fail = [number_fail t];
   end
   
   qhat(t) = qreturn;
%-----------------------------------------------------

% Getting posterior variance
   denom       = -1/sigsqold(t) - exp(qhat(t));
   sigsq(t)    = -1/denom;
   
end

% if isempty(number_fail) < 1
%    fprintf(2,'Newton convergence failed at times %d \n', number_fail)
% end
