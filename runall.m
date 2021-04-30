%spike data
I = [6     5     9     5     3     5     8     5     9     1     6     7     5     6     9 ...
     1    11     8     5    14     7    19    12     2    22    27    37    35    38    24 ... 
     29    23    25    25    17    32    35    24    30    27    33    27    35     9    38 ...
     28    28    29    30    29    14    31    24    40    27];
%make the spike data longer to see how long it takes to run..
% I = [I I I I I I I I I I I I I I I I I I];
% I = [I I I I I I I I I];
% I = [I I];
%time = 700;
time = 1000; 
%I=sum(bigI,1);
% I=[zeros(100) 0 1 0 1 0 1 0 ];      
%this run provides a better estimate of qstart for the start conditions
[qstart, sigstart] = runpoissonPRE(I, time); 
%estimates the count rate using the estimate of the initial condition from above
[qnew, signewsq, ppm] = runpoisson(I, qstart, sigstart, time);     


