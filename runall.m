%% Testing with mua counts
clear all
close all
clc

load('Testpoissonworkspace.mat')
binDur = 0.010; % 5 milisecond bins
bins = -2 : binDur : 2 - binDur;
muaRate = histcounts( MUA_so, bins );

time = 1000;
[ qRev, sigRev ] = runpoissonPRE( muaRate, time );
[ ppm, e1, e2 ] = runpoisson( muaRate, qRev, sigRev, time );

%% Plot the estimate and its CIs.
figure
t2plot = bins( 1 : end - 1 );
bar( t2plot, muaRate, 'FaceColor', [ 0.6 0.6 0.6 ], 'EdgeColor', 'none' )
hold on
% plot( t2plot, e1, 'b' )
% plot( t2plot, ppm, 'r' )
% plot( t2plot, e2, 'b' )

%     x_vector = [options.x_axis', fliplr(options.x_axis')];
t2patch = [ t2plot'; fliplr( t2plot )' ];
y2patch = [ e1'; fliplr( e2 )' ];
hPatch = fill( t2patch, y2patch, 'r' );
set( hPatch, 'LineStyle', 'none' )
set( hPatch, 'FaceAlpha', 0.25 );
hold on;
plot( t2plot, ppm, 'color', 'r', ...
    'LineWidth', 1.5 );
hold off;


%% Original testing from anne's paper. Will reporoduce figure eight.
% I = [6     5     9     5     3     5     8     5     9     1     6     7     5     6     9 ...
%      1    11     8     5    14     7    19    12     2    22    27    37    35    38    24 ... 
%      29    23    25    25    17    32    35    24    30    27    33    27    35     9    38 ...
%      28    28    29    30    29    14    31    24    40    27];
% %make the spike data longer to see how long it takes to run..
% % I = [I I I I I I I I I I I I I I I I I I];
% % I = [I I I I I I I I I];
% % I = [I I];
% %time = 700;
% time = 1000; 
% %I=sum(bigI,1);
% % I=[zeros(100) 0 1 0 1 0 1 0 ];      
% %this run provides a better estimate of qstart for the start conditions
% [qstart, sigstart] = runpoissonPRE(I, time); 
% %estimates the count rate using the estimate of the initial condition from above
% [qnew, signewsq, ppm] = runpoisson(I, qstart, sigstart, time);     


