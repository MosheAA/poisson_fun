function trialtotrial( stats )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code takes the output from runfiles.m
%and computes the within curve probability that
%one trial is different from another
%Anne Smith July 2004

%allps is a matrix that contains the final p-value comparisons
%between trial k and trial j

qnew = stats.qNew;
signewsq = stats.sigqNew;
a = stats.a;
muone = 0;
allps = [];

pcrit = 0.01;
numbermcs = 500;   %number of Monte Carlo samples

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%big loop through all the comparisons

for trial1 = 1:length(qnew)
    
    pvalue = NaN*zeros(1, length(qnew));
    
    for trial2 = trial1 + 1: length(qnew)
        
        shortvec = trial1: trial2-1;
        
        a1       = eye(size(a,1));
        
        % Appendix D in Smith, Stefani, Mogahaddam and Brown, 2004
        for ppp = trial2 - 1: -1: trial1
            a1  = mtimes(a(ppp), a1);
        end
        
        a1   = mtimes(a1, signewsq(trial2));
        covx = a1(1,1);
        covsave(trial1, trial2) = covx;
        
        mean1 = [ qnew(trial1) qnew(trial2)];
        cov1  = [ signewsq(trial1) covx;  ...
            covx signewsq(trial2) ];
        
        r         = mvnrnd(mean1, cov1, numbermcs);
        r1        = muone+ r(:,1);
        r2        = muone+ r(:,2);
        
        pp1       = exp(r1);
        pp2       = exp(r2);
        
        pvalue(trial2)  = length(find(r1>r2))/ numbermcs;
        
        %if(trial1 == 40 & trial2 == 41)
        %    cov1
        %    mean1
        %    figure(3); hist([pp1 pp2],100); length(find(pp1>pp2))/ numbermcs
        %end
        
    end
    
    allps = [allps; pvalue];
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

allps2 = allps;

for ii = 1:length(qnew)
    for jj = ii+1:length(qnew)
        allps2(jj,ii) = NaN; %  make the matrix lower triangular
    end
end

% figure(3);
% pno = 2;
% subplot(2, 1, pno)
imagesc(allps2,[0 1]); colormap('bone')

hold on;
[i,j] = find( allps2 < pcrit );

plot(j,i,'.r');
minj = min(j);
if(~isempty(minj))
    title(['Earliest trial signif above estimated start distribution '  num2str(minj(1)) ])
else
    title('No trials above start distribution')
end


[i,j] = find(allps2> (1-pcrit));

plot(j,i,'.b');

axis([0.5 length(qnew)-.5 0.5 length(qnew)-.5]);
line([0.5 length(qnew)-.5],[0.5 length(qnew)-.5]);

xlabel('Trial Number')
ylabel('Trial Number')


