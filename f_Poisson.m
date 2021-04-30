function [estimated] = f_Poisson(wave, bins)

%wave= x_1;
[qstart, sigstart] = runpoissonPRE(wave, bins);
if isnan(qstart) || isnan(sigstart)
    return, 
end
[~,~,estimated] = runpoisson(wave, qstart, sigstart, bins);

end