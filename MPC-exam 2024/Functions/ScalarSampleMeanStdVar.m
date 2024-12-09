function [xmean, s, xmeanp2s, xmeanm2s] = ScalarSampleMeanStdVar(x)
% Computes the mean, standard deviation, and 2-sigma bounds for each row of x.

% Compute the row-wise mean
xmean = mean(x, 2);  

% Compute the row-wise standard deviation (sample standard deviation)
s = std(x, 0, 2);  

% Compute the upper and lower 2-sigma bounds
xmeanp2s = xmean + 2 * s;  
xmeanm2s = xmean - 2 * s;  

end
