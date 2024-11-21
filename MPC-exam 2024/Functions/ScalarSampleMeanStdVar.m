function [xmean,s,xmeanp2s,xmeanm2s]=ScalarSampleMeanStdVar(x)

xmean = mean(x,2); % Gennemsnittet af dataene i 2. dimension (for hver række).
s = std(x,0,2); %Standardafvigelsen af dataene. Når du angiver 0, betyder 
% det, at du ønsker at beregne sample standard deviation.
% I dette tilfælde betyder 2, at standardafvigelsen skal beregnes langs den 
% anden dimension, dvs. for hver række i x

xmeanp2s = xmean + 2*s; %Værdien af gennemsnittet plus to gange standardafvigelsen.
xmeanm2s = xmean - 2*s; %Værdien af gennemsnittet minus to gange standardafvigelsen.
