function [rho, phi] = acfpacf (y, maxOrd, signLvl, plotIt, maOrder, includeZeroLag)

if nargin<2
    maxOrd = 50;
end
if nargin<3
    signLvl = 0.05;
end
if nargin<4
    plotIt = 1;
end
if nargin<5
    maOrder = 0;
end
if nargin<6
    includeZeroLag = 1;
end

figure;

subplot(211);
rho = acf(y, maxOrd, signLvl, plotIt, maOrder, includeZeroLag);
title("acf");

subplot(212);
phi = pacf(y, maxOrd, signLvl, plotIt, includeZeroLag);
title("pacf");

% subplot(122);
% normplot(y);
sgtitle('ACF and PACF for the one step predictor error')