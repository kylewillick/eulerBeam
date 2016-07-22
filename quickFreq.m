function [ resFreq, step, soln ] = quickFreq( L, diameter, T, stepCount )
% Calculates resonant frequency using defined Tension

% To Do: 

% Load fixed CNT and SMM parameters
load FixedParameters.mat

% Calculated parameters
rOut = diameter/2;
rIn = rOut - wallThickness;
mCNT = rhoA * 0.735*pi*diameter*L;
% momentInertia = (pi/4)*(rOut.^4 - rIn.^4);
momentInertia = (pi/4)*(rOut.^4);
epsilon = 1e-12; % For regular use, escape condition has < 1Hz accuracy
% epsilon = 0; % For Debug: shuts off escape condition

xi = sqrt(T / (E * momentInertia));

% Begin Calculations
step = 1:stepCount;
if xi*L < 1
    freqStart = sqrt(E*momentInertia/(mCNT/L))*(22.38/L^2 + 0.28*xi^2); % Start from low T equation
else
    freqStart = sqrt(E*momentInertia/(mCNT/L))*(2*pi/L^2 + pi * xi /L); % Start from high T equation
end
threshold = epsilon*freqStart;
escapeIn = 9999999;
for k = step
    lambda = sqrt(mCNT/(L*E*momentInertia))*freqStart;
    yPlus = (L/sqrt(2))*sqrt(sqrt(xi^4 +4*lambda^2)+xi^2);
    yMinus = (L/sqrt(2))*sqrt(sqrt(xi^4 +4*lambda^2)-xi^2);
%     soln = cos(yMinus) - (yPlus^2 - yMinus^2)/(2*yPlus*yMinus) * sin(yMinus) %F
    soln = cosh(yPlus)*cos(yMinus) - (yPlus^2 - yMinus^2)/(2*yPlus*yMinus) * sinh(yPlus)*sin(yMinus);
    resFreq(k) = freqStart*(1 - (0.05-0.04*k/stepCount)*(soln-1)/exp(yPlus));
    freqStart = resFreq(k);
    % If difference in values within a threshold, run a few more loops then
    % exit
    if k > 2 && abs(resFreq(k)-resFreq(k-1)) < threshold && escapeIn > 10
        escapeIn = 3;
    end
    % The escape
    if escapeIn < 1
        break
    end
    escapeIn = escapeIn - 1;
end