function [ maxDeltaF, deltaOmegas] = script_positionDependence( L, diameter, Q, Vg, VgAC, Vbias, h, cRatio, dBzdx, Sz, TPrime_0, pos, TStepCount, fStepCount, zVec)
%POSDEPENDENCESCRIPT Summary of this function goes here
%   Detailed explanation goes here
maxDeltaF = zeros(1,max(size(pos)));
deltaOmegas = cell(1,max(size(pos)));

for j=1:max(size(pos))
    % Find closest point in zVec to position
    currPos = pos(j)*L-L/2;
    [ ~, pullIndex ] = min(abs(zVec-currPos));
    [ ~, ~, pdeltaOmega] = VecNumeric_pos2( L, diameter, Q, Vg, VgAC, Vbias, h, cRatio, dBzdx(pullIndex), Sz, TPrime_0, pos(j), TStepCount, fStepCount);
    maxDeltaF(j) = max(abs(pdeltaOmega));
    deltaOmegas{j} = pdeltaOmega;
end

