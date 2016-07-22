function [ omega0, omega0F, deltaOmega, tension, tensionF, xMax, xMaxF, deltaX0, VgDeltaX0 ] = mainEulerBeam( L, diameter, Q, Vg, VgAC, Vbias, h, cRatio, dBdz, Sz, T_0, pos, TStepCount, fStepCount)
% Calculate resonant frequency and frequency shifts of CNT beam with small
% applied force (by way of Sz change in dBdz)
% Outputs will always be arrays with dimensions as needed in order of input
% variables
%
% All units are SI (eg, L is in m)
%
% INPUT
%   L - length of suspended CNT section
%   diameter - diameter of CNT
%   Q - quality factor of mechanical resonance
%   Vg - DC gate voltage applied
%   VgAC - the magnitude of the AC signal applied to drive CNT
%   vBias - voltage that would be applied to source of CNT
%   h - separation from gate to CNT
%   cRatio - the ratio between capacitances (1/leverArm). C_total/C_g
%   dBdz - the magnitude of magnetic field gradient along perpendicular
%      direction
%   Sz - the change in spin of the attached SMM (force change is g mu_b Sz
%      dBdz)
%   T_0 - the residual (built in) tension scaled for device (eg, buckling
%      is at T_0 = -39.4)
%   pos - the position of the applied point force (from SMM) as ratio of L
%      (ie, SMM in middle is at pos = 0.5)
%   TStepCount - the maximum number of iterations in Tension solver
%   fStepCount - the maximum number of iterations in frequency solver
%
% OUTPUT
%   omega0 - resonance frequency of CNT with no point force
%   omega0F - resonance frequency when point force is applied
%   deltaOmega - omega0F - omega0
%   tension - the tension in the CNT with no point force
%   tensionF - tension with applied point force
%   xMax - maximum displacement along CNT (no point force)
%   xMaxF - + point force
%   deltaX0 & VgdeltaX0 - related to max displacement or oscillation
%   amplitude when driving, exact meaning unknown (calculated with point
%   force)

fprintf('---Starting eulerBeam Solver---\nCheck Matlab status bar for progress\n');

% Load fixed CNT and SMM parameters
load FixedParameters.mat

% Measure size of all inputs
LDim = max(size(L));
dDim = max(size(diameter));
QDim = max(size(Q));
VgDim = max(size(Vg));
ACDim = max(size(VgAC));
VbDim = max(size(Vbias));
hDim = max(size(h));
cDim = max(size(cRatio));
dBDim = max(size(dBdz));
SzDim = max(size(Sz));
T0Dim = max(size(T_0));
posDim = max(size(pos));

% Make a giant array to hold all solutions, this will be flattened later
tmp1 = zeros(LDim, dDim, QDim, VgDim, ACDim, VbDim, hDim, cDim, dBDim, SzDim, T0Dim, posDim); %Stores omega0
tmp2 = zeros(LDim, dDim, QDim, VgDim, ACDim, VbDim, hDim, cDim, dBDim, SzDim, T0Dim, posDim); %Stores omega0F
tmp3 = zeros(LDim, dDim, QDim, VgDim, ACDim, VbDim, hDim, cDim, dBDim, SzDim, T0Dim, posDim); %Stores tension
tmp4 = zeros(LDim, dDim, QDim, VgDim, ACDim, VbDim, hDim, cDim, dBDim, SzDim, T0Dim, posDim); %Stores tensionF
tmp5 = zeros(LDim, dDim, QDim, VgDim, ACDim, VbDim, hDim, cDim, dBDim, SzDim, T0Dim, posDim); %Stores xMax
tmp6 = zeros(LDim, dDim, QDim, VgDim, ACDim, VbDim, hDim, cDim, dBDim, SzDim, T0Dim, posDim); %Stores xMaxF
STRdeltaX0 = zeros(LDim, dDim, QDim, VgDim, ACDim, VbDim, hDim, cDim, dBDim, SzDim, T0Dim, posDim); %Stores deltaX0
STRVgDeltaX0 =  zeros(LDim, dDim, QDim, VgDim, ACDim, VbDim, hDim, cDim, dBDim, SzDim, T0Dim, posDim); %Stores VgDeltaX0

% Loop over all
for Lk = 1:LDim
    for dk = 1:dDim
        for Qk = 1:QDim
            for Vgk = 1:VgDim
                for ACk = 1:ACDim
                    for Vbk = 1:VbDim
                        for hk = 1:hDim
                            for ck = 1:cDim
                                for dBk = 1:dBDim
                                    for Szk = 1:SzDim
                                        for T0k = 1:T0Dim
                                            for posk = 1:posDim
        %                                     fprintf('Starting Calculations: L-%d/%d, d-%d/%d, Q-%d/%d, Vg-%d/%d, Vb-%d/%d, h-%d/%d, c-%d/%d, dBdz-%d/%d, Sz-%d/%d\n',Lk,LDim,dk,dDim,Qk,QDim,Vgk,VgDim,Vbk,VbDim,hk,hDim,ck,cDim,dBk,dBDim,Szk,SzDim);
                                                setDesktopStatus(sprintf('Starting Calculations: L-%d/%d, d-%d/%d, Q-%d/%d, Vg-%d/%d, Vb-%d/%d, h-%d/%d, c-%d/%d, dBdz-%d/%d, Sz-%d/%d, T_0-%d/%d, pos-%d/%d\n',Lk,LDim,dk,dDim,Qk,QDim,Vgk,VgDim,Vbk,VbDim,hk,hDim,ck,cDim,dBk,dBDim,Szk,SzDim,T0k,T0Dim,posk,posDim));
                                                [ T, tmpXmax, ~, ~,~,~,K_elec, F_mag ] = eulerTension( L(Lk) , diameter(dk), Q(Qk), Vg(Vgk), Vbias(Vbk), h(hk),cRatio(ck), dBdz(dBk), 0, T_0(T0k), pos(posk), TStepCount );
                                                if max(size(T)) == TStepCount
                                                    fprintf('Warning: Hit max step count in omega0 T calc at %d %d %d %d %d %d %d %d %d %d %d %d\n',Lk,dk,Qk,Vgk,ACk, Vbk, hk,ck,dBk,Szk,T0k,posk);
                                                    fprintf('Current Setting: L = %e, d = %e, Q = %e, Vg = %e, VgAC = %e, Vb = %e, h = %e, cRatio = %e, dBdz = %e, Sz = 0(forced), T_0 = %e, pos = %e\n', L(Lk) , diameter(dk), Q(Qk), Vg(Vgk), VgAC(ACk), Vbias(Vbk), h(hk),cRatio(ck), dBdz(dBk),T_0(T0k),pos(posk));
                                                end
                                                tmp3(Lk,dk,Qk,Vgk,ACk,Vbk,hk,ck,dBk,Szk,T0k,posk) = T(end);
                                                tmp5(Lk,dk,Qk,Vgk,ACk,Vbk,hk,ck,dBk,Szk,T0k,posk) = tmpXmax;
                                                omega = eulerFreq( L(Lk), diameter(dk), T(end), K_elec, F_mag, pos(posk), fStepCount );
                                                if max(size(omega)) == fStepCount && max(abs(omega(end-5:end-1)-omega(end))) > 50
                                                    fprintf('Warning: Hit max step count in omega0 freq calc and large deviations at %d %d %d %d %d %d %d %d %d %d %d %d\n',Lk,dk,Qk,Vgk,ACk, Vbk, hk,ck,dBk,Szk,T0k,posk);
                                                    fprintf('Current Setting: L = %e, d = %e, Q = %e, Vg = %e, VgAC = %e, Vb = %e, h = %e, cRatio = %e, dBdz = %e, Sz = 0(forced), T_0 = %e, pos = %e\n', L(Lk) , diameter(dk), Q(Qk), Vg(Vgk), VgAC(ACk), Vbias(Vbk), h(hk),cRatio(ck), dBdz(dBk),T_0(T0k),pos(posk));
                                                end
                                                tmp1(Lk,dk,Qk,Vgk,ACk,Vbk,hk,ck,dBk,Szk,T0k,posk) = omega(end);
                                                [ T, tmpXmaxF, ~, ~,~,~,K_elec, F_mag ] = eulerTension( L(Lk) , diameter(dk), Q(Qk), Vg(Vgk), Vbias(Vbk), h(hk),cRatio(ck), dBdz(dBk), Sz(Szk), T_0(T0k), pos(posk), TStepCount );
                                                if max(size(T)) == TStepCount
                                                    fprintf('Warning: Hit max step count in omega0F T calc at %d %d %d %d %d %d %d %d %d %d %d %d\n',Lk,dk,Qk,Vgk,ACk, Vbk, hk,ck,dBk,Szk,T0k,posk);
                                                    fprintf('Current Setting: L = %e, d = %e, Q = %e, Vg = %e, VgAC = %e, Vb = %e, h = %e, cRatio = %e, dBdz = %e, Sz = %e, T_0 = %e, pos = %e\n', L(Lk) , diameter(dk), Q(Qk), Vg(Vgk), VgAC(ACk), Vbias(Vbk), h(hk),cRatio(ck), dBdz(dBk),Sz(Szk),T_0(T0k),pos(posk));
                                                end
                                                tmp4(Lk,dk,Qk,Vgk,ACk,Vbk,hk,ck,dBk,Szk,T0k,posk) = T(end);
                                                tmp6(Lk,dk,Qk,Vgk,ACk,Vbk,hk,ck,dBk,Szk,T0k,posk) = tmpXmaxF;
                                                omega = eulerFreq( L(Lk), diameter(dk), T(end), K_elec, F_mag, pos(posk), fStepCount );
                                                if max(size(omega)) == fStepCount  && max(abs(omega(end-5:end-1)-omega(end))) > 50
                                                    fprintf('Warning: Hit max step count in omega0F freq calc and large deviations at %d %d %d %d %d %d %d %d %d %d %d %d\n',Lk,dk,Qk,Vgk,ACk, Vbk, hk,ck,dBk,Szk,T0k,posk);
                                                    fprintf('Current Setting: L = %e, d = %e, Q = %e, Vg = %e, VgAC = %e, Vb = %e, h = %e, cRatio = %e, dBdz = %e, Sz = %e, T_0 = %e, pos = %e\n', L(Lk) , diameter(dk), Q(Qk), Vg(Vgk), VgAC(ACk), Vbias(Vbk), h(hk),cRatio(ck), dBdz(dBk),Sz(Szk),T_0(T0k),pos(posk));
                                                end
                                                tmp2(Lk,dk,Qk,Vgk,ACk,Vbk,hk,ck,dBk,Szk,T0k,posk) = omega(end);

                                                % Calculate the oscillation amplitude from analytic estimates
                                                dCg = 2*pi*epsilon0*L(Lk)/((log(2*(h(hk)-tmpXmaxF)/(diameter(dk)/2)))^2*(h(hk)-tmpXmaxF));
                                                mCNT = rhoA*0.735*pi*diameter(dk)*L(Lk);
                                                STRdeltaX0(Lk,dk,Qk,Vgk,ACk,Vbk,hk,ck,dBk,Szk,T0k,posk) = (Q(Qk) * Vg(Vgk) * VgAC(ACk) * dCg)/(mCNT * omega(end)^2);
                                                STRVgDeltaX0(Lk,dk,Qk,Vgk,ACk,Vbk,hk,ck,dBk,Szk,T0k,posk) = Vg(Vgk)*STRdeltaX0(Lk,dk,Qk,Vgk,ACk,Vbk,hk,ck,dBk,Szk,T0k);
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

omega0 = squeeze(tmp1);
omega0F = squeeze(tmp2);
deltaOmega = omega0F - omega0;
tension = squeeze(tmp3);
tensionF = squeeze(tmp4);
xMax = squeeze(tmp5);
xMaxF = squeeze(tmp6);
deltaX0 = squeeze(STRdeltaX0);
VgDeltaX0 = squeeze(STRVgDeltaX0);
fprintf('Complete on %s\n',datestr(now));