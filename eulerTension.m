function [ T, maxx, step, x, z, dxdz, K_electric, F_mag ] = eulerTension( L, diameter, Q, Vg, Vbias, h,cRatio, dBdz, Sz, TPrime_0, zPrime_0, stepCount, varargin )
%Calculates tension in CNT using Euler Beam analtyic expression with
%constant T and recursively approximating T
% Note TPrime_0 is residual, built in CNT tension (due to fabrication), and
% is given such that T_0 = TPrime_0 * E*momentInertia/L^2

% To Do: 
%  - Fix K_elec component for high dangerZ
%   ie, for xi=7.95e9, L=1e-6,
%   (sinh(xi*L)/(cosh(xi*L)-1))*(cosh(xi*myz)-1)-sinh(xi*myz)+xi*myz-xi*myz^2/L
%   throws NaN

% Numeric Settings
spatialElements = 1e5;
epsilon = 1e-10;
xiL_UpperLimit = 25;
xiL_LowerLimit = 3e-2;

% Load fixed CNT and SMM parameters
load FixedParameters.mat

% Read the varargin
effectiveCapLengthRatio = 1;
if ~isempty(varargin)
    nVarIn = length(varargin{1});
    if rem(nVarIn,2) == 1
        error('Expects an even number of varargin, using identifier followed by value');
    else
        for k = 1 : nVarIn/2
            tag = varargin{1}{2*k-1};
            if strcmp(tag,'effCapLenRatio')
                effectiveCapLengthRatio = varargin{1}{2*k};
            else
                error('Unrecognized varargin tag');
            end
        end
    end
end

% Calculated parameters
rOut = diameter/2;
rIn = rOut - wallThickness;
mCNT = rhoA * 0.735*pi*diameter*L;
% A = pi *(rOut.^2 - rIn.^2);
A = pi * (rOut.^2);
% momentInertia = (pi/4)*(rOut.^4 - rIn.^4);
momentInertia = (pi/4)*(rOut.^4);
lengthDensity = pi*diameter*rhoA;
Cg = 4*pi*epsilon0*L*effectiveCapLengthRatio./(2*log(2*h./rOut));
C = cRatio * Cg;
CL = (C-Cg)/2;

% The position of the point force
a = zPrime_0 * L;

K_electric = (1./(4*pi*epsilon0*L.^2*h)).*(1/cRatio)^2.*(C.*Vg-CL*Vbias).^2; %F/L
F_mag = g * muB * dBdz * Sz;
T_0 = TPrime_0*E*momentInertia/L^2;

% Begin Calculations
% **** Matlab can't handle cosh(xi*z)-1 - sinh(xi*z) for large z, so
% trickery is needed [The error occurs around xi*z = 37]
step = 1:stepCount;

TInit = T_0 + (E*A/24 * (K_electric^2 * L^2 + 3*K_electric*L*F_mag + 3*F_mag^2))^(1/3);  %Use the high T limit as a starting guess, unless:
% TInit = vpa(TInit);
% If TInit would give xi < 1, use the low T limit as a starting guess
if abs(sqrt(TInit / (E * momentInertia))*L) < 1
%     fprintf('Here\n');
    TInit = T_0 + K_electric^2*L^6*A/(60480*E*momentInertia^2) + F_mag^2*L^4*A/(30720*E*momentInertia^2); %0 will give a divide by zero error, so just get close
end

threshold = abs(epsilon*TInit);
z = linspace(0,L,spatialElements);
dz = z(2)-z(1);
dxdz = zeros(1,spatialElements);
counter = 1;
escapeIn = 99999999;
aIndex = floor(spatialElements * zPrime_0);
halfIndex = floor(spatialElements/2); % Used in F != 0 profiles
quarterIndex = floor(spatialElements/4); % Used in F!=0 profile with dangerz
threeQuarterIndex = spatialElements - quarterIndex;
for j = step
    xi = sqrt(TInit / (E * momentInertia));
    % **** Matlab also can't handle the cosh()... expression for small
    % values (xi*L < 2e-2), in this case we can use the Taylor series
    % approximations to calculate x
    % NB: The error from matlab is an oscillatory error, in that the
    % resulting x follows the right profile, but has small rapid
    % oscillations
%     if abs(xi) * L < xiL_LowerLimit
    if 0
        fprintf('Area0\n');
        cosh_xiL_minusOne = ((xi*L)^2)/2 + ((xi*L)^4)/24 + ((xi*L)^6)/720 + ((xi*L)^8)/40320;
        sinh_xiL = (xi*L) + ((xi*L)^3)/6 + ((xi*L)^5)/120 + ((xi*L)^7)/5040;
        cosh_xiZ_minusOne = ((xi*z).^2)/2 + ((xi*z).^4)/24 + ((xi*z).^6)/720 + ((xi*z).^8)/40320;
        sinh_xiZ = (xi*z) + ((xi*z).^3)/6 + ((xi*z).^5)/120 + ((xi*z).^7)/5040;
        % Determine the x(z) solution for the current Tension
        x = zeros(1,spatialElements);
        % First, the DC gate effect
        x = (K_electric * L / (2 * TInit * xi))*( (sinh_xiL*cosh_xiZ_minusOne/cosh_xiL_minusOne) - sinh_xiZ + xi*z - xi*z.^2/L);
        % Now solve the point force part and sum
        x(1:halfIndex) = x(1:halfIndex) + (F_mag/(2*E*momentInertia*xi^3) * (tanh(xi*L/4)*(cosh(xi*z(1:halfIndex))-1)-sinh(xi*z(1:halfIndex))+xi*z(1:halfIndex)));
        x((halfIndex+1):end) = x((halfIndex+1):end) + (F_mag/(2*E*momentInertia*xi^3) * (tanh(xi*L/4)*(cosh(xi*(L-z((halfIndex+1):end)))-1)-sinh(xi*(L-z((halfIndex+1):end)))+xi*(L-z((halfIndex+1):end))));
    else % xi*L is not too small to calculate
        dangerZ = xiL_UpperLimit/abs(xi); % z values past dangerZ will incorrect calculate the cosh..
%         if dangerZ > L*(1-1/spatialElements)
        if 1
%             fprintf('Area1\n');
            % Determine the profile from the DC gate force
            x = (K_electric * L / (2 * TInit * xi))*( (sinh(xi*L)/(cosh(xi*L)-1)) * (cosh(xi*z)-1) - sinh(xi*z) + xi*z - xi * z.^2 / L);
%             x = zeros(1,spatialElements);
            % Now determine the point force correction and sum
            sinhka = sinh(xi*a);
            sinhkL = sinh(xi*L);
            sinhkLa = sinh(xi*(L-a));
            coshka = cosh(xi*a);
            coshkL = cosh(xi*L);
            coshkLa = cosh(xi*(L-a));
            FPrime = F_mag/(E*momentInertia);
            k = xi;
            sigma1 = FPrime*(sinhka-sinhkL+sinhkLa+a*k+k*(L-a)*coshkL-L*k*coshkLa);
            sigma2 = L*k*sinhkL - 2*coshkL + 2;
            sigma3 = FPrime*(coshkL - coshka + coshkLa - k*(L-a)*sinhkL - 1);
            c1 = sigma3/(k*sigma2);
            c2 = sigma1/(k*sigma2);
            c3 = -sigma3/(k^2*sigma2);
            c4 = -sigma1/(k^3*sigma2);
            z1 = z(1:aIndex);
            z2 = z(aIndex+1:end);
            x(1:aIndex) = x(1:aIndex) + c1*sinh(k*z1)/k^2 + c2*cosh(k*z1)/k^2 + c3*z1+c4;
            x(aIndex+1:end) = x(aIndex+1:end) + FPrime*sinh(k*(z2-a))/k^3 - FPrime*(z2-a)/k^2 + c1*sinh(k*z2)/k^2 + c2*cosh(k*z2)/k^2 + c3*z2+c4;
            
            
%             coshka*coshkLa
%             FPrime = F_mag / (E*momentInertia*xi^3);
%             k=xi;
%             denom = (-2*sinhka*sinhkLa-2*coshka*coshkLa+sinhka*coshkLa*k*L+2+sinhkLa*coshka*k*L);
%             A1 = -FPrime * (sinhka + sinhkLa + sinhka*sinhkLa*(k*(L-a)) - sinhka*coshkLa - sinhkLa*coshka+k*a+coshka*coshkLa*k*(L-a)-coshkLa*k*L)/denom;
%             B1 = FPrime * (sinhka*coshkLa*k*(L-a)+1+sinhkLa*coshka*k*(L-a) - sinhka*sinhkLa - coshka*coshkLa-coshkLa+coshka)/denom;
%             A2 = -FPrime * (coshka*k*L - k*a*coshka*coshkLa + sinhkLa*coshka-sinhka*sinhkLa*k*a + sinhka*coshkLa+k*a-sinhka-k*L-sinhkLa)/denom;
% %             B2 = FPrime * (-sinhka*coshkLa*k*a + sinhka*sinhkLa - sinhkLa*coshka*k*a + coshkLa*coshka - 1 - coshkLa+coshka)/denom;
%             FPrime = F_mag/(E*momentInertia);
%             k = xi;
%             denom = k^4 * (sinhka*coshkLa*L-2*sinhka*coshkLa*a+sinhkLa*coshka*L-2*sinhkLa*coshka*a);
%             A1 = FPrime*(sinhka+sinhkLa-k*a-k*a*coshka*coshkLa+2*coshkLa*k*a+k*L*coshka*coshkLa-coshkLa*k*L+sinhka*sinhkLa*k*(L-a)-sinhkLa*coshka-sinhka*coshkLa)/denom;
%             B1 = -FPrime*(sinhka*coshkLa*k*(L-a)-1+sinhkLa*coshka*k*(L-a)-sinhka*sinhkLa-coshkLa*coshka+coshkLa+coshka)/denom;
%             A2 = FPrime*(-sinhkLa*coshka-2*k*a*coshka+k*a*coshka*coshkLa+k*L*coshka+sinhkLa+k*a-k*L-sinhka*coshkLa+sinhka+sinhka*sinhkLa*k*a)/denom;
%             B2 = FPrime*(1+coshkLa*coshka-coshkLa-coshka-sinhka*coshkLa*k*a+sinhka*sinhkLa-sinhkLa*coshka*k*a)/denom;
%             x(1:aIndex) = x(1:aIndex) + A1*(cosh(xi*z(1:aIndex))-1)+B1*(sinh(xi*z(1:aIndex)) - xi*z(1:aIndex));
%             x(aIndex+1:end) = x(aIndex+1:end) + A2*(cosh(xi*(L-z(aIndex+1:end)))-1) + B2*(sinh(xi*(L-z(aIndex+1:end)))-xi*(L-z(aIndex+1:end)));
            
            %x(1:halfIndex) = x(1:halfIndex) + (F_mag/(2*E*momentInertia*xi^3) * (tanh(xi*L/4)*(cosh(xi*z(1:halfIndex))-1)-sinh(xi*z(1:halfIndex))+xi*z(1:halfIndex)));
            %x((halfIndex+1):end) = x((halfIndex+1):end) + (F_mag/(2*E*momentInertia*xi^3) * (tanh(xi*L/4)*(cosh(xi*(L-z((halfIndex+1):end)))-1)-sinh(xi*(L-z((halfIndex+1):end)))+xi*(L-z((halfIndex+1):end))));
        else
            fprintf('Error: entered uncoded region\n');
            [ ~, dangerIndex ] = min(abs(z - dangerZ));
            dangerIndex2 = spatialElements - dangerIndex;
            x = zeros(1,spatialElements);
            x(1:dangerIndex) = (K_electric * L / (2 * TInit * xi))*( (sinh(xi*L)/(cosh(xi*L)-1)) * (cosh(xi*z(1:dangerIndex))-1) - sinh(xi*z(1:dangerIndex)) + xi*z(1:dangerIndex) - xi * z(1:dangerIndex).^2 / L);
            x((dangerIndex+1):(dangerIndex2-1)) = (K_electric * L / (2 * TInit * xi))*( -1 + xi*z((dangerIndex+1):(dangerIndex2-1)) - xi * z((dangerIndex+1):(dangerIndex2-1)).^2 / L);
            % Assumes symmetry of the trouble expression
            x(dangerIndex2:end) = (K_electric * L / (2 * TInit * xi))*( (sinh(xi*L)/(cosh(xi*L)-1)) * (cosh(xi*z(dangerIndex+1:-1:1))-1) - sinh(xi*z(dangerIndex+1:-1:1)) + xi*z(dangerIndex2:end) - xi * z(dangerIndex2:end).^2 / L);
            if halfIndex < dangerIndex
    %             fprintf('Area2\n');
                x(1:halfIndex) = x(1:halfIndex) + (F_mag/(2*E*momentInertia*xi^3) * (tanh(xi*L/4)*(cosh(xi*z(1:halfIndex))-1)-sinh(xi*z(1:halfIndex))+xi*z(1:halfIndex)));
                x((halfIndex+1):end) = x((halfIndex+1):end) + (F_mag/(2*E*momentInertia*xi^3) * (tanh(xi*L/4)*(cosh(xi*(L-z((halfIndex+1):end)))-1)-sinh(xi*(L-z((halfIndex+1):end)))+xi*(L-z((halfIndex+1):end))));
            else
                % The tanh(xi*L/4) formula needs to be examined in quarters
                % if danger will be reached within a half
                if quarterIndex < dangerIndex
%                     fprintf('Area3\n');
                    % In this case splitting to quarters will avoid danger.
                    % The "cosh - sinh" term's second quarter is a 180rot image of first
                    % quarter, about (z(1/4),-1)
                    x(1:quarterIndex) = x(1:quarterIndex) + (F_mag/(2*TInit*xi))*(tanh(xi*L/4)*(cosh(xi*z(1:quarterIndex))-1)-sinh(xi*z(1:quarterIndex))+xi*z(1:quarterIndex));
                    x(quarterIndex+1:halfIndex) = x(quarterIndex+1:halfIndex) + (F_mag/(2*TInit*xi))*(-2-tanh(xi*L/4)*(cosh(xi*z(quarterIndex:-1:1))-1)+sinh(xi*z(quarterIndex:-1:1))+xi*z(quarterIndex+1:halfIndex));
                    x(halfIndex+1:threeQuarterIndex) = x(halfIndex+1:threeQuarterIndex) + (F_mag/(2*TInit*xi))*(-2 - tanh(xi*L/4)*(cosh(xi*(L-z(end:-1:threeQuarterIndex+1)))-1)+sinh(xi*(L-z(end:-1:threeQuarterIndex+1))) + xi*(L-z(halfIndex+1:threeQuarterIndex)));
                    x(threeQuarterIndex+1:end) = x(threeQuarterIndex+1:end) + (F_mag/(2*TInit*xi))*(tanh(xi*L/4)*(cosh(xi*(L-z(threeQuarterIndex+1:end)))-1)-sinh(xi*(L-z(threeQuarterIndex+1:end)))+xi*(L-z(threeQuarterIndex+1:end)));
                else
%                     fprintf('Area4\n');
                    x(1:dangerIndex) = x(1:dangerIndex) + (F_mag/(2*E*momentInertia*xi^3) * (tanh(xi*L/4)*(cosh(xi*z(1:dangerIndex))-1)-sinh(xi*z(1:dangerIndex))+xi*z(1:dangerIndex)));
                    x(dangerIndex+1:(halfIndex-dangerIndex)) = x(dangerIndex+1:(halfIndex-dangerIndex)) + (F_mag/(2*TInit*xi))*(-1+xi*z(dangerIndex+1:(halfIndex-dangerIndex)));
                    x((halfIndex-dangerIndex+1):halfIndex) = x((halfIndex-dangerIndex+1):halfIndex) + (F_mag/(2*TInit*xi))*(-2 - tanh(xi*L/4)*(cosh(xi*z(dangerIndex:-1:1))-1) + sinh(xi*z(dangerIndex:-1:1)) +xi*z((halfIndex-dangerIndex+1):halfIndex));
                    x((halfIndex+1):(halfIndex+dangerIndex)) = x((halfIndex+1):(halfIndex+dangerIndex)) + (F_mag/(2*TInit*xi))*(-2 - tanh(xi*L/4)*(cosh(xi*(L-z(end:-1:dangerIndex2+1)))-1) + sinh(xi*(L-z(end:-1:dangerIndex2+1))) + xi*(L-z((halfIndex+1):(halfIndex+dangerIndex))));
                    x((halfIndex+dangerIndex+1):dangerIndex2) = x((halfIndex+dangerIndex+1):dangerIndex2) + (F_mag/(2*TInit*xi))*(-1+xi*(L-z((halfIndex+dangerIndex+1):dangerIndex2)));
                    x(dangerIndex2+1:end) = x(dangerIndex2+1:end) + (F_mag/(2*E*momentInertia*xi^3) * (tanh(xi*L/4)*(cosh(xi*(L-z(dangerIndex2+1:end))) - 1) - sinh(xi*(L-z(dangerIndex2+1:end))) + xi*(L-z(dangerIndex2+1:end))));
                end
            end
        end
    end
    % Calculate derivatives
    dxdz(1) = (x(2)-x(1))/dz;
    for count = 2:(spatialElements-1)
        dxdz(count) = (x(count+1)-x(count-1))/(2*dz);
    end
    dxdz(spatialElements) = (x(spatialElements)-x(spatialElements-1))/dz;
    T(j) = T_0 + E*A/(2) * sum(dxdz.^2)/spatialElements;
    % **** The algorithm diverges if we try T(k) as the next guess, so
    % average the current guess and the correction
%     TInit = 0.6*TInit + 0.4*T(k);
    TInit = (0.6+ 0.37*j/stepCount)*TInit + (0.4-0.37*j/stepCount)*T(j);
%    fprintf('Diff at step %f: %e \n',k,T(k)-TInit);
    T(j) = TInit;
    % If within error is within a threshold, do *3* more loops and exit
    if j > 2 && abs(T(j)-T(j-1)) < threshold && escapeIn > 10
        escapeIn = 3;
%         break
    end
    % The escape
    if escapeIn < 1
        break
    end
    % Helpful trick: If the algorithm is oscillating about a point,
    % try to force the solution as the average of the last ~10, and see if
    % that works
    if counter > 25 && abs(T(j)-T(j-24)) < (1e4*threshold)
        counter = 1;
        TInit = mean(T(j-10:j));
        T(j) = TInit;
    end
    counter = counter + 1;
    escapeIn = escapeIn - 1;
end

maxx = max(abs(x));