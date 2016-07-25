function [ resFreq, step, soln, tmp, solverMatrix] = eulerFreq( L, diameter, T, K_elec, F_mag, pos, stepCount )
% Calculates resonant frequency using defined Tension
% This function uses the extended model from Poot/Witkamp, with T_ac

% To Do: Check the "error catching" equations
% Figure out why HPrime seems to be double the calculated curvature (could
% just be numeric error)
%   Figure out errorCatching functions for new Tac equations (if needed)
%   Check +- on omegaPrimeInit correction (ie scaling on soln)
%   Onset of deltaOmega errors occurs as low a kPrime=19 (checked for
%   Vg=0.28, L=1.6e-6)

% Load fixed CNT and SMM parameters
load FixedParameters.mat

% Calculated parameters
rOut = diameter/2;
mCNT = rhoA * 0.735*pi*diameter*L;
% momentInertia = (pi/4)*(rOut.^4 - rIn.^4);
momentInertia = (pi/4)*(rOut.^4);
% epsilon = 1e-10; % For regular use, escape condition has < 1Hz accuracy
epsilon = 0; % For Debug: shuts off escape condition

escapeIn = -9999; %Counter Holder

spatialElements = 1e6; % Number of points to use in integral calculations
z = linspace(0,1,spatialElements);

% Begin Calculations
% Start with basic analytic estimate based on tension limit
step = 1:stepCount;
freqs = quickFreq( L, diameter, T, stepCount );
freqStart = freqs(end); % Always start low, so increasing will find soln
threshold = epsilon*freqStart;
threshold = 1;
slowFactor = 1; % A scaling factor that will slow down corrections when the solution starts to oscillate
errorSlow = 1;
slowFactorCountdown = 2000; % Number of loops before slowFactor can be used, this is reset in the slowFactor region to slowing "period"
slowFactorPeriod = 2000;

% Scale frequency
nu = (1/L^2)*sqrt(E*momentInertia/(mCNT/L));
omegaPrimeInit = freqStart/nu;
TPrime = T*L^2/(E*momentInertia);
kPrime = sqrt(TPrime);
lDcPrime = K_elec*L^4/(E*momentInertia*rOut);
lMagPrime = F_mag*L^3/(E*momentInertia*rOut);
% lMagPrime=0; lDcPrime=0; %Debug, turns off Tac

% Store some repeatedly used calculations
a = pos;
k = kPrime;
sinhka = sinh(k*a);
sinhk = sinh(k);
sinhkLa = sinh(k*(1-a));
coshka = cosh(k*a);
coshk = cosh(k);
coshkLa = cosh(k*(1-a));
sigma1 = (coshk - coshka + coshkLa - k*(1-a)*sinhk - 1);
sigma2 = (sinhka-sinhk+sinhkLa+a*k+k*(1-a)*coshk-k*coshkLa);
sigma3 = k*sinhk - 2*coshk + 2;

% Compute the DC function derivatives needed later
d2udz0 = ((lDcPrime/(2*kPrime^2))*( (kPrime*sinh(kPrime))/(cosh(kPrime)-1) - 2) + lMagPrime*sigma2/(kPrime*sigma3));
d2udz1 = ((lDcPrime/(2*kPrime^2))*( (kPrime*sinh(kPrime))/(cosh(kPrime)-1) - 2) + lMagPrime/(kPrime*sigma3) * (sigma1*sinhk + sigma2*coshk) + lMagPrime/kPrime * sinhkLa);
d3udz0 = -lDcPrime/2 + lMagPrime*sigma1/sigma3;
d3udz1 = lDcPrime/2 + lMagPrime/sigma3 * (sigma1*coshk + sigma2*sinhk) + lMagPrime*coshkLa;
d3udzAdiscont = lMagPrime; % d3udzAplus-d3udzAminus
d5udzAdiscont = kPrime^2*lMagPrime;

for j = step
    % Temporary Solution to Accuracy issue, if kPrime is too high, just
    % give up and exit
    if kPrime > 19
        resFreq(1) = 0;
        break
    end
    kPlus = 1/2 * sqrt(-2*TPrime + 2*sqrt(TPrime^2 + 4*omegaPrimeInit^2));
    kMinus = 1/2 * sqrt(2*TPrime + 2*sqrt(TPrime^2 + 4*omegaPrimeInit^2));
    
    % Store some repeatedly used calculations
    cospa = cos(kPlus*a);
    sinpa = sin(kPlus*a);
    cospLa = cos(kPlus*(1-a));
    sinpLa = sin(kPlus*(1-a));
    coshma = cosh(kMinus*a);
    sinhma = sinh(kMinus*a);
    coshmLa = cosh(kMinus*(1-a));
    sinhmLa = sinh(kMinus*(1-a));
    
    % Calculate T_ac components
    TacC1 = TacC1_calc(kPrime,kPlus,lDcPrime,lMagPrime,a);
    TacC2 = TacC2_calc(kPrime,kPlus,lDcPrime,lMagPrime,a);
    TacC3 = TacC3_calc(kPrime,kMinus,lDcPrime,lMagPrime,a);
    TacC4 = TacC4_calc(kPrime,kMinus,lDcPrime,lMagPrime,a);
    TacB1 = TacB1_calc(kPrime,kPlus,lDcPrime,lMagPrime,a);
    TacB2 = TacB2_calc(kPrime,kPlus,lDcPrime,lMagPrime,a);
    TacB3 = TacB3_calc(kPrime,kMinus,lDcPrime,lMagPrime,a);
    TacB4 = TacB4_calc(kPrime,kMinus,lDcPrime,lMagPrime,a);
%     if kPrime < 0.05 % If kPrime is too low, we need to use VPA
    TacB5 = TacE_calc(kPrime,lDcPrime,lMagPrime,a);
%     TacB5 = double(TacB5);
%     else
%     TacB5 = -(lDcPrime^2 / (2*kPrime^4))*(csch(kPrime/2))^2*(kPrime^2+kPrime*sinh(kPrime)-4*cosh(kPrime)+4); % DC term
%     TacB5 = TacB5 + TacEmag_calc(pos,kPrime,lMagPrime, lDcPrime );
%     end
    
    % Build the equation matrix
    line1 = [ 1 0 1 0 0 0 0 0 d2udz0 ];
    line2 = [ 0 kPlus 0 kMinus 0 0 0 0 d3udz0 ];
    line3 = [ 0 0 0 0 1 0 1 0 d2udz1 ];
    line4 = [ 0 0 0 0 0 -kPlus 0 -kMinus d3udz1 ];
    line5 = [ cospa sinpa coshma sinhma -cospLa -sinpLa -coshmLa -sinhmLa 0 ];
    line6 = [ -kPlus*sinpa kPlus*cospa kMinus*sinhma kMinus*coshma -kPlus*sinpLa kPlus*cospLa kMinus*sinhmLa kMinus*coshmLa -d3udzAdiscont ];
    line7 = [ -kPlus^2*cospa -kPlus^2*sinpa kMinus^2*coshma kMinus^2*sinhma kPlus^2*cospLa kPlus^2*sinpLa -kMinus^2*coshmLa -kMinus^2*sinhmLa 0 ];
    line8 = [ kPlus^3*sinpa -kPlus^3*cospa kMinus^3*sinhma kMinus^3*coshma kPlus^3*sinpLa -kPlus^3*cospLa kMinus^3*sinhmLa kMinus^3*coshmLa -d5udzAdiscont ];
    line9 = [ TacC1 TacC2 TacC3 TacC4 TacB1 TacB2 TacB3 TacB4 (TacB5+omegaPrimeInit^2)];
    solverMatrix = [line1;line2;line3;line4;line5;line6;line7;line8;line9];
    soln = det(solverMatrix);
    
    % Store the first solution value (for adjusting convergence)
    if j == 1
        initSoln = soln;
    end

    % The easy escape condition, if we directly hit a solution, just end
    if abs(soln) < epsilon
        resFreq(j) = omegaPrimeInit*nu;
        break;
    end
    
    % Below commented was semi-working as of 02-09-13 morning
    % If the solution is oscillating about a point, slow it down
    if slowFactorCountdown < 1 && 0.95*mean(abs(resFreq(j-50:j-45) - resFreq(j-1))) < mean(abs(resFreq(j-25:j-20)-resFreq(j-1)))
        slowFactor = slowFactor * 1e-1;
        slowFactorCountdown = slowFactorPeriod;
%         if slowFactor < 1e-5
%             fprintf('Warning: slowFactor critically low in numericFreq\n');
%         end
    end
    % Update the estimate of frequency
	if abs(TPrime) > 12
        resFreq(j) = omegaPrimeInit*nu * (1+(0.05-0.03*j/stepCount)*5e0*slowFactor*errorSlow*soln/(abs(initSoln)^(12/11)));
    elseif abs(TPrime) > 5
        resFreq(j) = omegaPrimeInit*nu * (1+(0.05-0.03*j/stepCount)*5e-1*slowFactor*errorSlow*soln/(abs(initSoln)^(12/11))); % THIS WAS A GOOD ONE
    elseif abs(TPrime) > 1
        resFreq(j) = omegaPrimeInit*nu * (1+(0.05-0.03*j/stepCount)*5e-2*slowFactor*errorSlow*soln/(abs(initSoln)^(12/11))); % THIS WAS A GOOD ONE
    elseif abs(TPrime) > 1e-5
        resFreq(j) = omegaPrimeInit*nu * (1+(0.05-0.03*j/stepCount)*slowFactor*errorSlow*1e-3*soln/max(abs(initSoln),1));
    else
        resFreq(j) = omegaPrimeInit*nu * (1+(0.05-0.03*j/stepCount)*slowFactor*errorSlow*1e-3*soln/max(abs(initSoln),1));
    end

    % Check condition, if the resFreq calculation throws a wild value, go
    % back a few steps and slow down corrections
    if abs(log10((resFreq(j)/freqStart))) > 1
        if j < 4
            resFreq(j) = freqStart;
        else
            resFreq(j) = resFreq(j-3);
        end
        errorSlow = errorSlow*1e-1;
%         if errorSlow < 1e-5
%             fprintf('Warning: Error slow critically low in NumericFreq\n');
%         end
    end
    omegaPrimeInit = resFreq(j)/nu;

    % If difference in values within a threshold, run a few more loops then
    % exit (ie, the soft escape condition)
    if j > 5 && mean(abs(resFreq(j-4:j-1)-resFreq(j))) < threshold && escapeIn < -1
        escapeIn = 3;
    end
    % The escape
    if escapeIn < 1 && escapeIn > -1
        break
    end
    escapeIn = escapeIn - 1;
    
    slowFactorCountdown = slowFactorCountdown-1;
    
    % For debugging: keep track of "soln" evolution
    tmp(j) = soln;    
% %     lastError = errorVal;
end
end

function TacB4 = TacB4_calc(k,m,fDc,fMag,a)
sigma10 = k*sinh(k)-2*cosh(k)+2;
sigma9 = k^2 - m^2;
sigma8 = cosh(m*(a-1));
sigma7 = sinh(m*(a-1));
sigma6 = k^2*m*sigma10;
sigma5 = k^2*sigma9*sigma10;
sigma3 = k*(sinh(k)-sinh(a*k)*sigma8)+m*sigma7*cosh(a*k);
sigma4 = k*(cosh(k)-cosh(a*k)*sigma8)+m*sigma7*sinh(a*k);
sigma2 = sinh(k*(a-1));
sigma1 = cosh(k*(a-1));

DCPart = sigma3/(2*k^2*sigma9) - (a*m*sigma7-sigma8+1)/(k^2*m^2)+sigma7/(2*k^2*m)-(sinh(k)*sigma4)/(2*k^2*(cosh(k)-1)*sigma9);
magLine1 = (k*sigma2-m*sigma7)/(k^2*sigma9) - sigma7/(k^2*m) + sigma3/sigma5 + sigma7/sigma6 + (sigma1*sigma4)/(k*sigma9*sigma10)+(sinh(k)*sigma3)/(k*sigma9*sigma10) - sigma1*sigma3/sigma5 + sigma2*sigma4/sigma5+sigma7*(cosh(a*k)-cosh(k))/sigma6 +sigma7*sinh(k)/(k*m*sigma10) - a*sigma4/(k*sigma9*sigma10) - sigma7*sigma1/sigma6 + cosh(a*k)*sigma3/sigma5 - sinh(a*k)*sigma4/sigma5 - cosh(k)*sigma4/(k*sigma9*sigma10) - cosh(k)*sigma3/sigma5 + sinh(k)*sigma4/sigma5;
magLine2 = (a*cosh(k)*sigma4)/(k*sigma9*sigma10) - (a*sinh(k)*sigma3)/(k*sigma9*sigma10) - a*sigma7*sinh(k)/(k*m*sigma10);

TacB4 = fDc*4*m*DCPart + fMag*4*m*(magLine1+magLine2);
end

function TacB3 = TacB3_calc(k,m,fDc,fMag,a)
s10 = k*sinh(k)-2*cosh(k)+2;
s9 = k^2 - m^2;
s8 = sinh(m*(a-1));
s7 = cosh(m*(a-1));
s6 = k^2*m*s10;
s5 = k^2*s9*s10;
s4 = m*(cosh(k)-cosh(a*k)*s7)+k*s8*sinh(a*k);
s3 = m*(sinh(k)-sinh(a*k)*s7)+k*s8*cosh(a*k);
s2 = cosh(k*(a-1));
s1 = (sinh(m*(a-1)/2))^2;

DcPart = (s8 - m*(a*s7-1))/(k^2*m^2) + s1/(k^2*m) - s4/(2*k^2*s9) + sinh(k)*s3/(2*k^2*(cosh(k)-1)*s9);
magLine1 = 2*s1/(k^2*m) + (s4-s2*s4+sinh(k)*s3+sinh(k*(a-1))*s3+cosh(a*k)*s4-sinh(a*k)*s3-cosh(k)*s4)/s5 - m*(s2-s7)/(k^2*s9)+(-2*s1-s1*2*cosh(a*k)+s1*cosh(k)*2+2*s1*s2)/s6;
magLine2 = (s2*s3-a*s3-cosh(k)*s3+sinh(k)*s4+a*cosh(k)*s3-a*sinh(k)*s4)/(k*s9*s10) +(-2*s1*sinh(k)+a*2*s1*sinh(k))/(k*m*s10);

TacB3 = -fDc*4*m*DcPart + fMag*4*m*(magLine1+magLine2);
end

function TacB2 = TacB2_calc(k,p,fDc,fMag,a)
C = fDc;
M = fMag;
s24 = cos(p*(1-a));
s23 = sin(p*(1-a));
s22 = 8*C*k^2*s24;
s21 = 8*C*p^2*s24;
s20 = 8*C*a*k^2*p*s23;
s19 = 4*M*a*k*p^3*s23;
s18 = 8*C*a*p^3*s23;
s17 = 2*C*k*p^3*s23;
s16 = 4*C*k^2*p*s23;
s15 = 4*C*k*p^2*s24;
s14 = 4*M*k*p^2*s24;
s13 = 4*C*p^3*s23;
s12 = 4*M*k^2*p*s23;
s11 = 8*C*k^2;
s10 = 8*C*p^2;
s9 = 2*k^4*p;
s8 = k-2*a*k;
s7 = 4*M*k*p^2;
s6 = 2*C*k^2*p^2;
s5 = 4*M*k^2*p^2;
s4 = 2*k^2*p^3;
s3 = 4*M*a*k^2*p^2;
s2 = cosh(k*(1-a));
s1 = sinh(k*(1-a));

denom = (-s9-s4)*cosh(k)+(k^5*p+k^3*p^3)*sinh(k)+s9+s4;
Ccoshk = s22-s10-s6-s11+s21+s13-s18+s16-s12-s3-2*M*k^2*p^2*s24-s20;
Csinhk = 4*C*k^3*(1-s24)+8*C*k*p^2 + s7-s15+s14-s17-2*C*k^3*p*s23 - 2*M*k*p^3*s23+4*C*a*k*p^3*s23+4*C*a*k^3*p*s23+s19+4*M*a*k^3*p*s23;
Ccoshak = s5+s13-s12+2*C*k^2*p^2*s24+4*M*a*k^2*p^2*s24;
Csinhak = -s7-s15-s14-s17-s19;
line01 = s11+s10-s6-s5-s22-s21-s13-4*C*p^3*s2*s23+s18-s16-4*M*k*p^2*s1+s12+s3+4*M*k^2*p^2*s2*s24-2*M*k^2*p^2*cosh(s8)*s24-4*C*k*p^2*s24*s1+s20-4*M*k*p^2*s24*s1;
line2 = 4*M*k^2*p*s2*s23 + 2*C*k*p^3*s1*s23 + 4*M*k*p^3*s1*s23-2*M*k*p^3*sinh(s8)*s23+2*C*k^2*p^2*s2*s24-4*M*a*k*p^3*s1*s23 - 4*M*a*k^2*p^2*s2*s24;

TacB2 = (Ccoshk*cosh(k)+Csinhk*sinh(k)+Ccoshak*cosh(a*k)+Csinhak*sinh(a*k)+line01+line2)/denom;
end

function TacB1 = TacB1_calc(k,p,fDc,fMag,a)
C = fDc;
M = fMag;
s20 = sin(p*(1-a));
s19 = cos(p*(1-a));
s18 = 8*C*k^2*s20;
s17 = 8*C*p^2*s20;
s16 = 8*C*a*p^3*s19;
s15 = 2*C*k*p^3*s19;
s14 = 4*C*k^2*p*s19;
s13 = 4*C*k*p^2*s20;
s12 = 4*M*k*p^2*s20;
s11 = 8*C*a*k^2*p*s19;
s10 = 4*M*a*k*p^3*s19;
s9 = 4*C*p^3*s19;
s8 = 4*M*k^2*p*s19;
s7 = 2*k^4*p;
s6 = k-2*a*k;
s5 = 4*C*k^2*p;
s4 = 2*k^2*p^3;
s3 = 4*M*k^2*p;
s2 = sinh(k*(1-a));
s1 = cosh(k*(1-a));

denom = (-s7-s4)*cosh(k)+(k^5*p+k^3*p^3)*sinh(k)+s7+s4;
Ccoshk = s9+s5-s18-s17+s3-s16+s14-s8+2*M*k^2*p^2*s20-s11;
Csinhk = 4*C*k^3*s20 - 2*C*k^3*p-s15-2*C*k^3*p*s19-4*M*a*k^3*p-2*M*k*p^3*s19+s13-s12+4*C*a*k*p^3*s19+4*C*a*k^3*p*s19+s10+4*M*a*k^3*p*s19;
Ccoshak = s9 + s3 -s8-2*C*k^2*p^2*s20-4*M*a*k^2*p^2*s20;
Csinhak = s13-s15+s12-s10;
line0 = s18-s5-s9+s17-s3;
line1 = -4*C*p^3*s1*s19 + s16 - s14 - 4*M*k^2*p*s1+s8-2*C*k^2*p^2*s1*s20-4*M*k^2*p^2*s1*s20+2*M*k^2*p^2*cosh(s6)*s20+s11+4*M*k^2*p*s1*s19+2*C*k*p^3*s19*s2+4*M*k*p^3*s19*s2;
line2 = -2*M*k*p^3*s19*sinh(s6)+4*C*k*p^2*s2*s20 + 4*M*k*p^2*s2*s20 - 4*M*a*k*p^3*s19*s2+4*M*a*k^2*p^2*s1*s20;

TacB1 = (Ccoshk*cosh(k)+Csinhk*sinh(k)+Ccoshak*cosh(a*k)+Csinhak*sinh(a*k)+line0+line1+line2)/denom;
end

function TacC4 = TacC4_calc(k,m,fDc,fMag,a)
s7 = k*sinh(k)-2*cosh(k)+2;
s6 = k^2-m^2;
s5 = k^2*m*s7;
s4 = k^2*s6*s7;
s3 = cosh(k*(a-1));
s2 = k*cosh(a*m)*sinh(a*k)-m*cosh(a*k)*sinh(a*m);
s1 = k*(cosh(a*k)*cosh(a*m)-1)-m*sinh(a*k)*sinh(a*m);

DcPart = (a*m*sinh(a*m)-cosh(a*m)+1)/(k^2*m^2) + s2/(2*k^2*s6) - sinh(a*m)/(2*k^2*m)-sinh(k)*s1/(2*k^2*(cosh(k)-1)*s6);
magLine1 = (-sinh(a*m)+sinh(a*m)*s3-cosh(a*k)*sinh(a*m)+sinh(a*m)*cosh(k))/s5 + (s2-sinh(a*k)*s1+cosh(a*k)*s2+sinh(k)*s1+sinh(k*(a-1))*s1-cosh(k)*s2-s2*s3)/s4;
magLine2 = (-sinh(a*m)*sinh(k)+a*sinh(a*m)*sinh(k))/(k*m*s7) + (-cosh(k)*s1 +s3*s1+sinh(k)*s2-a*s1-a*sinh(k)*s2+a*cosh(k)*s1)/(k*s6*s7);

TacC4 = fDc*(-4*m)*DcPart + fMag*(-4*m)*(magLine1+magLine2);
end

function TacC3 = TacC3_calc(k,m,fDc,fMag,a)
s8 = k*sinh(k)-2*cosh(k)+2;
s7 = k^2-m^2;
s6 = k^2*m*s8;
s5 = k^2*s7*s8;
s4 = cosh(k*(a-1));
s3 = (sinh(a*m/2))^2;
s2 = k*cosh(a*k)*sinh(a*m)-m*cosh(a*m)*sinh(a*k);
s1 = m*(cosh(a*k)*cosh(a*m)-1)-k*sinh(a*k)*sinh(a*m);

DcPart = (sinh(a*m)-a*m*cosh(a*m))/(k^2*m^2) + s1/(2*k^2*s7) + s3/(k^2*m) + sinh(k)*s2/(2*k^2*(cosh(k)-1)*s7);
magLine1 = (-s1-cosh(a*k)*s1+cosh(k)*s1-sinh(a*k)*s2+s4*s1+sinh(k)*s2+sinh(k*(a-1))*s2)/s5 + (-s3*2-cosh(a*k)*2*s3+s3*cosh(k)*2+s3*s4*2)/s6;
magLine2 = (-a*s2-sinh(k)*s1-cosh(k)*s2+s4*s2+a*sinh(k)*s1+a*cosh(k)*s2)/(k*s7*s8) + (-s3*2*sinh(k)+a*s3*2*sinh(k))/(k*m*s8);

TacC3 = fDc*4*m*DcPart + fMag*(-4*m)*(magLine1+magLine2);
end

function TacC2 = TacC2_calc(k,p,fDc,fMag,a)
C = fDc;
M = fMag;
s21 = sinh(k*(1-a));
s20 = cosh(k*(1-a));
s19 = 4*M*k*p^2*s21;
s18 = 4*M*k^2*p^2*s20;
s17 = 4*C*k^3;
s16 = 2*k^4*p;
s15 = k-2*a*k;
s14 = 8*C*a*p^3;
s13 = 2*C*k*p^3;
s12 = 4*C*k*p^2;
s11 = 4*C*k^2*p;
s10 = 8*C*a*k^2*p;
s9 = 4*M*a*k*p^3;
s8 = 2*k^2*p^3;
s7 = 4*C*p^3;
s6 = 4*M*k^2*p;
s5 = 8*C*k^2;
s4 = 8*C*p^2;
s3 = 4*M*k*p^2;
s2 = 2*C*k^2*p^2;
s1 = 4*M*a*k^2*p^2;

line1 = (s1+s2)*cos(a*p)*cosh(a*k)+(s3-s12-s17)*cos(a*p)*sinh(k)+(-s12-s3)*cos(a*p)*sinh(a*k)+(-2*M*k^2*p^2+s5+s4)*cos(a*p)*cosh(k);
line2 = (2*C*k^2*p^2*s20 - s4 - 4*C*k*p^2*s21-s19-s5+s18-2*M*k^2*p^2*cosh(s15)-4*M*a*k^2*p^2*s20)*cos(a*p)+(s6-s7)*sin(a*p)*cosh(a*k);
line3 = (s13+2*C*k^3*p+2*M*k*p^3+4*M*k^3*p-4*C*a*k*p^3-4*C*a*k^3*p-s9-4*M*a*k^3*p)*sin(a*p)*sinh(k)+(s13+s9)*sin(a*p)*sinh(a*k);
line4 = (s14-s7-s11-s6+s10)*sin(a*p)*cosh(k)+(s7+4*C*p^3*s20-s14+s11+s6-s10-4*M*k^2*p*s20-2*C*k*p^3*s21-4*M*k*p^3*s21+2*M*k*p^3*sinh(s15)+4*M*a*k*p^3*s21)*sin(a*p);
line5 = (s17+8*C*k*p^2+s3)*sinh(k)-s3*sinh(a*k)+(s1-s4-s2-4*M*k^2*p^2-s5)*cosh(k)+s5+s4-s2-s19+s18-s1;
denom = (k^5*p+k^3*p^3)*sinh(k)+(-s16-s8)*cosh(k)+s16+s8;

TacC2 = (line1+line2+line3+line4+line5)/denom;
end

function TacC1 = TacC1_calc(k,p,fDc,fMag,a)
s14 = sinh(k*(a-1));
s13 = 2*k*p^2*sin(a*p)*s14;
s12 = k*(2*a-1);
s11 = 2*k^2*p;
s10 = 2*k^2*p*cosh(k);
s9 = 2*k^2*p*cos(a*p);
s8 = k*p^3*cos(a*p)*sinh(k);
s7 = 2*k^2*p*cos(a*p)*cosh(k);
s6 = 2*k*p^2*sin(a*p)*sinh(k);
s5 = 2*a*k*p^3*cos(a*p)*sinh(k);
s4 = 2*a*k^3*p*cos(a*p)*sinh(k);
s3 = 2*k*p^2*sinh(a*k)*sin(a*p);
s2 = k^2*p*(k^2+p^2)*(k*sinh(k)-2*cosh(k)+2);
s1 = cosh(k*(a-1));

dcLine1 = 2*p^3*cos(a*p)-s11+4*k^2*sin(a*p)+4*p^2*sin(a*p)-4*a*p^3*cos(a*p)+s9-2*p^3*cosh(a*k)*cos(a*p)+s10-k^3*p*sinh(k)-2*p^3*cos(a*p)*cosh(k)-4*k^2*sin(a*p)*cosh(k)-4*p^2*sin(a*p)*cosh(k)+2*k^3*sin(a*p)*sinh(k)+2*p^3*cos(a*p)*s1;
dcLine2 = s8+k^3*p*cos(a*p)*sinh(k)+s6+k*p^3*cos(a*p)*s14-k^2*p^2*cosh(a*k)*sin(a*p)-s13-k^2*p^2*sin(a*p)*s1-4*a*k^2*p*cos(a*p)+k*p^3*cos(a*p)*sinh(a*k)+s3+4*a*p^3*cos(a*p)*cosh(k)-s7+4*a*k^2*p*cos(a*p)*cosh(k)-s5-s4;
magLine1 = s11+2*k^2*p*cosh(a*k)-s9-s10+2*k^3*p*sinh(k)-2*k^2*p*s1-s8-2*k^3*p*cos(a*p)*sinh(k)+2*k^2*p*cos(a*p)*s1+s6-2*k*p^3*cos(a*p)*s14+s13-k^2*p^2*sin(a*p)*cosh(k)+k*p^3*sinh(s12)*cos(a*p);
magLine2 = 2*k^2*p^2*sin(a*p)*s1-2*k^2*p*cosh(a*k)*cos(a*p)-2*a*k^3*p*sinh(k)-s3-k^2*p^2*cosh(s12)*sin(a*p)+s7-2*a*k*p^3*cos(a*p)*sinh(a*k)+s5+s4+2*a*k*p^3*cos(a*p)*s14+2*a*k^2*p^2*cosh(a*k)*sin(a*p)-2*a*k^2*p^2*sin(a*p)*s1;

TacC1 = fDc*2*(dcLine1+dcLine2)/s2 - fMag*2*(magLine1+magLine2)/s2;
end

function TacE = TacE_calc(k,fDc,fMag,a)
%C = vpa(fDc);
C = fDc;
M = fMag;
%M = vpa(fMag);
%k = vpa(k);
%a = vpa(a);
TacE = (-(32*cosh(k) - 8*cosh(2*k) + 10*k*sinh(2*k) - 3*k^3*sinh(k) - 4*k^2*cosh(2*k) + (k^3*sinh(2*k))/2 - 20*k*sinh(k) + k^4 + 4*k^2*cosh(k) + k^4*cosh(k) - 24)/(k^4*(cosh(k) - 1)*(4*cosh(k) - 4*k*sinh(k) + k^2 + k^2*cosh(k) - 4)))*C^2 + (-(4*k*sinh(2*k) - 6*k^3*sinh(k) - 8*k^2*cosh(k*(a - 1)) + 4*k^2*cosh(k*(a - 2)) - 6*k^3*sinh(k*(a - 1)) + 2*k^3*sinh(k*(a + 1)) + 2*k^3*sinh(k*(a - 2)) - (k^3*sinh(2*k*(a - 1)))/2 - 4*k^2*cosh(2*k) + (3*k^3*sinh(2*k))/2 + 12*k*sinh(a*k) + 2*k^3*sinh(k*(2*a - 1)) - 4*a*k^4 - 8*k*sinh(k) + 4*k^2*cosh(a*k) - 4*k^4*cosh(a*k) + k^4*cosh(2*a*k) - 12*k*sinh(k*(a - 1)) - 4*k*sinh(k*(a + 1)) + 4*k*sinh(k*(a - 2)) + 2*k^3*sinh(a*k) - (3*k^3*sinh(2*a*k))/2 - 4*k^2 + 3*k^4 + 4*a^2*k^4 + 8*k^2*cosh(k) - 12*a*k^2*cosh(a*k) + 6*a*k^4*cosh(a*k) - a*k^4*cosh(2*a*k) + 2*a^2*k^3*sinh(2*k) - 8*a*k^3*sinh(a*k) + 2*a*k^3*sinh(2*a*k) + 4*a*k^4*cosh(k) + 4*a*k^3*sinh(k) + 12*a*k^2*cosh(k*(a - 1)) + 4*a*k^2*cosh(k*(a + 1)) - 4*a*k^2*cosh(k*(a - 2)) - 2*a*k^4*cosh(k*(a - 1)) - 2*a*k^4*cosh(k*(a + 1)) - 2*a*k^4*cosh(k*(a - 2)) + a*k^4*cosh(2*k*(a - 1)) - 2*a^2*k^4*cosh(a*k) + 16*a*k^3*sinh(k*(a - 1)) - 8*a*k^3*sinh(k*(a - 2)) + 2*a*k^3*sinh(2*k*(a - 1)) + 12*a^2*k^3*sinh(a*k) - 4*a^2*k^4*cosh(k) - 2*a*k^3*sinh(2*k) - 4*a^2*k^3*sinh(k) - 2*a^2*k^4*cosh(k*(a - 1)) + 2*a^2*k^4*cosh(k*(a + 1)) + 2*a^2*k^4*cosh(k*(a - 2)) - 4*a*k^3*sinh(k*(2*a - 1)) - 12*a^2*k^3*sinh(k*(a - 1)) - 4*a^2*k^3*sinh(k*(a + 1)) + 4*a^2*k^3*sinh(k*(a - 2)))/(k^4*(cosh(k) - 1)*(4*cosh(k) - 4*k*sinh(k) + k^2 + k^2*cosh(k) - 4)))*M^2 + (4*C*M*k^2 - 2*C*M*k^4 + 8*C*M*k*sinh(k) - 4*C*M*k^2*cosh(a*k) + 2*C*M*k^4*cosh(a*k) + 12*C*M*k*sinh(k*(a - 1)) + 4*C*M*k*sinh(k*(a + 1)) - 4*C*M*k*sinh(k*(a - 2)) + C*M*k^3*sinh(a*k) - 8*C*M*k^2*cosh(k) - 2*C*M*k^4*cosh(k) - 4*C*M*k*sinh(2*k) + 6*C*M*k^3*sinh(k) + 8*C*M*k^2*cosh(k*(a - 1)) - 4*C*M*k^2*cosh(k*(a - 2)) + C*M*k^4*cosh(k*(a - 1)) + C*M*k^4*cosh(k*(a + 1)) + 3*C*M*k^3*sinh(k*(a - 1)) - 3*C*M*k^3*sinh(k*(a + 1)) - C*M*k^3*sinh(k*(a - 2)) + 4*C*M*k^2*cosh(2*k) - C*M*k^3*sinh(2*k) - 12*C*M*k*sinh(a*k) - 12*C*M*a*k^2*cosh(k*(a - 1)) - 4*C*M*a*k^2*cosh(k*(a + 1)) + 4*C*M*a*k^2*cosh(k*(a - 2)) + C*M*a*k^4*cosh(k*(a - 1)) - C*M*a*k^4*cosh(k*(a + 1)) + C*M*a*k^4*cosh(k*(a - 2)) - 4*C*M*a*k^3*sinh(k*(a - 1)) + 4*C*M*a*k^3*sinh(k*(a + 1)) + 4*C*M*a*k^3*sinh(k*(a - 2)) + 12*C*M*a*k^2*cosh(a*k) - C*M*a*k^4*cosh(a*k) - 4*C*M*a*k^3*sinh(a*k))/(k^4*(cosh(k) - 1)*(4*cosh(k) - 4*k*sinh(k) + k^2 + k^2*cosh(k) - 4));
end

function TacC1 = TacC1_calc_old(k,p,fDc,fMag,a)
C = fDc;
M = fMag;
s20 = cosh(k*(1-a));
s19 = 4*M*k^2*p*s20;
s18 = 8*C*k^2;
s17 = 8*C*p^2;
s16 = 2*k^4*p;
s15 = k-2*a*k;
s14 = 8*C*a*p^3;
s13 = 2*C*k*p^3;
s12 = 2*C*k^3*p;
s11 = 4*C*k*p^2;
s10 = 4*M*k*p^2;
s9 = 4*M*k^3*p;
s8 = 8*C*a*k^2*p;
s7 = 4*M*a*k*p^3;
s6 = 4*M*a*k^3*p;
s5 = 2*k^2*p^3;
s4 = 4*C*p^3;
s3 = 4*C*k^2*p;
s2 = sinh(k*(1-a));
s1 = 4*M*k^2*p;

denom = (k^5*p+k^3*p^3)*sinh(k)+(-s16-s5)*cosh(k)+s16+s5;
line1 = (s1-s4)*cos(a*p)*cosh(a*k)+(s13+s12+2*M*k*p^3+s9-4*C*a*k*p^3-4*C*a*k^3*p-s7-s6)*cos(a*p)*sinh(k)+(s13+s7)*cos(a*p)*sinh(a*k);
line2 = (s14-s4-s3-s1-s8)*cos(a*p)*cosh(k)+(s4+4*C*p^3*s20-s14+s3+s1-s8-s19-2*C*k*p^3*s2-4*M*k*p^3*s2+2*M*k*p^3*sinh(s15)+4*M*a*k*p^3*s2)*cos(a*p);
line3 = (-2*C*k^2*p^2-4*M*a*k^2*p^2)*sin(a*p)*cosh(a*k)+(4*C*k^3+s11-s10)*sin(a*p)*sinh(k)+(s11+s10)*sin(a*p)*sinh(a*k)+(2*M*k^2*p^2-s18-s17)*sin(a*p)*cosh(k);
line4 = (s18+s17+4*C*k*p^2*s2+4*M*k*p^2*s2-2*C*k^2*p^2*s20-4*M*k^2*p^2*s20+2*M*k^2*p^2*cosh(s15)+4*M*a*k^2*p^2*s20)*sin(a*p)-s1*cosh(a*k)+(s6-s9-s12)*sinh(k);
line5 = (s3+s1)*cosh(k)+s19-s1-s3;

TacC1 = (line1+line2+line3+line4+line5)/denom;
end

function [ TacEmag ] = TacEmag_calc( a, k, lMagPrime, lDcPrime )
%TACEMAG_calc Uses results from OffsetPointForceSolve MuPad to calculate
%magnetic force contribution to TacE, with general position of force
% a is position (in primed coordinates)
% k is kPrime
    sigma1 = sinh(k*(a-1));
    sigma2 = sinh(k*(a-2));
    sigma3 = sinh(k*(a+1));
    sigma4 = cosh(k*(a-2));
    sigma5 = cosh(k*(a-1));
    sigma6 = cosh(k*(a+1));
    sigma7 = sinh(k*(2*a-1));
    sigma8 = sinh(2*a*k);
    sigma9 = cosh(2*a*k);
    sigma10 = 2*k*(a-1);
    sinhk = sinh(k);
    sinhak = sinh(a*k);
    sinh2k = sinh(2*k);
    coshk = cosh(k);
    coshak = cosh(a*k);
    cosh2k = cosh(2*k);
    line1 = 8*k-8*sinh2k-24*sinhak+16*sinhk+24*sigma1+8*sigma3-8*sigma2+12*k^2*sinhk + 12*k^2*sigma1-4*k^2*sigma3-4*k^2*sigma2+k^2*sinh(sigma10)-8*k*coshak-3*k^2*sinh2k-4*k^2*sigma7+8*a*k^3-16*k*cosh(k)+16*k*sigma5-8*k*sigma4;
    line2 = 8*k^3*coshak-2*k^3*sigma9-4*k^2*sinhak+3*k^2*sigma8-6*k^3-8*a^2*k^3+8*k*cosh2k-24*a*k*sigma5-8*a*k*sigma6+8*a*k*sigma4-12*a*k^3*coshak+2*a*k^3*sigma9-4*a^2*k^2*sinh2k+16*a*k^2*sinhak-4*a*k^2*sigma8-8*a*k^3*coshk;
    line3 = -8*a*k^2*sinhk+4*a*k^3*sigma5+4*a*k^3*sigma6+4*a*k^3*sigma4-2*a*k^3*cosh(sigma10)+4*a^2*k^3*coshak-32*a*k^2*sigma1+16*a*k^2*sigma2-4*a*k^2*sinh(sigma10)-24*a^2*k^2*sinhak+8*a^2*k^3*coshk+24*a*k*coshak+4*a*k^2*sinh2k;
    line4 = 8*a^2*k^2*sinhk + 4*a^2*k^3*(sigma5-sigma6-sigma4)+8*a*k^2*sigma7+24*a^2*k^2*sigma1+8*a^2*k^2*(sigma3-sigma2);
    if k > 0.05
        denom = 2*k^3*(coshk-1)*(4*coshk-4*k*sinhk+k^2+k^2*coshk-4);
    else
        denom = k^6/72 + k^8/1440 + k^10/67200; % Series expansion up to O(k^11)
    end
    part1 = lMagPrime^2 * (line1+line2+line3+line4)/denom;
    if k > 0.05
        line5 = 2*k-2*sinh2k-6*sinhak+4*sinhk+6*sigma1+2*sigma3-2*sigma2+2*k^2*sinhk+k^2*(sigma1-sigma3)-3*k*coshak-4*k*coshk+3*k*sigma5+k*sigma6-k*sigma4+2*k*coshk^2-6*a*k*sigma5-2*a*k*sigma6+2*a*k*sigma4-a*k^2*sinhak-a*k^2*sigma1+a*k^2*sigma3+a*k^2*sigma2+6*a*k*coshak;
    else
        line5 = a^2*(a-1)^2*(k^9/72 + k^11*(4*a^2-4*a+7)/4320); % Power series to O(k^13)
    end
    denom2 = k^3*(coshk-1)*(k*sinhk-2*coshk+2);
    part2 = -lMagPrime*lDcPrime*line5/denom2;
    TacEmag = part1+part2;
end
