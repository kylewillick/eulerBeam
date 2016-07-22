function [ z, u ] = ModeShape_pos2( mat, Tx, Tz ,L, diameter,pos)
%MODESHAPE_POS Plots the resonant mode shape given from NumericFreq_pos

% Toss out a row of the solution matrix (as det(mat)=0) and add an
% arbitrary value
tmp = mat(1:8,:);
tmp = [tmp;[1 0 0 0 0 0 0 0 0]];
% Solve system of equations
coeffs = linsolve(tmp,[0;0;0;0;0;0;0;0;1]);
% Determine d2u/dz2 in Prime coordinates
newu = Tx/(diameter/2);
deltaz = 1/max(size(Tz));
du = diff(newu)/deltaz;
d2u = diff(du)/deltaz;
% Determine the modeshape (Note the k+ and k- are pulled from mat)
zLength = max(size(Tz))-2;
z = linspace(0,1,zLength);
posIndex = floor(zLength*pos);
z1 = z(1:posIndex);
z2 = z(posIndex+1:end);
u = zeros(1,zLength);
u(1:posIndex) = coeffs(1)*cos(mat(2,2)*z1)+coeffs(2)*sin(mat(2,2)*z1)+coeffs(3)*cosh(mat(2,4)*z1)+coeffs(4)*sinh(mat(2,4)*z1)+coeffs(9)*d2u(1:posIndex);
u(posIndex+1:end) = coeffs(5)*cos(mat(2,2)*(1-z2))+coeffs(6)*sin(mat(2,2)*(1-z2))+coeffs(7)*cosh(mat(2,4)*(1-z2))+coeffs(8)*sinh(mat(2,4)*(1-z2))+coeffs(9)*d2u(posIndex+1:end);

end

