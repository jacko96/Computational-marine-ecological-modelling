function [dYdt,pI,pN] = odefun3(t,Y,param)
% Setup parameters
n = param.n; % no. cells
dz = param.dz;
%z = param.z;
u = param.u;
d = param.d;
l = param.l;
HI = param.HI;
HN = param.HN;
pmax = param.pmax;
Iin = param.Iin;
k = param.k;
Kbg = param.Kbg;
y = param.y;
m = param.m;
ND = param.ND;

% Splitting the Y variable into P (plankton) and N (nutrients):
P = Y(1:n);
N = Y((n+1):end);

% Advection+diffusion flux for plankton:
Ja = zeros(n+1,1);
Jd = zeros(n+1,1);
for i = 2:n
    Ja(i) = u*P(i-1);
    Jd(i) = -d*(P(i)-P(i-1))/dz;
    
end
JP = Ja+Jd; % Total plankton flux
% Neumann boundary conditons:
JP(1) = 0;
JP(end) = 0;

% Advection+diffusion flux for nutrients:
JaN = zeros(n+1,1);
JdN = zeros(n+1,1);
for i = 2:n
    JaN(i) = 0;
    JdN(i) = -d*(N(i)-N(i-1))/dz;
    
end
JN = JaN+JdN; % Total nutrient flux
% Neumann boundary conditions:
JN(1) = 0;
JN(end) = -d*(ND-N(n-1))/dz;

% Light intensity as a function of depth:
I = Iin*exp(-cumsum( (P*k + Kbg)*dz));

for i = 1:n
    
    % Production due to light or minerals (limiting factor):
    pI(i) = I(i)/(HI+I(i));
    pN(i) = N(i)/(HN+N(i));
    p(i) = pmax*min( pI(i) , pN(i) );
    g(i) = p(i) - l;
    
    % PDEs
    dPdt(i) = g(i)*P(i) - m*P(i) - (JP(i+1)-JP(i))/dz; % Plankton
    dNdt(i) = - y*g(i)*P(i) + y*m*P(i) - (JN(i+1)-JN(i))/dz; % Nutrients
    
end
pI = pI';
pN = pN';

dPdt = dPdt';
dNdt = dNdt';
dYdt = [dPdt;dNdt];

end