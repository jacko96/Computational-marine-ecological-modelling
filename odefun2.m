function dPdt = odefun2(t,P,param)
% Setup parameters
n = param.n; % no. cells
%D = param.D; % depth
dz = param.dz;
%z = param.z;
u = param.u;
d = param.d;
l = param.l;
H = param.H;
pmax = param.pmax;
Iin = param.Iin;
k = param.k;
Kbg = param.Kbg;

Ja = zeros(n+1,1);
Jd = zeros(n+1,1);
% Advective and diffusion flux:
for i = 2:n
    Ja(i) = u*P(i-1);
    Jd(i) = -d*(P(i)-P(i-1))/dz;
end
J = Ja+Jd; % Total flux

% cumsum for int
I = Iin*exp(-cumsum( (P*k + Kbg)*dz));

for i = 1:n
    p(i) = pmax*I(i)/(H+I(i));
    g(i) = p(i) - l;
    dPdt(i) = g(i)*P(i) - (J(i+1)-J(i))/dz;
    
end
p = p';
g = g';
dPdt = dPdt';
end