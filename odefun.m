function dphidt = odefun(t,phi,param)
% Setup parameters
n = param.n; % no. cells
%D = param.D; % depth
dz = param.dz;
%z = param.z;
u = param.u;
d = param.d;

Ja = zeros(n+1,1);
Jd = zeros(n+1,1);

% Advective and diffusion flux:
for i = 2:n
    Ja(i) = u*phi(i-1);
    Jd(i) = -d*(phi(i)-phi(i-1))/dz;
end
J = Ja+Jd; % Total flux

for i = 1:n
    dphidt(i) = - (J(i+1)-J(i))/dz;
end
dphidt = dphidt';
end