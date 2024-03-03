% Exercise 1
clear;close all;clc;
% 1. Define parameters
param.n = 100; % no. cells
param.D = 100; % depth
param.dz = param.D/param.n;

% 2. Set up grid
param.z =  param.dz/2:param.dz:param.D;

% 3. Set in. cond.
phi0 = normpdf(param.z,param.D/2,param.D/10)';

% 4. Derivative function (time, phi)
param.u = 1; % current
param.d = 10; % Diffusion constant
%Ja = zeros(n+1,1);
%Jd = zeros(n+1,1);

% Set up odefun

% 5. solve
tspan = 0:100;
[t,phi] = ode45(@(t,phi)odefun(t,phi,param),tspan,phi0);
phi = phi';

% 6. plot
image(phi,'CDataMapping','scaled')
colorbar
xlabel('Time')
ylabel('Depth')
title('Conc. over time')