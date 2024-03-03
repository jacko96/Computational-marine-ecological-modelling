% Exercise 3
clear;close all;clc;
% 1. Define parameters
param.n = 100; % no. cells
param.D = 100; % depth
param.dz = param.D/param.n; % grid size
param.z =  param.dz/2:param.dz:param.D; % grid
param.u = 0.04; % current (gravitational flux) / sinking velocity (m/h)
param.d = 1; % Diffusion constant (cm2/s)
param.Iin = 100; % Light surface input (micromol photons/m2 s)
param.k = 0.01*6*10^(-10); % light absorbtion coefficient of plankton (m2/cell)
param.Kbg = 0.045; % background turbidity (1/m)
param.l = 0.01; % plankton loss rate (1/h)
param.HI = 20; % Half saturation constant of light limited growth (micromol photons/m2 s)
param.HN = 0.0425; % Half saturation constant of nutrient limited growth (mmol nutrient/m3)
param.pmax = 0.04; % maximum growth rate (both plankton and nutrients) (1/h)
param.y = 1*10^(-9); % Nutrient content of phytoplankton (mmol nutrient/cell)
param.m = 0.01; % specific loss rate (plankton mortality) (1/h)
param.ND = 5; % Nutrient conc. at bottom (mmol nutrient/m3)

% 5. solve
P0 = ones(param.n,1)*10; % Initial conditions:
N0 = ones(param.n,1)*5;
Y0 = [P0;N0];
tspan = 0:100:10000;
[t,Y] = ode45(@(t,Y)odefun3(t,Y,param),tspan,Y0);
[~, pI,pN] = odefun3(t(end),Y(end,:),param); % Finding the limiting factors (light or minerals)

Y = Y';
P = Y(1:param.n,:);
N = Y(param.n+1:end,:);

% 6. plot
% Plankton
figure(1)
image(P,'CDataMapping','scaled')
colorbar
xlabel('Time [100 hours]')
ylabel('Depth [meters]')
title('Phytoplankton conc. [cells/m3] over time')

% Nutrients
figure(2)
image(N,'CDataMapping','scaled')
colorbar
xlabel('Time [100 hours]')
ylabel('Depth [meters]')
title('Nutrient conc. [mmol nutrient/m3]')

% Limiting factors
figure(3)
plot(pI,-param.z)
hold on
plot(pN,-param.z)
hold off
legend('Light','Minerals')
title('Limiting factors')
ylabel('Depth [meters]')
