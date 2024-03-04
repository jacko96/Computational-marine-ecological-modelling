% Convergence plot
clear;close all;clc;

% Define parameters
run('Parameters.m')
%param.n = 100; % no. cells
n_max = 200; % no. cells
dz_max = param.D/n_max; % grid size
z_max =  dz_max/2:dz_max:param.D; % grid
max_length = length(z_max);

dz_vals = [];
Ysols = [];
Psols = [];
Nsols = [];
Dsols = [];

n_vals = 10:10:n_max;

%for i = 1:length(n_vals)
for i = 1:length(n_vals)
% Defining grid
param.n = n_vals(i); % no. cells
param.dz = param.D/param.n; % grid size
param.z =  param.dz/2:param.dz:param.D; % grid

dz_vals(i) = param.dz;

% 5. solve
P0 = ones(param.n,1)*10; % Initial conditions:
N0 = ones(param.n,1)*5;
D0 = ones(param.n,1);
Y0 = [P0;N0;D0];
tspan = 0:1:2500;

[t,Y] = ode45(@(t,Y)odefun4(t,Y,param),tspan,Y0);
[~, pI,pN] = odefun4(t(end),Y(end,:),param); % Finding the limiting factors (light or minerals)

Y = Y';
P = Y(1:param.n,:);
N = Y(param.n+1:(2*param.n),:);
D = Y((2*param.n+1):end,:);

Pint = interp1(param.z,P(:,end),z_max)';
Nint = interp1(param.z,N(:,end),z_max)';
Dint = interp1(param.z,D(:,end),z_max)';
Yint = [Pint;Nint;Dint];

Psols(:,i) = Pint;
Nsols(:,i) = Nint;
Dsols(:,i) = Dint;
Ysols(:,i) = Yint;
end

%% Error
% Plankton error
PE = [];
for i = 1:length(n_vals)-1
    Perror = abs(Psols(:,i) - Psols(:,end));
    Perror = Perror(~isnan(Perror));
    PE(i) = norm(Perror,1); 
    
    Nerror = abs(Nsols(:,i) - Nsols(:,end));
    Nerror = Nerror(~isnan(Nerror));
    NE(i) = norm(Nerror,1);
   
    Derror = abs(Dsols(:,i) - Dsols(:,end));
    Derror = Derror(~isnan(Derror));
    DE(i) = norm(Derror,'inf'); 
end

figure()
subplot(3,1,1)
plot(dz_vals(1:end-1),PE)
%xlabel('\Delta z')
title('Error with grid size')
ylabel('Plankton error')
subplot(3,1,2)
plot(dz_vals(1:end-1),NE)
%xlabel('\Delta z')
ylabel('Nutrient error')
subplot(3,1,3)
plot(dz_vals(1:end-1),DE)
xlabel('\Delta z')
ylabel('Detritus error')
