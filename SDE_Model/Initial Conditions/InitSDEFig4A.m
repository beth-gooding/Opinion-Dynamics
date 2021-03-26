
%% InitSDEFig4A: Parameters for the SDE model for use in SDE_SolverSave.m
%This set up is used for simulations of the SDE model used to replicate
%Figure 4A of the paper by Wang et. al.

T = 500; %Simulation end time
N = 5*10^3; %Number of time steps
h = T/N; %Step size
nParticles=300; %Number of individuals
x = zeros(nParticles, N+1); %Vector to save opinion paths in

%Sample the initial condition from a uniform distribution
x(:,1) = mod(rand(nParticles, 1),1);
IC = 'U'; %Uniformly distributed initial condition indicator

%Set the initial condition so that all opinions are 0.5.
% x(:,1) = 0.5; 
% IC = 'H'; %Initial condition indicator

Rmin = 0.05; Rmax = 0.05; %Set min and max confidence interval widths
bmin = 0.01; bmax = 0.01; %Set min and max noise strengths
c=1; %Sets the number of linearily spaced values in the intervals
     %[Rmin, Rmax] and [bmin, bmax], here we just have one pair.