
%% InitSDE: Parameters for the SDE model for use in SDE_SolverSave.m
%This set up is used for general simulations of the SDE model, doing
%both individual simulations and parameter experiments depending on
%the value of c.

T = 10^4; %End time of simulation
N = 5*10^4; %Number of time steps
h = T/N; %Size of the time steps
nParticles=300; %Number of individuals
x = zeros(nParticles, N+1); %Vector for saving opinion paths in 

%Sample the initial condition from a uniform distribution
x(:,1) = mod(rand(nParticles, 1),1);
IC = 'U'; %Uniformly distributed initial condition indicator

%Set the initial condition so that all opinions are 0.5.
% x(:,1) = 0.5; 
% IC = 'H'; %Initial condition indicator

%% Parameters of the confidence interval widths and noise strength
%Comment out the block which you would not like to use.

%Set up for simulation using one parameter pair

Rmin = 0.2; Rmax = 0.2; %Set min and max confidence interval widths
bmin = 0.09; bmax = 0.09; %Set min and max noise strengths
c=1; %Sets the number of linearily spaced values in the intervals
     %[Rmin, Rmax] and [bmin, bmax], here we just have one pair.

%Set up for parameter experiments using multiple pairs

Rmin = 0.01; %minimum confidence interval width
Rmax = 0.35; %maximum confidence interval width
bmin = 0;  %minimum noise strength
bmax = 0.17; %maximum noise strength
c=30; %Sets the number of linearily spaced values in the intervals
      %[Rmin, Rmax] and [bmin, bmax]. Our parameter experiment 
      %will have c^2 pairs.
