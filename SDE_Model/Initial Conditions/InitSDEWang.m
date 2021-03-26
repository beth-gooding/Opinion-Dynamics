
%% InitSDEWang: Parameters for the SDE model for use in SDE_SolverSave.m
%This set up is used for simulations of the SDE model used to replicate
%Figure 3 of the paper by Wang et. al.

T = 10^4; %End time of simulation
N = 5*10^4; %Number of time steps
h = T/N; %Size of the time steps
nParticles=300; %Number of individuals
x = zeros(nParticles, N+1); %Vector for saving opinion paths in 

%Set the initial condition so that all opinions are 0.5.
x(:,1) = 0.5; 
IC = 'H'; %Initial condition indicator

%Initial condition for parameter experiments using multiple pairs

Rmin = 0.01; %minimum half confidence interval width
Rmax = 0.35; %maximum half confidence interval width
bmin = 0;  %minimum noise strength
bmax = 0.17; %maximum noise strength
c=30; %Sets the number of linearily spaced values in the intervals
      %[Rmin, Rmax] and [bmin, bmax]. Our parameter experiment 
      %will have c^2 pairs.



