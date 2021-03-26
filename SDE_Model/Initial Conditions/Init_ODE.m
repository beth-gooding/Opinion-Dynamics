
%% Init_ODE: Parameters for the ODE model for use in SDE_SolverSave.m

T = 500; %End time of simulation
N = 5*10^4; %Number of time steps
h = T/N; %Size of time step
nParticles=625; %Number of individuals

%Vector to save the opinion paths in
x = zeros(nParticles, N+1);

IC = 'U'; %Uniformly distributed initial condition indicator

%Use the following line to sample the same initial condition as for
%the HK model simulations by setting the seed of the random number
%generator
ra = 1; rng(ra);

x(:,1) = rand(nParticles, 1); %Sample the initial condition from
                              %the uniform distribution

Rmin = 0.01; Rmax = 0.01; %Set the minimum and maximum size for 
                          %half the width of the confidence interval
bmin = 0; bmax = 0; %Set the minimum and maximum noise strength.
                    %Here it is zero as the ODE model has no noise.
c=1; %The number of linearily spaced values we would like in the 
     %intervals [Rmin, Rmax] and [bmin, bmax]