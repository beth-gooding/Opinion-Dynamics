%% InitSDERadExperiment.m: Parameters for the SDE model with radicals.
%For use in the script Rad_Solver_SDESave_Proportion.m

T = 10^4; %End time of simulation
N = 10^5; %Number of time steps
h = T/N;  %Size of the time steps

nParticles=50; %Number of normals
x = zeros(nParticles, N+1); %Vector to save the opinion paths in

%Setting the initial condition to be 50 points that are equally
%spaced across the interval [0,1]
i = linspace(1,nParticles,nParticles);
x(:,1) = i/(nParticles+1);
IC = 'Uni'; %Uniformly spaced initial condition indicator

Rmin = 0.01; Rmax = 0.5; %Set the min and max confidence interval widths
bmin = 0.01; bmax = 0.01; %Set the min and max noise strengths
