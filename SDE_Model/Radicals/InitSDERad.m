
%% InitSDERad.m : Parameters for the SDE Model including radicals
%For use in the script Rad_SDE_SolverSave.m

T = 10^4; %Simulation end time
N = 10^5; %Number of time steps
h = T/N;  %Size of time step

nParticles=50; %Number of normals

x = zeros(nParticles, N+1); %Vector to save opinion paths in

%Sampling the initial condition from the uniform distribution
x(:,1) = mod(rand(nParticles, 1),1);
IC = 'U'; %Uniformly distributed initial condition indicator

Rmin = 0.05; Rmax = 0.05; %Set min and max confidence interval widths
bmin = 0.01; bmax = 0.01; %Set min and max noise strengths
c=1; %Sets the number of linearly spaced values in the intervals
     %[Rmin, Rmax] and [bmin,bmax], here we have just one pair.


