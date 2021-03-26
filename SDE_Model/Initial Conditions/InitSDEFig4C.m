T = 500; 
%T = 10^4; 
N = 5*10^3; 
h = T/N;
nParticles=300;
x = zeros(nParticles, N+1);

IC = 'U'; %Uniformly distributed initial condition
x(:,1) = mod(rand(nParticles, 1),1);

%IC = 'N'; %Normally distributed initial condition
%x(:,1) = mod(randn(nParticles, 1),1);

% IC = 'H'; %Initial condition when every particle starts at 0.5.
% x(:,1) = 0.5;


Rmin = 0.05; Rmax = 0.05;
bmin = 0.03; bmax = 0.03;
c=1;