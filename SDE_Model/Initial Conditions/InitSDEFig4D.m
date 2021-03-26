T = 400; 
%T = 10^4; 
N = 4*10^3; 
h = T/N;
nParticles=300;
x = zeros(nParticles, N+1);
IC = 'U'; %Uniformly distributed initial condition
x(:,1) = mod(rand(nParticles, 1),1);

%IC = 'N'; %Normally distributed initial condition
%x(:,1) = mod(randn(nParticles, 1),1);

% IC = 'H'; %Initial condition when every particle starts at 0.5.
% x(:,1) = 0.5;


Rmin = 0.2; Rmax = 0.2;
bmin = 0.07; bmax = 0.07;
c=1;