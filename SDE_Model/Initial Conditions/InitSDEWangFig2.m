
T = 50; 
N = 10^4; 
h = T/N;
nParticles=100;
x = zeros(nParticles, N+1);

IC = 'H'; %Initial condition when every particle starts at 0.5.
x(:,1) = 0.5;


Rmin = 0.1; Rmax = 0.1;
bmin = 0.07; bmax = 0.07;
c=1;