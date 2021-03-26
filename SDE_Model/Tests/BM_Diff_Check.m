
%% Brownian Motion/Diffusion Equation Check
%Testing to see whether the SDE model with no interaction is solved
%correctly by comparing it to the 1D diffusion equation, which we know has
%the same solution by the Fokker-Planck equation

clear all
close all
set(0, 'DefaultAxesFontSize', 24)

%% Parameters
nParticles=10000; %Number of individuals
T = 1; %End time of simulation
N = 10^4; %Number of time steps
h = T/N; %Stepsize
Ones = ones(1, nParticles); %For the function F_SDE
W=randn(nParticles, N+1); %Realisations of Brownian Increments

x = zeros(nParticles, N+1); %Vector to save opinion paths in
x(:,1) = randn(nParticles, 1); %Initial condition is a sample from the
                               %standard normal distribution

R = 0; %Turning off interaction term
b = 1; %Strength of noise

for n = 1:N %For each time step
   %Update the opinions of the individuals using Euler $-$
   %we don't need F_SDE as the interaction term is turned off.
   x(:,n+1) = x(:,n) + b*sqrt(h)*W(:,n); 
end
                          
% t = 0:h:T;
% figure()
%     for i = 1:nParticles
%         plot(t,x(i,:),'.');
%         hold on
%     end

%% Plotting
%Evaluate the exact solution
y0 = 5;
y = -y0:y0/45:y0; %Setting y values to evaluate exact solution at
p = (1/sqrt(4*pi*T)).*exp(-(y.^2)./(4*T)); %Calculates the exact solution 
                                           %for sigma =1
                                           
%create the histogram of the distribution of opinions in the SDE solution                                           
edges = y; %Sets the edges of the histogram boxes to be the y values
[bincounts,edges] = histcounts(x(:,N+1),edges); 
%counts the number of particles in each bin, and resaves the edges
dx = abs(edges(2)-edges(1)); %Finds the width of each bin
bincounts = bincounts/(nParticles*dx); %Normalises the bincounts so that 
                                     %the total area of the histogram is 1.

Area = sum(bincounts*dx) %Check that the area of the histogram is 1.

figure()
%Plots the normalised bincounts at the y coordinate of the middle of each
%bin except for the count in the very last bin
plot(edges(1:end-1)+dx/2,bincounts, 'o', 'markersize', 12) 
hold on
plot(y,p,'.k', 'markersize', 24)
legend('Observed','Exact')
xlabel('$x$', 'Interpreter', 'latex'); 
ylabel('$\rho(x,1)$', 'Interpreter', 'latex');
