%% Shifted start distribution test
%We use this script to test whether shifting the initial start distribution
%by a constant (with everything else left the same) causes the dynamics to 
%evolve in the same way apart from also being shifted by the same constant.
clear all
close all
set(0, 'DefaultAxesFontSize', 30)

%% Parameters
nParticles=50; %Number of individuals
T = 10^3; %End time of the simulation
N = 10^5; %Number of time steps
h = T/N;  %Step size

W=randn(nParticles, N+1); %Realisations of Brownian increments
x = zeros(nParticles, N+1); %Vector to save the opinions in
x(:,1) = rand(nParticles, 1); %Sampling the initial distribution from the
                              %uniform distribution
Ones = ones(1, nParticles); %For use in the function F_SDE

R = 0.2; %Half the width of the confidence intervals
%b=0; %Noise strength of 0 means the model is deterministic i.e the
      %ODE model
b = 0.01; %Noise strength


x0 = 0.2; %The size of the shift
xs = zeros(nParticles, N+1); %Vector to save the shifted dynamics in
xs(:,1) = mod(x(:,1)+x0,1); %Shifting the initial distribution by 
                            %x_0 and applies the boundary condition

%% Run both simulations at the same time

for n = 1:N %for each time step
    
   %Calculate forcing term for the original distribution dynamics
   [F, ~] = F_SDE(nParticles, x, R, n, Ones);
   %Calculate forcing term for the shifted distribution dynamics
   [FS, ~] = F_SDE(nParticles, xs, R, n, Ones);
    
   %Update the original dynamics using Euler-Maruyama
   x(:,n+1) = x(:,n) -h*F + b*sqrt(h).*W(:,n);
   x(:,n+1) = mod(x(:,n+1),1);  %Apply boundary condition
   
   %Update the shifted dynamics using Euler-Maruyama
   xs(:,n+1) = xs(:,n) -h*FS + b*sqrt(h).*W(:,n);
   xs(:,n+1) = mod(xs(:,n+1),1);  %Apply the boundary condition
end

%% Plots - use this section to plot

t = 0:h:T; %Creating the plot times in a vector

%Plot the opinion paths
figure() 
subplot(2,1,1) %Plot paths of the original dynamics
    for i = 1:nParticles
        plot(t,x(i,:),'.');
        hold on
        xlabel('Time'); ylabel('Opinion');
    end
    
subplot(2,1,2) %Plot paths of the shifted dynamics
    for i = 1:nParticles
        plot(t,xs(i,:),'.');
        hold on
        xlabel('Time'); ylabel('Opinion');
    end

%Plot the histograms of the distributions of opinions
figure() 
subplot(2,1,1) %Plot for the original dynamics
histogram(x(:,N+1),50) %Histogram using 50 bins
axis([0,1,0,max(histcounts(x(:,N+1),50))+1])
xlabel('Opinion'); 
ylabel('Number of Individuals');

subplot(2,1,2) %Plot for the shifted dynamics
histogram(xs(:,N+1),50) %Histogram using 50 bins
axis([0,1,0,max(histcounts(x(:,N+1),50))+1])
xlabel('Opinion'); 
ylabel('Number of Individuals');

%% Functions

 function [F,nRij] = F_SDE(nParticles,x,R,n, Ones)
 %F_SDE Calculates F for SDE solver
 
    xn = x(:,n)'; %Extracts all individuals opinions at time n in a row
    yn = xn';     %Column vector of all the opinions at time n
    X = xn(Ones,:); %Makes a square matrix where each column m is a load of
                    %copies of the opinion of individual m.
    Y = yn(:,Ones); %Square matrix, each row m is copies of the opinion of
                    %individual m.
    Rij = X-Y;      %Calculates the difference in the opinions of all the
                    %individuals with the sign


    mask = ( abs(Rij) > 1/2 ); % deal with periodicity
    Rij(mask) = Rij(mask) - sign(Rij(mask))*1;
    
    nRij = abs(Rij)<=R;
    nRij = sum(nRij,1); %Finds the number of individuals who are
                        %within R of each individual in a row vector
    Rij(abs(Rij)>R) = 0; %Works out which opinions (i,j) aren't within R
                         %and sets all these Rij to be zero
    sumRij = sum(Rij,1); %Each column m is the sum of opinions in row m,
                         %which is the opinions that influence individual m
    F = (nParticles.^(-1) .* sumRij)';
 end