
%% HK.m: This script implements the Hegselmann-Krause BC Model
clear all
close all
set(0, 'DefaultAxesFontSize', 30)
set(0, 'DefaultLineLineWidth', 1.5)
%% Initial Setup - Simulating the Hegselmann-Krause Model

%Half the width of the Confidence Interval:
% R = 0.01;
%R = 0.15;
R = 0.25;
%R = [0.01,0.15,0.25];


%HK Investigations
nParticles=625; %Number of individuals in simulations
N = 15; %Number of timesteps
T = N; %Final time
h = T/N; %Size of timestep


x = zeros(nParticles, N+1); %For storing opinion values
rng(1) %Makes sure the start distribution is always the same
x(:,1) = rand(nParticles, 1); %Sampling the initial condition 
Ones = ones(1, nParticles); %For the function F_HK

%% Running the model

%If the data exists already, load data and plot the results
filename = sprintf('test_dis_data_%d_%d_%d_%d.mat', T,N,nParticles,R*100);
if exist(filename, 'file')==2
    load(filename)
    
      figure() 
      subplot(2,1,1)
    for i = 1:nParticles %Plot opinion paths for every individual
            plot(t,x(i,:),'.');
            hold on
            xlabel('Time')
            ylabel('Opinion'); 
    end

    subplot(2,1,2)
    histogram(x(:,N+1),100) %plot the distribution of opinions at the end
                            %time of the simulation
    xlim([0,1]); xlabel('Opinion');
    ylabel('Count');
    
    disp('loaded data')
    
else %If the data does not exist, run the simulation

    for k = 1:length(R)
        
        for n = 1:N %Updating the opinions for every time step

        [F, nRij] = F_HK(x,R(k),n, Ones); %Calculates the updated opinions
                                          %for each individual

        x(:,n+1) = F; %Saves the updated opinions

        end
        t = 0:h:T; %calculates the times we calculated the opinions at
        figure()
        for i = 1:nParticles %plots the opinion paths
            plot(t,x(i,:), '.');
            hold on
            xlabel('Time');
            ylabel('Opinion'); 
        end

    figure() 
    histogram(x(:,N+1),100) %plots the distribution of opinions at the end
                            %time step
    xlim([0,1]); xlabel('Opinion');
    ylabel('Count');

    end
    %save the new data so that it can be reused in the future
    newfilename = sprintf('test_dis_data_%d_%d_%d_%d', T,N,nParticles,R*100);
    save(newfilename)
end

%% Functions
 function [F, nRij] = F_HK(x,R,n, Ones)
 
    xn = x(:,n)'; %Extracts all individuals opinions at time n in a row
    yn = xn';     %Column vector of all the opinions at time n
    X = xn(Ones,:); %Makes a square matrix where column m contains N
                    %copies of the opinion of individual m.
    Y = yn(:,Ones); %Square matrix, where row m contains N copies of 
                    %the opinion of individual m.
    Rij = X-Y;      %Calculates the distances between all pairs of opinions 
                    %with entry (i,j) equal to xj-xi.

    mask = (abs(Rij)>R); %Identifies which pairs of opinions are more than 
                         %R apart
    nRij = abs(Rij)<=R;  %Identifies which pairs of opinions are considered
                         %close
    nRij = sum(nRij,1);  %Creates a row vector where entry i tells us how
                         %many individuals are close to individual i.
    X(mask) = 0;         %Sets all entries corresponding to distances that
                         %are bigger than R to be zero.
    sumRij = sum(X,2)';  %Creates a row vector where entry i is the sum of
                         %all opinions that are within R of opinion xi.
    
    F = (nRij.^(-1) .* sumRij)'; %Calculates the mean of all opinions that
                                 %are close to each individual in a vector.
                                 %This gives the updated opinions for the
                                 %next time step.
 
 end


