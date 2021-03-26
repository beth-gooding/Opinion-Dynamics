function PDE_SDE_Compare()
%This script compares the solutions of the SDE and PDE models for
%opinion dynamics by running the same set up in both models. We
%compare the output by plotting the solutions rho of the PDE model
%against the histogram of the final distribution of opinions in the
%SDE model.

close all
clear all
%% Common parameters

%Noise Parameters
sigma1 = 0.20/sqrt(2*pi); %For a single simulation
%For the parameter experiment
%sigma1 = linspace(0.05, 0.25, 15)/sqrt(2*pi); 

%Width of confidence intervals
Rfull = 0.2; %For a single simulation
%For the parameter experiment
%Rfull = linspace(0.01,0.37, 15);

%Normal Distribution Parameters for Initial Condition
mu = 0.5; %Mean
Sigma = 1/sqrt(40); %Standard Deviation

%End time of simulations 
T = 2*pi;
% T = 20*2*pi; %For the sigma = 0.25 case

%Order Parameters
QC = zeros(length(Rfull),length(sigma1));
Q = QC;

%% PDE Solver
N_p = 300; %Number of grid points
h_p = 1/N_p; %Stepsize of equispaced points
x_p = (0:N_p-1)'*h_p; %Equispaced periodic points on [0,1]

% Construct spectral differentiation matrix on [0,2pi]:
h_D = 2*pi/N_p; %Equispaced points on [0,2pi]
column = [0 .5*(-1).^(1:N_p-1).*cot((1:N_p-1)*h_D/2)];
D = toeplitz(column,column([1 N_p:-1:2])); %First derivative

% D on [0,1]
D = 2*pi*D; %Transform D so that it is differentiating on [0,1]

D2 = D^2; %Second derivative, found by squaring D.

Int = h_p*ones(1,N_p); %integration weights on [0,1], used to 
                       %normalise initial conditions and calculate
                       %the order parameter in the PDE model

tic
for i = 1:length(Rfull) %For each CI width
     R = Rfull(i); %For a specific width
     
    for j = 1:length(sigma1) %For each noise strength
      sigma = sigma1(j); %For a specific noise strength
        
      %convolution matrix on [0,1]
      C = ConvolutionMatrix(x_p,R); %Calculates the matrix that will do 
                                    %integration and interpolation
      C_Q = ConvolutionMatrixQ(x_p,R); %Needed to calculate Q, the order   
                                       %parameter
        
        
      plotTimes=[0:T/100:T]; %Times to plot our solution on

      %Creates the initial condition for our problem
      rhoIC = initialCondition(x_p); 

      %Solves the resulting ODE now that we are able to
      %approximate the right hand side of the PDE which
      %deals with all the derivatives with respect to x
      [plotTimes,rho] = ode15s(@rhs, plotTimes, rhoIC); 
                                             
                                                          
       QC(i,j) = Qc(rho); %Calculates the order parameter
 
%% SDE Solver

N = 2000; %Number of timesteps
h = T/N; %Size of timesteps
nParticles=300; %Number of individuals
nSim = 100; %Number of simulations
Ones = ones(1, nParticles); %For F_SDE
W=randn(nParticles, N+1); %Realisations of Brownian Motion

Results = zeros(nParticles*nSim, N+1); %For saving opinions
Q_sim = zeros(nSim,1);
for p = 1:nSim %for each simulation
    x = zeros(nParticles, N+1);
    %Sample an initial condition for a standard normal
    %distribution and scale it to be from the right one
    x_1 = randn(nParticles,1);
    x_1 = Sigma*x_1 + mu;
    x(:,1) = mod(x_1,1); %Apply the boundary condition
   
    p/nSim; %Shows progress if we want an indicator
    
    for n = 1:N %for each timestep
        %Run the SDE dynamics
        [F, nRij] = F_SDE(nParticles, x, R, n, Ones);
        x(:,n+1) = x(:,n) -h*F + sigma*sqrt(h).*W(:,n);
        x(:,n+1) = mod(x(:,n+1),1); %Apply the boundary condition  
                    
    end
    %save the value of the order parameter for each simulation
    Q_sim(p) = nParticles^(-2)*sum(nRij);  
            
    %Save the results of the simulation into some rows of the 
    %matrix so that each block of nParticles rows is the 
    %results of a simulation.
    Results((1:nParticles)+nParticles*(p-1),:) = x;
end
    Q(i,j) = mean(Q_sim);
    end
end
time = toc; %finds the run time of all the simulations

%% Plot Section

% Plot the results of the parameter experiment
% figure()
% pcolor(Rfull,sigma1, Q'); shading interp;

%Plot density functions and histograms

        timesteps = h*(0:N); %calculate plot times
        xhist = 0:0.02:1; %Calculate bin edges
        %Sets the edges of the histogram boxes to be the x values
        edges = xhist;
        figure() %Plots the solution to the problem as time evolves
        plots = [1,21,41,61,81,101];
        colour = ['k','b','r','g','m','c'];
        s = 0;
        %Create a vector which only saves the results we use to plot
        Results_save = zeros(nParticles*nSim,length(plots));
        for iPlot = plots
            %Find the closest plot times of the results from
            %the SDE model to the PDE plot times
            [~,I] = min(abs(plotTimes(iPlot) - timesteps));
            s = s+1;
            %hold off
            plot(x_p, rho(iPlot,:)',colour(s)) %plot PDE solution
            hold on
            
            %counts the number of particles in each bin, 
            %and resaves the edges
            [bincounts,edges] = histcounts(Results(:,I),edges);
          
            %Finds the width of each bin
            dx = abs(edges(2)-edges(1));
             
            %Normalises the bincounts so that the total area of the
            %histogram is 1.
            bincounts = bincounts/(nParticles*nSim*dx); 

            Area = sum(bincounts*dx); %Check that the area of the 
                                       %histogram is 1.

            %Plots the histogram bin counts
            scatter(edges(1:end-1)+dx/2,bincounts, colour(s))  
            ylim([0,14]);
            
            %Saves the results we used to plot the SDE histogram
            Results_save(:,s) = Results(:,I);
        end
        
%Choose a sensible filename and save the data
clear Results %clear all results to keep file size small(ish)
save('PDE_SDE_Compare_4')

%% PDE Functions
 function [drhodt] = rhs(t,rho)
      %This is the PDE Model
        drhodt = (sigma^2/2)*D2*rho + D*(rho.*(C*rho));
    end

    function Q = Qc(rho)
         %Calculate the order parameter
       Q = Int*(rho(end,:)'.*(C_Q*(rho(end,:)'))); 
    end

    function rhoIC = initialCondition(y) 
        %Sets up the initial condition
        rhoIC = exp(-0.5*((y-mu)/Sigma).^2);
        Z = Int*rhoIC;
        rhoIC = rhoIC/Z;
    
    end

end