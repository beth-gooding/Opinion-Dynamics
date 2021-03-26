function PDE_SDE_Compare_Rad()
%This script compares the solutions of the SDE and PDE models for
%opinion dynamics by running the same set up in both models. We
%compare the output by plotting the solutions rho of the PDE model
%against the histogram of the final distribution of opinions in the
%SDE model. Here we include the influence of radicals.
clear all
close all
%% Common parameters

%Noise Parameters
sigma1 = 0.20/sqrt(2*pi);

%Width of confidence intervals
R = 0.2;

%Normal Distribution Parameters for initial condition
mu = 0.5; %Mean
Sigma = 1/sqrt(40); %Standard deviation

%Radical Distribution Parameters
mu_r = 0.8; %Mean
Sigma_r = 1/sqrt(40); %Standard deviation

%End time of simulations
T = 2*pi;

nParticles=300; %Number of normals in the SDE model
nrParticles = nParticles/10; %Number of radicals in the SDE model

%% PDE Solver
N_p = 300; %Number of grid points
h_p = 1/N_p; %Stepsize of periodic points
x_p = (0:N_p-1)'*h_p; %Equispaced periodic points on [0,1]


% Construct spectral differentiation matrix on [0,2*pi]:
h_D = 2*pi/N_p; %Equispaced on [0,2pi]
column = [0 .5*(-1).^(1:N_p-1).*cot((1:N_p-1)*h_D/2)];
D = toeplitz(column,column([1 N_p:-1:2])); %First derivative

% D on [0,1]
D = 2*pi*D; %Transform D so that it is differentiating on [0,1]

D2 = D^2; %Second derivative, found by squaring D.

Int = h_p*ones(1,N_p); % integation weights on [0,1] for 
                       % normalisation tasks

for i = 1:length(R) %for each CI width
     R = R(i);
    for j = 1:length(sigma1) %FOr each noise strength
      sigma = sigma1(j);
        
      %convolution matrix on [0,1]
      C = ConvolutionMatrix(x_p,R); %Calculates the matrix that will do 
                                    %integration and interpolation
%     C_Q = ConvolutionMatrixQ(x_p,R); %Needed to calculate Q, the order   
%                                      %parameter
        
        
      plotTimes=[0:T/100:T]; %Times to plot our solution on

      %Creates the initial condition for our problem
      rhoIC = initialCondition(x_p); 

      %Solves the resulting ODE now that we are able to
      %approximate the right hand side of the PDE which
      %deals with all the derivatives with respect to x
      [plotTimes,rho] = ode15s(@rhs, plotTimes, rhoIC); 
                                             
                                                          
%       QC(i,j) = Qc(rho); %Calculates the order parameter
        

    end
end

%% SDE Solver
N = 2000; %Number of timesteps
h = T/N; %Size of timesteps

nSim = 100; %Number of simulations
Ones = ones(1, nParticles); %For the function F_SDE
W=randn(nParticles, N+1); %Realisation of Brownian Increments
rOnes = ones(nParticles+nrParticles, 1); %For F_SDE
Results = zeros(nParticles*nSim, N+1); %To save opinions in

%For use when we do sigma1 = 0.25/sqrt(2pi) to be able to save
%all the data
% oops = [1,2001,4001,8001,16001,24001,40001];
% Results = zeros(nParticles*nSim, length(oops));

for i = 1:nSim %for each simulation
    %Sample the radical distribution for the SDE
    r = zeros(nrParticles, N+1);
    r(:,1) = randn(nrParticles, 1);
    r(:,1) = mod(Sigma_r*r(:,1) + mu_r,1); %Apply BC
    
    %Sample the initial distribution for the SDE Model
    x = zeros(nParticles, N+1);
    x_1 = randn(nParticles,1);
    x_1 = Sigma*x_1 + mu;
    x(:,1) = mod(x_1,1); %apply boundary condition
    
    xr = x; %For saving the opinions of normals
    z = [xr;r]; %For using in the function F_SDE
   
    i/nSim; %indicating how many simulations we have done
    
    for n = 1:N %for each time step
        %Run the dynamics of the SDE model
         [Fr, ~] = F_SDE(nParticles, z, R, n, rOnes);
            xr(:,n+1) = xr(:,n) -h*Fr(1:nParticles)+...
                        sigma1*sqrt(h).*W(:,n);
            xr(:,n+1) = mod(xr(:,n+1),1);  
            z(1:nParticles,n+1) = xr(:,n+1);
            z(nParticles+1:end,n+1) = r(:,1);           
    end
    %Save the results of the simulation into some rows of a 
    %matrix so that each block of nParticles rows is the results 
    %of a simulation.
    Results((1:nParticles)+nParticles*(i-1),:) = z(1:nParticles,:);
    
    %For the case where the noise strength is 0.25/sqrt(2pi)
    Results((1:nParticles)+nParticles*(i-1),:) = z(1:nParticles,oops);
end

%Plot histograms


%% Plot Section
%Sets the edges of the histogram boxes to be the x values
        timesteps = h*(0:N); %SDE Model plot times
        xhist = 0:0.01:1; %Define edges for histogram
        edges = xhist; %set the edges
        figure() %Plots the solution to the problem as time evolves
        plots = [1,21,41,61,81,101];%Plots for PDE solution

        %vector to save SDE results we use to plot
        Results_save = zeros(nParticles*nSim,length(plots));
        for iPlot = plots
            %Find which time step from the SDE model is
            %closest to the PDE plot time
            [~,I] = min(abs(plotTimes(iPlot) - timesteps));
           
            plot(x_p, rho(iPlot,:)') %Plot PDE solution
            hold on
            
            %counts the number of particles in each bin, 
            %and resaves the edges
            [bincounts,edges] = histcounts(Results(:,I),edges);
          
            %Finds the width of each bin
            dx = abs(edges(2)-edges(1));
             
            %Normalises the bincounts so that the total area
            %of the histogram is 1.
            bincounts = bincounts/(nParticles*nSim*dx); 

            Area = sum(bincounts*dx); %Check that the area of 
                                       %the histogram is 1.

            %Plots the histogram bin counts
            plot(edges(1:end-1)+dx/2,bincounts, 'o')             
            %save the results we used for the SDE histogram
            Results_save(:,s) = Results(:,I);
        end
        
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$\rho(x,t)$','Interpreter','latex')
          
          %Various legends
          
% legend('$t = 0$','$t = 0$', '$t = 1.2566$', '$t = 1.2566$',...
%        '$t = 2.5133$', '$t = 2.5133$', '$t = 3.7699$', ...
%         '$t = 3.7699$', '$t = 5.0265$',  '$t = 5.0265$', ...
%         '$t = 6.2832$', '$t = 6.2832$', 'Interpreter', 'latex')

%legend('$t = 0$', '$t = 6.2832$', '$t = 12.5664$',...
%       '$t = 15.7080$', '$t = 31.4159$', '$t = 47.1239$',...
%        '$t = 62.8319$', 'Interpreter', 'latex')

%legend('$t = 0$', '$t = 3.1416$', '$t = 6.2832$',...
%        '$t = 12.5664$','$t = 25.1327$', '$t = 37.6991$',...
%        '$t = 62.8319$','Interpreter', 'latex')

legend('$t = 0$','$t = 0$','$t = 6.2832$','$t = 6.2832$',...
'$t = 12.5664$','$t = 12.5664$', '$t = 25.1327$', '$t = 25.1327$',...
'$t = 50.2655$', '$t = 50.2655$', '$t = 75.3982$',...
'$t = 75.3982$','$t = 125.6637$','$t = 125.6637$',...
'Interpreter', 'latex')
      
clear Results
%Choose a suitable filename and save the data          
newfilename = 'PDE_SDE_Compare_s02';
save(newfilename)

%% PDE Functions
    function [drhodt] = rhs(t,rho)
        %This is the PDE Model with radicals
        drhodt = (sigma^2/2)*D2*rho + D*(rho.*(C*(rho+rho_rad(x_p))));
    end

    function Q = Qc(rho)
        %Calculates the order parameter
       Q = Int*(rho(end,:)'.*(C_Q*(rho(end,:)'))); 
    end

    function rhoIC = initialCondition(y)
         %Sets up the initial condition
        rhoIC = exp(-0.5*((y-mu)/Sigma).^2);
        Z = Int*rhoIC; %Integrates the initial condition
        rhoIC = rhoIC/Z;%divides the initial condition by its
                         %integral to normalise so all mass is
                         %in [0,1]
    
    end

    function rho_r = rho_rad(y)
        %Defines the radical distribution for the PDE model
             rho_r = exp(-0.5*((y-mu_r)/Sigma_r).^2);
             Z = Int*rho_r;
             %normalise rho_r so that the mass of the radicals
             %in the PDE model is the same proportion of the
             %normal mass in the PDE model and the SDE model
             rho_r = (nrParticles/nParticles)*rho_r/Z;
    end
%% SDE Functions
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
    Rij(mask) = Rij(mask) - sign(Rij(mask));
    
    nRij = abs(Rij)<=R;
    nRij = sum(nRij,1); %Finds the number of individuals who are
                        %within R of each individual in a row vector
    Rij(abs(Rij)>R) = 0; %Works out which opinions (i,j) aren't within R
                         %and sets all these Rij to be zero
    sumRij = sum(Rij,1); %Each column m is the sum of opinions in row m,
                         %which is the opinions that influence individual m
    F = (nParticles.^(-1) .* sumRij)';

 end
end