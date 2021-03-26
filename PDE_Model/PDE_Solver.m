function PDE_Solver()
%This script implements the PDE model for opinion dynamics. The 
%boundary condition is a periodic boundary condition over the 
%interval [0,1], which is automatically implemented as all the
%parts used are designed to work with periodic functions
%% Initial Setup
clear all
close all

N = 300; %Number of grid points
h = 1/N; %Stepsize of periodic points
x = (0:N-1)'*h; %Equispaced periodic points on [0,1]


% Construct spectral differentiation matrix on [0,2*pi]:
hD = 2*pi/N;
column = [0 .5*(-1).^(1:N-1).*cot((1:N-1)*hD/2)];
D = toeplitz(column,column([1 N:-1:2])); %First derivative

% D on [0,1]
D = 2*pi*D; %Transform D so that it is differentiating on [0,1]

D2 = D^2; %Second derivative, found by squaring D.

Int = h*ones(1,N); % integration weights on [0,1], used in the
                   %calculation of the order parameter and for
                   %normalising the initial condition.

%% Noise Parameters

%For the set up given in the Wang paper
%sigma1 = 0.05;
sigma1 = [0.1, 0.15,0.2, 0.25];

%To actually replicate the Wang results
sigmaw = sigma1/sqrt(2*pi);

%For further investigations
%sigma2 = linspace(0.05, 0.25, 30)/sqrt(2*pi);
%sigma2 = sigma2(2);

%% Width of confidence intervals
Rfull = 0.2; %For replicating the Wang results
%Rfull = linspace(0.01,0.37, 30);
%Rfull = Rfull(27);

%For saving the order parameter in for parameter experiment
%Q = zeros(length(sigma2),length(Rfull));
%% Solve ODE
for i = 1:length(sigmaw) %for each noise strength
                         %change as appropriate

%Choose a specific noise strength                         
%     sigma = sigma2(i);
%     sigma = sigma1(i);
    sigma = sigmaw(i);

    for j = 1:length(Rfull) %for each confidence interval width
        
        R = Rfull(j); %choose specific width
        
        %Calculates the matrix that will do integration and 
        %interpolation
        C = ConvolutionMatrix(x,R); 
        
        %Needed to calculate Qc, the order %parameter  
        C_Q = ConvolutionMatrixQ(x,R); 
        
        plotTimes=[0:T/100:T]; %Times to plot our solution on

        %Creates the initial condition for our problem
        rhoIC = initialCondition(x); 
        
        %Solves the resulting ODE now that we are able to
        %approximate the right hand side of the PDE which
        %deals with all the derivatives with respect to x
        [plotTimes,rho] = ode15s(@rhs, plotTimes, rhoIC); 
        
        Q(i,j) = Qc(rho); %Calculates the order parameter
        
        
        figure() %Plots the solution to the problem as time evolves
         for iPlot = [1,21,41,61,81,101]
%           for iPlot = 1:length(plotTimes)
            
            %hold off
            plot(x, rho(iPlot,:)')
%             plot(x,rho(end,:)')
            hold on
            ylim([0,14]);
            pause(0.1)
    
         end
          xlabel('$x$', 'Interpreter', 'latex');
          ylabel('$\rho(x,t)$','Interpreter','latex')
          legend('$t = 0$', '$t = 1.2566$', '$t = 2.5133$',...
              '$t = 3.7699$','$t = 5.0265$', '$t = 6.2832$',...
                  'Interpreter', 'latex')
    end
end

%choose a suitable filename to save using
newfilename = sprintf('parameter_plot_particular_3');
save(newfilename)

%For plotting the results of the parameter experiment
% figure()
% pcolor(Rfull,sigma2,Q); shading interp;
%% Functions

    
    function [drhodt] = rhs(t,rho)
        %This defines our PDE model
        drhodt = (sigma^2/2)*D2*rho + D*(rho.*(C*rho));
    end

    function Q = Qc(rho) 
        %Calculates the order parameter
       Q = Int*(rho(end,:)'.*(C_Q*(rho(end,:)'))); 
    end

    function rhoIC = initialCondition(y)
        %Sets up the initial condition
        rhoIC = exp(-20*(y-0.5).^2); %Wang's initial condition
        Z = Int*rhoIC; %Integrates the initial condition
        rhoIC = rhoIC/Z; %divides the initial condition by its
                         %integral to normalise so all mass is
                         %in [0,1]
    end


end