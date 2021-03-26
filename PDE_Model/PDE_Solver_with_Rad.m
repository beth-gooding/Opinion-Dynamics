function PDE_Solver_with_Rad()
%This script implements the PDE model for opinion dynamics 
%including a radical influence. The boundary condition
%is a periodic boundary condition over the interval [0,1],
%which is automatically implemented as all the parts 
%used are designed to work with periodic functions
%% Initial Setup
clear all
close all

N = 300; %Number of grid points
h = 1/N; %Stepsize of periodic points
x = (0:N-1)'*h; %Equispaced periodic points on [0,1]
x_rad = x; %Equispaced points to evaluate the radicals at

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

%% Width of confidence intervals
Rfull = 0.2; %For replicating the Wang results
%Rfull = linspace(0.01,0.37, 30);
%Rfull = Rfull(27);

%For saving the order parameter in for parameter experiment
%Q = zeros(length(sigma2),length(Rfull));
a = 0.5; %Used for calculating the total mass of the radical 
       %distribution, as a proportion of the mass of the
       %normal distribution. E.g a = 1 means the mass of
       %radicals and normals is equal, whilst a = 0.5
       %means that the mass of the radicals is half the 
       %mass of the normals.
%% Solve ODE

for i = 1:length(Rfull) %For each confidence interval width
     R = Rfull(i); %set a particular width
    
    for j = 1:length(sigmaw) %For each noise strength
     
        sigma = sigmaw(j); %Set a particular noise strength
        
        %Calculates the matrix that will do integration and 
        %interpolation
        C = ConvolutionMatrix(x,R); 
        
        %Needed to calculate Qc, the order %parameter  
        C_Q = ConvolutionMatrixQ(x,R); 
        
        plotTimes=[0:T/100:T]; %Times to plot our solution on

        %Creates the initial condition for our problem
        rhoIC = initialCondition(x); 
        
        %Calculate the solution without the radicals
        %Solves the resulting ODE now that we are able to
        %approximate the right hand side of the PDE which
        %deals with all the derivatives with respect to x
        [plotTimes,rho] = ode15s(@rhs, plotTimes, rhoIC); 
        
        
        %Calculate the solution with the radicals
        %Solves the resulting ODE now that we are able to
        %approximate the right hand side of the PDE which
        %deals with all the derivatives with respect to x
        [plotTimes,rho_Rad] = ode15s(@rhs_rad, plotTimes, rhoIC);                                                  
                                                          
                                             
        Q(i,j) = Qc(rho_Rad) %Qc for model with radicals
        
        %% Plotting
        %Various different plot times depending on the final time
        plots = [1,21,41,61,81,101]; %Short Time
%       plots = [1,6,11,21,41,61,101]; %Short Time other ones
%       plots = [1,26,51,76,101]; %Long Time

        figure() %Plot the solutions as time evolves
        subplot(1,2,1) %Plot without radicals
        for iPlot = plots
%         for iPlot = 1:length(plotTimes)

            %hold off
            plot(x, rho(iPlot,:)')
    
            ylim([0,18]);
            pause(0.1)
            hold on
    
        end
          xlabel('$x$', 'Interpreter', 'latex');
          ylabel('$\rho(x,t)$','Interpreter','latex')
          
%Legends for different plot times

%legend('$t = 0$', '$t = 1.2566$', '$t = 2.5133$', '$t = 3.7699$',...
%         '$t = 5.0265$', '$t = 6.2832$', 'Interpreter', 'latex')
%legend('$t = 0$', '$t = 3.1416$', '$t = 6.2832$', '$t = 12.5664$',...
%         '$t = 25.1327$', '$t = 37.6991$', '$t = 62.8319$',...
%         'Interpreter', 'latex')
legend('$t = 0$','$t = 6.2832$','$t = 12.5664$',...
         '$t = 25.1327$','$t = 50.2655$','$t = 75.3982$',...
              '$t = 125.6637$','Interpreter', 'latex')
             
        subplot(1,2,2) %Plot with radical influence
%         for iPlot = 1:length(plotTimes)
        for iPlot = plots
%             hold off
            plot(x, rho_Rad(iPlot,:)')
    
            ylim([0,18]);
            pause(0.5)
            hold on
    
        end
        xline(0.5, 'k--');
        xline(0.8,'k--');
        xlabel('$x$', 'Interpreter', 'latex');
        ylabel('$\rho(x,t)$','Interpreter','latex')
        
 %Various legends for different plot times
%legend('$t = 0$', '$t = 1.2566$', '$t = 2.5133$', '$t = 3.7699$',...
%         '$t = 5.0265$', '$t = 6.2832$', 'Interpreter', 'latex')
%legend('$t = 0$', '$t = 3.1416$', '$t = 6.2832$', '$t = 12.5664$',...
%         '$t = 25.1327$', '$t = 37.6991$', '$t = 62.8319$',...
%         'Interpreter', 'latex')
legend('$t = 0$','$t = 6.2832$','$t = 12.5664$',...
         '$t = 25.1327$','$t = 50.2655$','$t = 75.3982$',...
              '$t = 125.6637$','Interpreter', 'latex')
% legend('$t = 0$', '$t = 3.1416$', '$t = 6.2832$',...
%         '$t = 12.5664$','$t = 25.1327$',...
%          '$t = 37.6991$', '$t = 62.8319$',...
%         'Interpreter', 'latex')   


    end
end
 %Choose a suitable filename and save the data for the simulation
 newfilename = sprintf('PDE_With_Rad_4')
 save(newfilename)
%% Functions
 
    function [drhodt] = rhs(t,rho) 
        %This defines our PDE model without a radical influence
        drhodt = (sigma^2/2)*D2*rho + D*(rho.*(C*(rho)));
    end

function [drhodt] = rhs_rad(t,rho)
    %This defines our PDE model with a radical influence
        drhodt = (sigma^2/2)*D2*rho +...
            D*(rho.*(C*(rho+rho_rad(x_rad))));
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

    function rho_r = rho_rad(y)
        %Sets up the distribution of radicals
         rho_r = exp((-1000*(y-0.8).^2)); %Initial distribution
         Z = Int*rho_r; %Integrates the initial distribution
         rho_r = a*rho_r/Z; %divides the initial distribution
                            %by its integral over [0,1] so that
                            %all of the mass is in [0,1], and 
                            %multiplies by a so that the final
                            %mass of the distribution is the
                            %desired proportion of the normals
                            %mass
    end


end