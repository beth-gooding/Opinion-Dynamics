
%Calculate an approximate solution to the full SDE model using the
% Euler-Maruyama method. We need to call a script that contains 
%the variables defining number of individuals, simulation end time, 
%number of time steps, noise strengths, confidence interval widths
%and initial condition.
clear all
close all
set(0, 'DefaultAxesFontSize', 30)
set(0, 'DefaultLineLineWidth', 2)

tic

%% Initial Conditions
InitSDE; % For parameter experiment sampling from a uniform
         %distribution
%InitSDEWang; % For parameter experiment replicating Wang Figure 3

%InitSDEFig4A; %Replicating Wang Figure 4A 
%InitSDEFig4B; %Replicating Wang Figure 4B 
%InitSDEFig4C; %Replicating Wang Figure 4C 
%InitSDEFig4D; %Replicating Wang Figure 4D 
%InitSDEFig4E; %Replicating Wang Figure 4E
%InitSDEFig4F; %Replicating Wang Figure 4F

%InitSDEWangFig2; %Replicating Wang Figure 2

%Init_ODE; %For simulating the ODE model

Ones = ones(1, nParticles); %For use in F_SDE
W=randn(nParticles, N+1); %Generating realisations of Brownian 
                          %increments

%% Load existing simulation or run a new one
%Generate potential filename
filename = [sprintf('data_%d_%d_%d_%d_%d_%d_%d_%d_%s.mat',... 
        T,N,nParticles,Rmin*100,Rmax*100,bmin*100,bmax*100,c,IC)];
%filename = [sprintf('ode_data_%d_%d_%d_%d_%d_%d_%d_%d_%d_%s.mat',...
       %T,N,nParticles,Rmin*100,Rmax*100,bmin*100,bmax*100,c,ra,IC)];

if exist(filename, 'file')==2 %If file exists load and plot data
    load(filename)
    disp('loaded data')
    t = 0:h:T;
         
        figure()
        subplot(2,1,1) %Plot the paths of opinions
        for i = 1:nParticles %plots all opinions
        %for i = 1:3:nParticles %Plots one third for the Wang figures
           plot(t,x(i,:),'.');
            hold on
            xlabel('Time'); 
            ylabel('Opinion');
        end   
%       xlim([0,500]);
        xlim([200,500]);
        xlim([200,400]);
        
        subplot(2,1,2) %Plot the histogram of the final distribution
                       %of opinions
        histogram(x(:,N+1),50) %Use 50 bins
        axis([0,1,0,max(histcounts(x(:,N+1),50))+2])
        xlabel('Opinion'); 
        ylabel('Number of Individuals');
        
%         figure() %Plot the evolution of the order parameters
%          %subplot(2,1,2)
%          for d = 1:3
%          plot(t,Q(d,:))
%          xlabel('Time');ylabel('$Q$', 'Interpreter','latex');
%          %ylim([0,0.6])
%          hold on 
%          end
    
%     figure() %Plot the colour map of the parameter experiment
%     pcolor(R,b,Q'); shading interp;
%     colorbar
%     xlabel('$R$','Interpreter','latex');
%      ylabel('$\sigma$', 'Interpreter','latex');
%     

else %If data doesn't exist, run simulation
    b = linspace(bmin,bmax,c); %create vector of all noise strengths
    R = linspace(Rmin,Rmax,c); %create vector of confidence intervals
    Q = zeros(length(R),length(b)); %for saving the order parameter 
                                    %values
    
    %Use for calculating the order parameter for multiple simulations 
    %using the same parameter pair.
%     Q = zeros(5,N+1);
%     for a = 1:5
%         rng(a)
%         W=randn(nParticles, N+1);
    for j = 1:length(R) %For all lengths of the confidence intervals
        for k = 1:length(b) %For all values of noise strength
%             
                for n = 1:N %For each time step
                  %Calculate the forcing and the number of pairs of
                  %opinions that are close to each other
                  [F, nRij] = F_SDE(nParticles, x, R(j), n, Ones);
                    
                  %Update the values of opinions using Euler-Maruyama
                  x(:,n+1) = x(:,n) -h*F + b(k)*sqrt(h).*W(:,n);
                    
                  %Apply periodic boundary conditions
                  x(:,n+1) = mod(x(:,n+1),1);
                  %Q(a,n+1) = (nParticles^-2)*sum(nRij);
                end
        end
        %Calculate the value of the order parameter at the final time
        Q(j,k) = (nParticles^-2)*sum(nRij);
         
        t = 0:h:T; %calculate plot times
%         figure() %plot the opinion paths
%         subplot(1,2,1)
%         for i = 1:nParticles %All opinions
%          %for i = 1:3:nParticles %Every third opinion
%            plot(t,x(i,:),'.');
%             hold on
%             xlabel('Time'); ylabel('Opinion');
%         end   
%         %xlim([0,500])

%         subplot(1,2,2) %plot the histogram of the distribution of
%                        %opinions
%         histogram(x(:,N+1),50)
%         axis([0,1,0,max(histcounts(x(:,N+1),50))+2])
%         xlabel('Opinion'); ylabel('Number of Individuals');
%         
%         figure() %Plot the evolution of the order parameter in time
%         plot(t,Q, '.k')
%         xlabel('Time');ylabel('$Q$', 'Interpreter','latex');

        %end
%     end
   figure() %Plot the colour map of the order parameter
   pcolor(R,b,Q'); shading interp;
   colorbar
   xlabel('$R$','Interpreter','latex');
   ylabel('$b$', 'Interpreter','latex');

    end
%% Save the results if we run a new simulation
%newfilename = sprintf('data_%d_%d_%d_%d_%d_%d_%d_%d_%s.mat',...
        %T,N,nParticles,Rmin*100,Rmax*100,bmin*100,bmax*100,c,IC);
newfilename=[sprintf('ode_data_%d_%d_%d_%d_%d_%d_%d_%d_%d_%s.mat',... 
     T,N,nParticles,Rmin*100,Rmax*100,bmin*100,bmax*100,c,ra,IC)];
time = toc;
save(newfilename)

end

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
    mask = ( abs(Rij) > 1/2 ); %Find the distances which are not the
                               %shortest distances
                               
    %Adjust the distances so that they are in [-0.5,0.5]
    Rij(mask) = Rij(mask) - sign(Rij(mask)); 
    
    nRij = abs(Rij)<=R; %Logical matrix with one in (i,j) if |xj-xi|<=R
    nRij = sum(nRij,1); %Finds the number of individuals who are
                        %within R of each individual in a row vector
    Rij(abs(Rij)>R) = 0; %Works out which opinions (i,j) aren't within R
                         %and sets all these Rij to be zero
    sumRij = sum(Rij,1); %Each column m is the sum of opinions in row m,
                         %which is the opinions that influence individual m
    %Calculate the forcing term that is applied to all opinions
    F = (nParticles.^(-1) .* sumRij)';
 end
