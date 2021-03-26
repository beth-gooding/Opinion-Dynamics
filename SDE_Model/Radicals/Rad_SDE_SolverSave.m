%Calculate an approximate solution to the SDE model including 
%radicals using the Euler Maruyama method. We need to call a 
%script that contains the variables defining number of 
%individuals, simulation end time, number of time steps, noise 
%strengths, confidence interval widths and initial condition.
%This script can calculate the dynamics both with and without
%radicals to allow us to see what difference the radical 
%influence makes.

clear all
close all
set(0, 'DefaultAxesFontSize', 30)
set(0, 'DefaultLineLineWidth', 2)
tic

%% Initial Conditions
InitSDERad; %Loading all of the initial inputs we need

Ones = ones(1, nParticles); %For the function F_SDE
W=randn(nParticles, N+1); %Generating realsiations of Brownian
                          %increments

%Setting the number of radicals as a proportion of the number
%of normals
% nrParticles = nParticles; %Equal numbers of radicals and
                            %normals
nrParticles = nParticles/5; %The radicals are one fifth of the
                            %size of the normals
r = zeros(nrParticles,N+1); %Setting the radicals opinions to
                            %zero at every time step

%For the function F_SDE when radicals are included                            
rOnes = ones(nParticles+nrParticles, 1); 


xr = x; %Vector for updating the normals opinion when radicals
        %are included
z = [xr;r]; %Vector for saving both the normal and radical 
            %opinions in so that the updating of the normals
            %opinions can be done easily.

%% Load exisiting simulation data or run a new one

%Generate potential filename
filename = [sprintf('rad3_data_%d_%d_%d_%d_%d_%d_%d_%d_%d_%s.mat',...
            T,N,nParticles,nrParticles,Rmin*100,Rmax*100,...
            bmin*100,bmax*100,c,IC)];
if exist(filename, 'file')==2 %If file exists load and plot data
    load(filename)
    disp('loaded data')
    
    figure()
    subplot(2,1,1) %Plot opinion paths of normals without
                   %radical influence
        for i = 1:nParticles
            plot(t,x(i,:),'.');
            hold on
            xlabel('Time'); ylabel('Opinion');
        end
    subplot(2,1,2) %Plot opinion paths of normals with
                   %radical influence
         for i = 1:nParticles
            plot(t,xr(i,:),'.');
            hold on
            xlabel('Time'); ylabel('Opinion');
        end
        ylim([0,1])
%       xlim([200,T]);

        figure() %Plot histograms of final distribution
                 %of normals' opinions without and with
                 %the radical influence
                 
        subplot(2,1,1) %Without radical influence
        histogram(x(:,N+1),50)
        axis([0,1,0,max(histcounts(x(:,N+1),50))+2])
        xlabel('Opinion'); ylabel('Count');

        subplot(2,1,2) %With radical influence
        histogram(xr(:,N+1),50)
        axis([0,1,0,max(histcounts(xr(:,N+1),50))+2])
        xlabel('Opinion'); ylabel('Count');
        
        
else %if data doesn't exist
    b = linspace(bmin,bmax,c); %vector of noise strengths
    R = linspace(Rmin,Rmax,c); %vector of CI widths
     Q = zeros(1,N+1); %Store the evolution of the order
                       %parameter without radical influence
     Qr = Q;           %Store the evolution of the order
                       %parameter with radical influence
    for j = 1:length(R) %For each confidence interval width
       
        %Fix initial value of the order parameter at 2R
        Q(1,1) = 2*R(j);
        Qr(1,1) = Q(1,1);
        
        for k = 1:length(b) %For each noise strength
            for n = 1:N
                
              %Update dynamics including the radicals
              [Fr, nrRij] = F_SDE(nParticles, z, R(j), n, rOnes);
              
              %Only use the first nParticles entries of Fr as
              %we only want to update the opinions of the 
              %normals
              xr(:,n+1) = xr(:,n) -h*Fr(1:nParticles)+...
                  b(k)*sqrt(h).*W(:,n);
              
              %Apply boundary condition
              xr(:,n+1) = mod(xr(:,n+1),1);
              
              %Update the xr entries of z to be used at the next
              %time step in F_SDE
              z(1:nParticles,n+1) = xr(:,n+1);
              
              %Calculate the new value of the order parameter
              Qr(1,n+1) = ((nParticles)^-2)*sum(nrRij(1:nParticles));
                
              %Calculate the dynamics without radical influence
              [F, nRij] = F_SDE(nParticles, x, R(j), n, Ones);
              x(:,n+1) = x(:,n) -h*F+ b(k)*sqrt(h).*W(:,n);
              
              %Apply boundary condition
              x(:,n+1) = mod(x(:,n+1),1);
              
              %Update the value of the order parameter
              Q(1,n+1) = (nParticles^-2)*sum(nRij);
            end
            
    t = 0:h:T; %Calculate the plot times
    figure() %Plot the opinion paths of the dynamics both
             %with and without the radical influence
    subplot(2,1,1) %Without the radical influence
        for i = 1:nParticles
            plot(t,x(i,:),'.');
            hold on
            xlabel('Time'); ylabel('Opinion');
        end
    subplot(2,1,2) %With the radical influence
         for i = 1:nParticles
            plot(t,xr(i,:),'.');
            hold on
            xlabel('Time'); ylabel('Opinion');
        end
        ylim([0,1])
%       xlim([200,T]);
        
        %Plot the evolution of the order parameter over time
        %both with and without radicals
        figure()
        plot(t,Q, '.r')
        xlabel('Time');ylabel('$Q$', 'Interpreter','latex');
        hold on
        plot(t,Qr, '.k')
        legend('Without Radicals', 'With Radicals')
    
        %Plot the histograms of the final distributions of 
        %opinions both with and without the radicals
        figure() 
        subplot(2,1,1) %without radical influence
        histogram(x(:,N+1),50)
        axis([0,1,0,max(histcounts(x(:,N+1),50))+2])
        xlabel('Opinion'); ylabel('Count');

        subplot(2,1,2) %with radical influence
        histogram(xr(:,N+1),50)
        axis([0,1,0,max(histcounts(xr(:,N+1),50))+2])
        xlabel('Opinion'); ylabel('Count');
        end
    end
end
    
%% Save the results if we run a new simulation
newfilename = sprintf('rad3_data_%d_%d_%d_%d_%d_%d_%d_%d_%d_%s',...
                  T,N,nParticles,nrParticles,Rmin*100,Rmax*100,...
                  bmin*100,bmax*100,c,IC);
time = toc;
save(newfilename)

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
    Rij(mask) = Rij(mask) - sign(Rij(mask));
    
    nRij = abs(Rij(1:nParticles,1:nParticles))<=R;
    nRij = sum(nRij,1); %Finds the number of individuals who are
                        %within R of each individual in a row vector
    Rij(abs(Rij)>R) = 0; %Works out which opinions (i,j) aren't within R
                         %and sets all these Rij to be zero
    sumRij = sum(Rij,1); %Each column m is the sum of opinions in row m,
                         %which is the opinions that influence individual m
    F = (nParticles.^(-1) .* sumRij)';

 end