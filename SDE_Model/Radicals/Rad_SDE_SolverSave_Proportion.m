%% Rad_SDE_SolverSave_Proportion.m
%This script calculates the number of normals who are
%radicalised for different confidence interval widths and 
%numbers of radicals, giving the results in the form of a 
%colour map. We need to call a script that contains the variables 
%defining number of individuals, simulation end time, number of 
%time steps, noise strengths, confidence interval widths and
%initial condition.

clear all
close all
set(0, 'DefaultAxesFontSize', 30)
set(0, 'DefaultLineLineWidth', 2)
tic

%% Initial Conditions
InitSDERadExperiment; %Loading the initial inputs we need

Ones = ones(1, nParticles); %For the function F_SDE
w=randn(nParticles, N+1); %Realisations of Brownian 
                          %increments
b = bmin; %Noise Strength
R = linspace(Rmin,Rmax,50); %Confidence Interval Widths

%Generating three further distinct realisations of Brownian
%increments
rng(1)
W1=randn(nParticles, N+1);
rng(2)
W2=randn(nParticles, N+1);
rng(3)
W3=randn(nParticles, N+1);

nrParticles = linspace(1,50,50); %Number of radicals
nSim = 4; %Number of different realisations of noise

%For saving the number of normals who end up with the radical opinion
N2R = zeros(length(R), length(nrParticles));
R_r = 0.02; %What we use to define close to the radical opinion.

%% Run the experiment

 for m = 1:length(nrParticles) %Run for each number of radicals
   
    r = zeros(nrParticles(m),N+1);  %Setting the radical opinion to
                                    %be zero
    
    %For the function F_SDE_Rad                               
    rOnes = ones(nParticles+nrParticles(m), 1);

    xr = x; %For saving the normals opinion in
    z = [xr;r]; %For updating the opinions of the normals
    
    
     for j = 1:length(R) %for each width
            
         %for saving the number of radicalised normals
         %for each realisation of Brownian increments
         nRad = zeros(nSim,1); 
            
         for k = 1:nSim %choose the realisation of Brownian
                           %increments
             if k == 1
                 W = W1;
             elseif k == 2
                 W = W2;
             elseif k == 3
                 W = W3;
             elseif k == 4
                 W = w;
             end
             for n = 1:N %for each time step
                 
                 %Calculates the forcing and the number of
                 %radicalised normals at the previous time step
                 [Fr, ~, nRij_Rad] = F_SDE_Rad(nParticles, z,...
                                            R(j), n, rOnes, R_r);
                 %Update the normals opinions using Euler
                 %Maruyama method
                 xr(:,n+1) = xr(:,n) -h*Fr(1:nParticles)+...
                                b*sqrt(h).*W(:,n);
                 %Apply the boundary condition          
                 xr(:,n+1) = mod(xr(:,n+1),1); 
                 %Update opinions to be used in the next timestep
                 z(1:nParticles,n+1) = xr(:,n+1);
             end
             %Calculate the final number of radicalised normals
             [~,~,nRij_Rad] = F_SDE_Rad(nParticles, z, R(j),...
                                        N+1, rOnes, R_r);
             nRad(k) = nRij_Rad;
                   
         end
         %Calculate the average number of radicalised normals
         %to be plotted in the colour map
         N2R(j,m) = mean(nRad);
    end  
 end

%% Plotting and saving data
%Plot the results of the experiment using a colour map
figure()
pcolor(R,nrParticles,N2R');
colorbar
xlabel('$R$','Interpreter','latex');
ylabel('Number of Radicals');

time = toc;
save('uniformspread_radical_experiment_2')

%% Functions
 function [F,nRij, nRij_Rad] = F_SDE_Rad(nParticles,x,R,n, Ones, R_r)
 %F_SDE Calculates F for SDE solver
 
    xn = x(:,n)'; %Extracts all individuals opinions at time n in a row
    yn = xn';     %Column vector of all the opinions at time n
    X = xn(Ones,:); %Makes a square matrix where each column m is a load of
                    %copies of the opinion of individual m.
    Y = yn(:,Ones); %Square matrix, each row m is copies of the opinion of
                    %individual m.
                    
    Z = ones(size(xn)); %Calculates ones for the radicals
    
    Rij = X-Y;      %Calculates the difference in the opinions of all the
                    %individuals with the sign
    mask = ( abs(Rij) > 1/2 ); % deal with periodicity
    Rij(mask) = Rij(mask) - sign(Rij(mask)); %Makes sure distance in the
                                             %closest distance i.e in 
                                             % [-0.5,0.5].
    %Finds the pairs of normals which are within R of each other
    nRij = abs(Rij(1:nParticles,1:nParticles))<=R; 
    nRij = sum(nRij,1); %Finds the number of individuals who are
                        %within R of each individual in a row vector
    Rij(abs(Rij)>R) = 0; %Works out which opinions (i,j) aren't within R
                         %and sets all these Rij to be zero
    sumRij = sum(Rij,1); %Each column m is the sum of opinions in row m,
                         %which is the opinions that influence individual m
    
    F = (nParticles.^(-1) .* sumRij)'; %Calculates the forcing for the SDE
    
    %Calculates the distance between radicals and normals
    Rij_Rad = xn-Z; 
    
    %Checking the distances to make sure they lie in [-0.5,0.5]
    mask_Rad = (abs(Rij_Rad) > 1/2); 
    Rij_Rad(mask_Rad) = Rij_Rad(mask_Rad) - sign(Rij_Rad(mask_Rad));
    
    %finding which normals are close to the radical opinion
    nRij_Rad = abs(Rij_Rad(1:50))<=R_r; %50 is number of normals 
    
    nRij_Rad = sum(nRij_Rad); %Finds the total number of normals 
                              %close to the radical opinion
 end