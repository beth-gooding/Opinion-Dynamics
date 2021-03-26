%% Time_Comparison.m 
%We use this script to compare the efficiency of algorithms 1 and 2, to see
%which one is better to use when we implement our full SDE model.

clear all
close all

%% Parameters
nParticles = 2.^(2:10); %Number of individuals
T = 50; %End time of simulation
N = 10^3; %Number of time steps
h = T/N; %Size of timestep

R = 0.2; %Half the width of the confidence intervals
b= 0.02; %Strength of the noise

%Preallocate vectors to save simulation run times in
time  = zeros(length(nParticles),1);
time_s = zeros(100,1);
mean_time_s = time;
for m = 1:length(nParticles)
    m;
%% Original SDE Calculations
    W=randn(N+1, nParticles(m)); %Realisations of brownian increments
    W_s = W';   %Take the transpose due to different set ups for the  
                %two calculations of the forcing term, this is for the
                %improved solver
    
    x = zeros(N+1,nParticles(m)); %To save the opinion paths in
    x(1,:) = rand(1,nParticles(m)); %Sample the initial condition from the
                                    %uniform distribution
    x_s = x'; %Transpose of the paths vector to be used in the improved
              %solver.
    
    tic %start timing the length of the simulation
    for n = 1:N %for each time step

        for i = 1:nParticles(m) %for each individual's opinion
            J = []; %List to save close opinions in
            J2 = []; %Second list to save close opinions in
            for j = 1:nParticles(m) %for every opinion check whether they
                                    %are close to individual i
                                    
                %First make sure we are using the shortest distance
                if abs(x(n,i)-x(n,j)) > 0.5
                    %adjust the distance so it lies in [-0.5,0.5]
                    sm = x(n,i)-x(n,j)-sign(x(n,i)-x(n,j));
                    if abs(sm)<=R %now if opinion i and j are close
                        J2 = [J2, x(n,j)]; %add to our list
                    end
                elseif abs(x(n,i) - x(n,j)) <= R %if i and j are close
                    J = [J,x(n,j)]; %add to our list
                end
                %calculate sum for the forcing term making sure to use the
                %shortest distance for the second list of opinions
                Jsum = sum(x(n,i)-J);
                J2sum = sum(x(n,i)-J2-sign(x(n,i)-J2));
                %calculate the full forcing term
                F = (1/nParticles(m))*(Jsum+J2sum);
            end

            %update the opinions using the Euler-Maruyama method
            x(n+1,i) = x(n,i) -h*F + b*sqrt(h)*W(n,i);

            x(n+1,:) = mod(x(n+1,:),1); %Apply the boundary condition    

        end
    end
    time(m) = toc; %save the simulation run time

%% Improved SDE Calculations
%Now run a simulation using exactly the same set up, but with the improved
%forcing calculation.

 Ones = ones(1, nParticles(m)); %For use in the F_SDE function

    for k = 1:100 %run 100 simulations for each number of individuals
    tic %start timing one simulation
        for n = 1:N %for each time step
            %calculate the forcing term using F_SDE
           [F, nRij] = F_SDE(nParticles(m), x_s, R, n, Ones);

           %update the opinions using the Euler-Maruyama method
           x_s(:,n+1) = x_s(:,n) -h*F + b*sqrt(h).*W_s(:,n);
           
           %Apply the boundary condition
           x_s(:,n+1) = mod(x_s(:,n+1),1);  
        end
    time_s(k) = toc; %record the run time of this simulation
    end
    mean_time_s(m) = mean(time_s); %calculate the mean run time for
                                   %each number of individuals
    
%A visual check to make sure that the two methods are giving us the same
%results, which they should be as none of the set up is different.
%     %Plotting paths to check they match
%     t = 0:h:T;
%     figure()
%     for i = 1:nParticles(m)
%         plot(t,x(:,i),'.');
%         hold on
%         xlabel('Time'); ylabel('Opinion');
%     end
% 
%     figure()
%     for i = 1:nParticles(m)
%         plot(t,x_s(i,:),'.');
%         hold on
%         xlabel('Time'); ylabel('Opinion');
%     end
end
%% Comparison of Time Taken
%Plots a loglog plot which allows us to compare the time taken to run the
%simulations which allows us to see the the power relationship between the
%simulation run time and the number of individuals.
figure()
loglog(nParticles, time, 'k.', 'markersize', 20)
hold on
loglog(nParticles, mean_time_s, 'r.', 'markersize', 20)
xlabel('Number of Individuals'); 
ylabel('Mean Simulation Run time in Seconds');
legend('Original Solver', 'Improved Solver');



%% Functions
 function [F,nRij] = F_SDE(nParticles,x,R,n,Ones)
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