function PDE_Diffusion()
 %Solves the diffusion equation and also plots an exact solution
 %to find the accuracy of our numerical solution.
%% Initial Setup
clear all
close all

set(0, 'DefaultAxesFontSize', 30);
set(0, 'DefaultLineLineWidth', 3);

N = 100; %Number of equispaced grid points
h1 = 1/N; %Spacing for grid points on [0,1]
hD = 2*pi/N; %Spacing for grid points on [0,2pi]
x1 = (0:N-1)'*h1; %grid points for evaluating functions at
xD = (0:N-1)'*hD; %grid points used to construct the spectral
                  %differentiation matrix

                
% Construct spectral differentiation matrix using code from P2 by 
% Trefethen:
column = [0 .5*(-1).^(1:N-1).*cot((1:N-1)*hD/2)];
D = toeplitz(column,column([1 N:-1:2])); %First derivative
D2 = D^2; %Second derivative

%% Solve ODE

plotTimes=[0,0.5,1] %Sets which times we would like the approximate
                    %solution to be calculated at

%Calculates the exact solution at t = 0, which is also the initial
%condition for our numerical approximation.
rhoIC1 = exactSolution1(x1,0); 

%Solves the resulting ODE with respect to time, giving us the solution
%evaluated at the plot times we asked for.
[plotTimes1,rho1] = ode15s(@rhs, plotTimes, rhoIC1);


figure()
%plot the exact and approximate solutions at the three plot times
for iPlot = 1:length(plotTimes1) % on [0,1]

    plot(x1, rho1(iPlot,:)') %Plot the approximate solution
    hold on
    %Calculate the exact solution at the plot times
    rhoE(iPlot,:) = exactSolution1(x1,plotTimes1(iPlot));
    if iPlot == 1
        plot(x1, rhoE,'--k');
    elseif iPlot == 2
        plot(x1, rhoE,':k');
    elseif iPlot == 3
        plot(x1, rhoE,'-.k');
    end
    ylim([0.9,3]);
    xlabel('$x$', 'Interpreter', 'latex'); 
    ylabel('$\rho(x,t)$', 'Interpreter', 'latex')
    legend('$\hat{\rho}(x,0)$','$\rho(x,0)$',...
        '$\hat{\rho}(x,0.5)$', '$\rho(x,0.5)$', ...
        '$\hat{\rho}(x,1)$','$\rho(x,1)$','Interpreter',...
        'latex')
end
%maximum error at t = 0.5
error05 = norm(abs(rho1(2,:)-rhoE(2,:)), inf);
%maximum error at t = 1
error1 = norm(abs(rho1(3,:)-rhoE(3,:)), inf);

%% Function

%Defines the right hand side of our PDE which reduces to an ODE with
%respect to time
function [drhodt] = rhs(t,rho)
    drhodt = D2*rho;
end

function rhoE = exactSolution1(y,t) %Exact solution on [0,1]
  rhoE = exp(-t)*sin(2*pi*y) + 2;
end

end