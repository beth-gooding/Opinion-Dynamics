
%% This script illustrates a range of pseudospectral methods


%% Finite Difference Scheme vs Spectral Differentition
%Made using the files P1.m and P2.m by Trefethen
%We differentiate the function u = (sin(x)^2)cos(x) over [0,2pi]
%using both a fourth order finite difference scheme and 
%spectral differentiation and compare the accuracy of both
%methods

%Fourth order finite difference scheme
Nvec = 2.^(3:12); %A vector of the number of grid points
  figure()
  for N = Nvec %for all the numbers of grid points
    h = 2*pi/N; %Set the distance between gridpoints
    x = (1:N)'*h; %Define the grid points
    %Evaluate the function we want to differentiate
    u = (sin(x).^2).*cos(x);
    %Evaluate the exact derivative of u
    uprime = 2.*sin(x).*cos(x).^2 - sin(x).^3;

    % Construct sparse fourth-order differentiation matrix:
    e = ones(N,1);
    D =   sparse(1:N,[2:N 1],2*e/3,N,N)...
        - sparse(1:N,[3:N 1 2],e/12,N,N);
    D = (D-D')/h;

    %Approximate the derivative of u by calculating D*u
    %Calculate the maximum error in the derivative for
    %each number of grid points
    error = norm(D*u-uprime,inf);
    subplot(1,2,1) %Plot the maximum errors on a loglog plot
    loglog(N,error,'k.','markersize',15), hold on
  end
  grid on, xlabel N, ylabel error
  title('Convergence of fourth-order finite differences')
  %plot the line with gradient -4 for reference as the errors 
  %of the fourth order finite difference scheme should also 
  %decrease with that gradient
  semilogy(Nvec,Nvec.^(-4),'--') 
  text(105,5e-8,'N^{-4}','fontsize',18)

  %Spectral Differentiation
  for N = 2:2:100 %For various numbers of grid points
    h = 2*pi/N; %Define the distance between grid points
    x =(1:N)'*h; %Define each of the grid points
   
    % Construct spectral differentiation matrix:
    column = [0 .5*(-1).^(1:N-1).*cot((1:N-1)*h/2)];
    D_s = toeplitz(column,column([1 N:-1:2]));

    %Approximate the derivative of u by calculating D_s*u
    %Calculate the maximum error in the derivative for
    %each number of grid points
    error = norm(D_s*u-uprime,inf);
    subplot(1,2,2) %plot the maximum errors on a loglog plot
    loglog(N,error,'k.','markersize',15), hold on
    
  end
  grid on, xlabel N, ylabel error
  title('Convergence of spectral differentiation')

  %% Accuracy of Spectral Differentiation
 %This script explores the accuracy of spectral differentiation for 
 %smooth and non smooth functions. This script is made using the file
 %P4.m by Trefethen. We differentiate two functions over [0,2pi].

% Set up grid:
  N = 24; %Number of grid points 
  h = 2*pi/N; %Grid spacing
  x = h*(1:N)'; %Define the grid points
  
  %Set up the spectral differentiation matrix:
  column = [0 .5*(-1).^(1:N-1).*cot((1:N-1)*h/2)]';
  D = toeplitz(column,column([1 N:-1:2]));

  %Differentiation of |x-pi|, not smooth
  v = abs(x-pi); %Define the function
  %Plot the function v
  subplot(2,2,1), plot(x,v,'.-','markersize',13)
  axis([0 2*pi -.5 3.5]), grid on, title('function')
  xlabel('$x$', 'interpreter', 'latex');
  ylabel('$f(x)$', 'interpreter', 'latex')
  
  %Approximate the derivative of v by calculating D*v and
  %plot this approximation
  subplot(2,2,2), plot(x,D*v,'.-','markersize',13)
  axis([0 2*pi -2 2]), grid on, title('Spectral Derivative')
  xlabel('$x$', 'interpreter', 'latex'); 
  ylabel('Derivative');

  
  %Differentiation of (sin(x)^2)*(cos(x)), smooth
  v = (sin(x).^2).*cos(x); %Evaluate the function at the grid points
  %Evaluate the exact derivative at the grid points
  vprime = 2.*sin(x).*cos(x).^2 - sin(x).^3;
  
  %Plot the function v
  subplot(2,2,3), plot(x,v,'.-','markersize',13)
  axis([0 2*pi -1 1]), grid on
  xlabel('$x$', 'interpreter', 'latex'); 
  ylabel('$g(x)$', 'interpreter', 'latex')
  
  %Approximate the derivative using D*v and plot
  subplot(2,2,4), plot(x,D*v,'.-','markersize',13)
  axis([0 2*pi -2 2]), grid on
  xlabel('$x$', 'interpreter', 'latex'); 
  ylabel('Derivative');
  
  %Calculate the maximum error in the approximate derivative
  max_error = norm(D*v-vprime,inf);
  text(2.2,1.4,['max error = ' num2str(max_error)])
  
  %% Second Derivative Check
  %This script checks that the square of the spectral 
  %differentiation matrix does indeed calculate the second 
  %derivative and that this calculation is accurate. We 
  %differentiate a function over the interval [0,1].
  
    N = 100; %Number of grid points
    h = 1/N; %Sets the spacing of the grid points
    x = (0:N-1)'*h; %Defines the grid points

    %Construct the differentiation matrix valid on [0,2pi]
    hD = 2*pi/N; %Defines equispaced points to generate
                 %the differentiation matrix with
    column = [0 .5*(-1).^(1:N-1).*cot((1:N-1)*hD/2)];
    D = toeplitz(column,column([1 N:-1:2])); %First derivative
    
    D = 2*pi*D; %Scales the differentiation matrix so that
                %it is defined for differentiation over [0,1]
    D2 = D^2; %Square the differentiation matrix to calculate
              %the second derivative
  
    comp = D2*(sin(2*pi*x)); %The approximate derivative
    true = -4*pi^2*sin(2*pi*x); %The exact derivative
    
    %plot the exact and approximate derivatives to compare them
    figure()
    plot(x,comp, 'r-');
    hold on
    plot(x,true, 'b:')

    %calculate the maximum error of the approximate solution
    max_error = norm(abs(comp-true),inf)

  %% Interpolation Comparison
  %This script compares interpolation using both equispaced and
  %Chebyshev interpolation points. This script was made using 
  %the MATLAB file P9.m by Trefethen.

  N = 16; %Set the number of grid points
  xx = -1.01:.005:1.01; %makes xx for polynomial to be evaluated at. 
  clf
  for i = 1:2
    if i==1 %carry out the interpolation using equispaced points
        s = 'Equispaced points'; 
        x = -1 + 2*(0:N)/N; %generate equispaced points on [-1,1]
    end
    if i==2 %carry out the interpolation using Chebyshev points
        s = 'Chebyshev points';  
        x = cos(pi*(0:N)/N); %generate Chebyshev points on [-1,1]
    end
    
    %Evaluate function at interpolation points
    u = 1./(1+0.5*sin(pi*x).^6);
    %'Exact' value of function to be used for calculating errors
    uu = 1./(1+0.5*sin(pi*xx).^6); 
    
    p = polyfit(x,u,N);  %carry out interpolation using polyfit
    pp = polyval(p,xx);  %evaluate the interpolating polynomial
                         %using polyval at the same points that
                         %the 'exact' function is evaluated at
    
    subplot(1,2,i) 
    %plot the interpolation points and the exact values of the
    %function at the interpolation points
    plot(x,u,'.','markersize',24)
    line(xx,pp)           %Plots the interpolating polynomial
    axis([-1.1 1.1 0 2]), title(s)
    
    %calculate the maximum error
    error = norm(uu-pp,inf);
    text(-.6,.2,['max error = ' num2str(error)])
    xlabel('$x$', 'Interpreter', 'latex'); 
    ylabel('$u(x)$', 'Interpreter', 'latex');
  end

  %% Integration
  %This script tests the accuracy of integration carried out using
  %Clenshaw-Curtis quadrature
   [x,w] = clencurt(50); %generate Chebyshev points and integration
                         %weights using 50 points. Both are
                         %valid for calculations on [-1,1]
    x = map(0,1,x);      %map the Chebyshev points to [0,1]
    w = w/2;             %map the integration weights to [0,1]
    int_S = w*sin(2*pi*x); %integrate the function sin(2pi*x)
    int_C = w*cos(2*pi*x); %integrate the function cos(2pi*x)
   
function [x,w] = clencurt(N)
  
  %Generate the Chebyshev points
  theta = pi*(0:N)'/N; 
  x = cos(theta);
  
  %Generate the integration weights
  w = zeros(1,N+1); ii = 2:N; v = ones(N-1,1);
  if mod(N,2)==0 
    w(1) = 1/(N^2-1); w(N+1) = w(1);
    for k=1:N/2-1, v = v - 2*cos(2*k*theta(ii))/(4*k^2-1); end
    v = v - cos(N*theta(ii))/(N^2-1);
  else
    w(1) = 1/N^2; w(N+1) = w(1);
    for k=1:(N-1)/2, v = v - 2*cos(2*k*theta(ii))/(4*k^2-1); end
  end
  w(ii) = 2*v/N;
end
   function [x2] = map(a,b,x)
        x2 = 0.5*(b-a)*x + 0.5*(a+b);
      end
  