function PDE_Integral_Matrix()
%This script tests whether the convolution integral can be 
%evaluated accurately using the convolution matrix.
set(0, 'DefaultAxesFontSize', 30);
set(0, 'DefaultLineLineWidth', 3);
close all
%% Parameters and Initial Things
R = 0.2;  %Half the width of the confidence interval
N = 50; %Number of equispaced grid points
h = 1/N; %Grid point spacing
x = (0:N-1)'*h; %Equispaced grid points on [0,1]

%Test functions evaluated at the grid points
rho_s = sin(2*pi*x);
rho_c = cos(2*pi*x);

Ncheb = 16; %Number of Chebyshev points
[~,xcFull] = cheb(Ncheb); %Generate Chebyshev points on
                          %[$-$1,1]
wFull = integrationWeights(); %Generate integration weights
                              %on [$-$1,1]


%% Convolution Integral
C = zeros(length(x), length(x)); %Set up convolution matrix
for i = 1:length(x) %for each grid point
    
    xi = x(i); %One grid point at a time
    
    % Find the interval [xi $-$ R, xi + R]
    [a,b] = interval(xi);
    
    % Construct Chebyshev points in [xi $-$ R, xi + R]
    xc = map(xcFull);
 
    %Construct interpolation matrix for [xi $-$ R, xi + R]
    Interp = fftInterpMatrix(x, xc);
    %map integration weights to [xi $-$ R, xi + R]
    w = wFull*(b-a)/2; 
    %Calculate ith row of convolution matrix
    C(i,:) = (w.*(xi-xc)')*Interp; 
end

int_s = C*rho_s; %Approximate solution
exact_s = exactIntegral_s(xc); %Exact solution
error_s = norm(abs(int_s-exact_s), inf); %max error

figure()
subplot(1,2,1) %Plot the results relating to sin(2pix)
plot(x,int_s,'r')
hold on
plot(x,exact_s,'k--')
xlabel('$x$', 'Interpreter', 'latex'); 
ylabel('$I_s$', 'Interpreter', 'latex');
legend('$\hat{f}(x)$', '$f(x)$', 'Interpreter', 'latex')

int_c = C*rho_c; %Approximate solution
exact_c = exactIntegral_c(xc); %Exact solution
error_c = norm(abs(int_c-exact_c), inf); %max error

subplot(1,2,2) %Plot the results relating to cos(2pix)
plot(x, int_c, 'r')
hold on
plot(x, exact_c, 'k--')
xlabel('$x$', 'Interpreter', 'latex'); 
ylabel('$I_c$', 'Interpreter', 'latex');
legend('$\hat{f}(x)$', '$f(x)$', 'Interpreter', 'latex')

save('interpolation_check_data')

%% Functions
     function [D,x] = cheb(N) %Cheb by Trefethen
      if N==0, D=0; x=1; return, end
        x = cos(pi*(0:N)/N)'; 
        c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
        X = repmat(x,1,N+1);
        dX = X-X';                  
        D  = (c*(1./c)')./(dX+(eye(N+1)));      % off-diagonal entries
        D  = D - diag(sum(D'));    % diagonal entries
      end

      function [x2] = map(x)
        x2 = 0.5*(b-a)*x + 0.5*(a+b);
      end

      function [xm,xp] = interval(x)
          xm = x - R;
          xp = x + R;

    %     X = x(j) - x;
    %     mask = abs(X) > 0.5;
    %     X(mask) = X(mask)-sign(X(mask));
    %     x_ivl = x(abs(X)<=R);
      end

      function I = exactIntegral_s(x)
            yp = x+R;
            ym = x-R;
            alpha = 2*pi;
            I = -( sin(alpha*yp)/alpha^2 - R.*cos(alpha*yp)/alpha ) ...
                + ( sin(alpha*ym)/alpha^2 + R.*cos(alpha*ym)/alpha ) ;
      end
  
    function I = exactIntegral_c(x)
           yp = x+R;
           ym = x-R;
           alpha = 2*pi;
           I = -R*sin(alpha*yp)/alpha -R*sin(alpha*ym)/alpha ...
               -cos(alpha*yp)/alpha^2 + cos(alpha*ym)/alpha^2;
        
    end

    function  Inter  = fftInterpMatrix( x, interp )
        %%FFT Interpolation
        LL        = 2*pi;
        phi       = interp*LL;
        N         = length(x);
        M     = floor((N+1)/2);

        n1   = 1:(M-1);       % first half of points, excluding zero
        n2   = (M+1):(N-1);   % second half

        interpol  = [exp(2*pi*1i*phi/LL*0) ,...
            exp(2*pi*1i*phi/LL*n1) , ...
            cos(2*pi*M*phi/LL) ,...
            exp(2*pi*1i*phi/LL*(n2-N))]/N;

        n = (0:N-1);
        k = (0:N-1);

        FFTMatrix = exp( -1i * 2*pi *n'*k /N);
        Inter = real(interpol*FFTMatrix);
    end

  function w = integrationWeights()
        theta = pi*(0:Ncheb)'/Ncheb;
        w    = zeros(1,Ncheb+1);    % Initialise w, the weights
        ii   = 2:Ncheb;             
        v    = ones(Ncheb-1,1);     

        % For even numbers
        if mod(Ncheb,2)==0          % Check for even total cheb points
            w(1)    = 1/(Ncheb^2-1);   
            w(Ncheb+1)  = w(1);     
            for k=1:Ncheb/2-1       
                v = v - 2*cos(2*k*theta(ii))/(4*k^2 - 1);   
            end
            v = v - cos(Ncheb*theta(ii))/(Ncheb^2-1);               
        else
            w(1)   = 1/Ncheb^2;      
            w(Ncheb+1) = w(1);     
            for k=1:(Ncheb-1)/2    
                v = v - 2*cos(2*k*theta(ii)) / (4*k^2 - 1);
            end
        end

        w(ii) = 2*v/Ncheb;
    end
        
end