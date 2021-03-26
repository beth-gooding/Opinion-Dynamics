function C = ConvolutionMatrixQ(x,R)
%Constructs a matrix which carries out both the interpolation from 
%equispaced to Chebyshev points and numerical integration needed
%to calculate the convolution integral in the calculation of the 
%order parameter Qc.

Ncheb = 50; %Number of Chebyshev points

%Generates Chebyshev points on [$-$1,1]
[~,xcFull] = cheb(Ncheb); 

%Generates integration weights on [$-$1,1]
wFull = integrationWeights();

%% Convolution Integral
C = zeros(length(x), length(x)); %Initialises the convolution matrix

for i = 1:length(x) %For each grid point
    
    xi = x(i); %One point at a time
    
    % Find the interval [xi $-$ R, xi + R]
    [a,b] = interval(xi);
    
    % Construct Chebyshev points in [xi $-$ R, xi + R]
    xc = map(xcFull);
 
    %Construct interpolation matrix for [xi $-$ R, xi + R]
    Interp = fftInterpMatrix(x, xc);
    
    %Map the integration weights to [xi $-$ R, xi + R]
    w = wFull*(b-a)/2;
    
    %Calculate ith row of convolution matrix
    C(i,:) = (w.*ones(length(xc),1)')*Interp; 
end

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