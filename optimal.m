function [flowin,flowot,gg]=optimal(A,T,M,ak2,iflag,imaglow);
%
% This function computes the initial flow structure which 
% achieves the maximum transient energy growth
%
% INPUT
% A       = 3D Orr-Sommerfeld operator 
% T       = time 
% M       = energy weight matrix
% ak2     = alpha^2+beta^2
% iflag   = flag
%           iflag = 1:  compute the maximum growth and 
%                       initial condition in time 
%                       interval [0,T]
%           iflag = 2:  compute the initial disturbance 
%                       yielding maximum growth at time T
% OUTPUT 
% flowin  =  optimal disturbance
%           flowin(1:Nos)         = normal velocity 
%                                   
%           flowin(Nos+1:Nos+Nsq) = normal vorticity 
%                                   
% flowot  = field at optimal time
%           flowot(1:Nos)         = normal velocity 
%                                  
%           flowot(Nos+1:Nos+Nsq) = normal vorticity 
%                                  
 
    global qb;

% Phase 1: Compute eigenvalues and eigenfunctions of 
% Orr-Sommerfeld matrix and sort in order of descending 
% imaginary part. The function nlize normalizes the 
% eigenfunctions with respect to the weight matrix M.

    [xs,es]=iord2(A);
    xs=normalize(xs,M);

% Phase 2: Choose which modes are to be used in optimal 
% calculation. Modes with imaginary part > 1 are neglected. 
% Modes with imaginary part < imin are neglected as well.

    ishift=1; 
    imin=imaglow;

    if imin==-inf
      cols=1:length(es);
    else
      while imag(es(ishift))>1, ishift=ishift+1; end;
      [n1,n2]=n9(es,imin);
      cols=(ishift:n2);
    end
    xu=xs(:,cols);
    eu=es(cols); 
    ncols=length(cols);
    fprintf('Number of modes used: %1.0f \n',ncols); 
    
    % Phase 3: Compute the reduced Orr-Sommerfeld operator
    
    [qb,invF]=QBmat(M,xu,eu);

    % Phase 4: Compute the time for the maximum growth using 
    % the built-in Matlab routine 'fmin'
    
    ts = 0;
    tf = T;
    
    if iflag==1
      % this should be done with numerical range instead:
      % we want to know if there can be initial growth or if everything decays
      gcheck=maxer(1/100); 
      gcheck=gcheck^2;
      if gcheck<1,
        tformax=0;
        mgrowth=1;
      else
	
       %%% options for the optimal control using fminbnd,
       %%% modified by 2004.08.10 Espen in order to match
       %%% the new version of builtin fminbnd.
       options.disp='off';
       options.tolx=1.e-3;
       tformax=fmy(ts,tf,options);
       mgrowth=maxer(tformax);
       mgrowth=mgrowth^2;
      end;
      fprintf('Time for maximum growth:  %e \n',tformax);
    else
      tformax=T;
      ts=0;
      %%% tf was not implemented originally for this option.
      %%% this definition of tf is coincidential, espen
      tf=tformax+tformax/10;
    end; 
    
    % Phase 5: Compute the initial condition that yields the 
    % maximum growth. This is obtained by
    % (1) computing the matrix exponential evaluated at the 
    %     optimal time;
    % (2) computing the SVD of the matrix exponential 
    %     exp(-i*A*t)=USV. 
    % The initial condition that yields the maximum growth is 
    % the first column of V. To convert the initial condition 
    % to a vector of coefficients in the eigenfunction basis 
    % multiply by the matrix of eigenfunctions and inv(F) 
    
    evol=expm(tformax*qb);
    [U,S,V]=svd(evol);
    mgrowth=S(1,1)^2;
    fprintf('Maximum growth in energy:  %e \n',mgrowth);
    
    
    %%%%%% states at time ts and tformax
    flowin=sqrt(2*ak2)*xu*invF*V(:,1);
    flowot=sqrt(2*ak2)*xu*invF*U(:,1);
    flowin_modes = V(:,1);
    
    ntimes = 80;
    gg = zeros(ntimes, 5);
    for i=1:ntimes,
      tid = ts + (tf-ts)*(i-1)/(ntimes-1);
      gg(i,1) = tid;
      % max amplification
      evol = expm(tid*qb);
      gg(i,2) = norm(evol)^2;
      % amplification with the optimal perturbation
      flow = sqrt(2*ak2)*xu*invF*evol*flowin_modes;
      gg(i,3) = (flow'*M*flow) / (flowin'*M*flowin);
      % part of the energy in normal velocity vs normal vorticity
      [m, ~] = size(xu);
      mask1 = ones(m, 1);
      mask2 = ones(m, 1);
      mask1(m/2+1:end) = 0;
      mask2(1:m/2) = 0;
      flow1 = mask1 .* flow;
      flow2 = mask2 .* flow;
      gg(i,4) = (flow1'*M*flow1) / (flowin'*M*flowin);
      gg(i,5) = (flow2'*M*flow2) / (flowin'*M*flowin);
    end 
    









