function [flowin,flowot,gg]=opt_response(A,omega,M,ak2,iflag,imaglow);
%
% This function computes the initial flow structure which 
% achieves the maximum transient energy growth
%
% INPUT
% A       = 3D Orr-Sommerfeld operator 
% omega   = frequency matrix 
% M       = energy weight matrix
% ak2     = alpha^2+beta^2
% iflag   = flag
%           iflag = 1:  compute the maximum response 
%                       in frequency interval [0,omega] 
%                       
%           iflag = 2:  compute the forcing profile
%                       yielding maximum energy at frequency omega(2)
%
% OUTPUT 
% flowin  =  optimal forcing
%           flowin(1:Nos)         = normal velocity 
%                                   
%           flowin(Nos+1:Nos+Nsq) = normal vorticity 
%                                   
% flowot  = field at optimal forcing
%           flowot(1:Nos)         = normal velocity 
%                                  
%           flowot(Nos+1:Nos+Nsq) = normal vorticity 
%                                  

% resolution
neps = 100;
nw_f = 80;
nw_r = 30;
nw_i = 30;
% vertical limits for resolvent
wimin = -1.;
wimax = 0.2;

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

    omegas=omega(1);
    omegaf=omega(2);
    
    % Phase 3: Compute the reduced Orr-Sommerfeld resolvant
    
%     [qb,invF]=Rmat(M,xu,eu,omegas);

    % Phase 4: Compute the time for the maximum growth using 
    % the built-in Matlab routine 'fmin'
    
    if iflag==1,
      %%% options for the optimal control using fminbnd,
      %%% modified by 2004.08.10 Espen in order to match
      %%% the new version of builtin fminbnd.
      options.disp='off';
      options.tolx=1.e-3;
      ommax=fmy2(omegas,omegaf,options,M,eu,xu);
      mgrowth=maxres(ommax,M,xu,eu);
      fprintf('Frequency for maximum growth:  %e \n',ommax);
    else
      ommax=omegaf;
    end; 
    
    % Phase 5: Compute the forcing profile that yields the 
    % maximum growth. This is obtained by
    % (1) computing the matrix resolvant evaluated at the 
    %     optimal frequency;
    % (2) computing the SVD of the resolvant matrix 
    %     (A-i*omega*I)=USV. 
    % The profile that yields the maximum growth is 
    % the first column of V. To convert the initial condition 
    % to a vector of coefficients in the eigenfunction basis 
    % multiply by the matrix of eigenfunctions and inv(F) 
    
    [qb,invF]=Rmat(M,xu,eu,ommax);
    [U,S,V]=svd(qb);
    mgrowth=S(1,1);
    if iflag==1
      fprintf('Maximum growth in energy:  %e \n',mgrowth);
    else
      fprintf('Growth in energy at omega %e:  %e \n',ommax,mgrowth);
    end
    
    %%%%%% states at time ts and tformax
    flowin=sqrt(2*ak2)*xu*invF*V(:,1);
    flowot=sqrt(2*ak2)*xu*invF*U(:,1);


    omegas=omega(1);
    omegaf=omega(2);
    for i=1:nw_f
      om = omegas + (omegaf-omegas)*(i-1)/(nw_f-1);
      [qb,invF]=Rmat(M,xu,eu,om);
      gg(i,2) = norm(qb);
      gg(i,1) = om;
    end 

    

    %PSEUDOSPECTRA: eigenvalues from perturbed matrix
    epsilon=1e-5;
    figure(5);plot(eu,'*'); hold on; axis([omegas omegaf wimin wimax])
    for i=1:neps,
       E=randn(size(A));
       ee=E/norm(E)*epsilon;%renormalise pert. to norm epsilon
       [xs,es]=iord2(ee+A);
       plot(es,'.g')%plot perturbed eigenvalues
       grid on
       axis ([omegas omegaf wimin wimax])
       set(gca,'YScale','lin')
       set(gca,'XScale','lin')
       xlabel('\omega_{r}')
       ylabel('\omega_{i}')
    end

    
    %Resolvent norm: loop over values of c_i and c_r
    %If cr+i*ci close enough to an eigenvalue the program would stop
    figure(4);plot(eu,'*'); hold on; axis([omegas omegaf wimin wimax])
    cr=linspace(omegas,omegaf, nw_r);
    ci=linspace(wimin, wimax, nw_i);;
    for i=1:length(cr),
      disp(['c_r = ', num2str(cr(i))]);
      for j=1:length(ci),
	      [qb,invF]=Rmat(M,xu,eu,cr(i)+sqrt(-1)*ci(j));
	      reno(i,j) = norm(qb);
      end
    end
    level=[1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6];
    contour(cr,ci,log10(reno)')


