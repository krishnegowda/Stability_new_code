function output = Neutral_alpha_beta()

% compute the Orr-Sommerfeld matrix for three-
% dimensional Poiseuille,
% compute the energy weight matrix and determines
% maximum transient growth in the
% alpha-beta-plane
% Plot the maximum optimal growth (fig 1)
% and the least stable eigenvalue (fig 2)
% in the alpha-beta plane
%
% INPUT
%
% Re        = Reynolds number
% N         = total number of modes for normal velocity
% T         = compute maximum growth in time interval [0 T]


clear

global D0 D1 D2 D4
global qb

zi = sqrt(-1);

%...input data
N      = input('Enter the number of Chebyshev polynomials: ');
Re     = input('Enter Reynolds number: ');
Tmax   = input('Enter Tmax: ');
T      = [0 Tmax];

%...generate Chebyshev differentiation matrices
[D0,D1,D2,D4] = ChebMat(N);

nreso       = 100;
alpha_min   = 0;   alpha_max = 2;
beta_min    = 0.1; beta_max  = 4;
beta_range  = linspace(beta_min,beta_max,nreso);
alpha_range = linspace(alpha_min,alpha_max,nreso);

for i=1:nreso
  for j=1:nreso
    beta  = beta_range(i);
    alpha = alpha_range(j);
    [A,B] = CouetteMatrix(N,alpha,beta,Re);%PoiseuilleMatrix(N,alpha,beta,Re);
    
    %...generate energy weight matrix
    k2 = alpha^2 + beta^2;
    M  = EnergyMatrix(N+1,N+1,k2);
    
    %...compute the Orr-Sommerfeld matrix (by inverting B)
    OS = inv(B)*A;
    eOS = eig(OS);
    
    %...compute the optimal
    [flowin,flowot,gg] = Optimal(OS,T,M,k2,1);
    Gmax(i,j) = max(gg(:,2));
    emax(i,j) = max(imag(eOS));
  end
end

%...graphics
figure(1)
title('Optimal growth')
contour(alpha_range,beta_range,log10(Gmax),30);colorbar
figure(2)
title('Neutral curve')
contour(alpha_range,beta_range,emax,30);colorbar
hold on
contour(alpha_range,beta_range,emax,[0 0],'k','Linewidth',3)

output = {};

end