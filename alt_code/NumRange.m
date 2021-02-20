
function NumRange()

%  compute the spectrum and numerical range
%
% INPUT
%
% Re        = Reynolds number
% alpha     = alpha (streamwise wave number)
% beta      = beta  (spanwise wave number)
% iflow     = type of flow  (Poiseuille=1, Couette=2)
% N         = total number of modes for normal velocity

clear

global D0 D1 D2 D4

zi=sqrt(-1);
%...input data
% iflow  = input('Poiseuille (1) or Couette flow (2) ');
% N      = input('Enter the number of Chebyshev polynomials: ');
% Re     = input('Enter the Reynolds number: ');
% alpha  = input('Enter alpha: ');
% beta   = input('Enter beta: ');

iflow  = 2 % input('Poiseuille (1) or Couette flow (2) ');
N      = 100 %input('Enter the number of Chebyshev polynomials: ');
Re     = 1000 %input('Enter the Reynolds number: ');
alpha  = 0 %input('Enter alpha: ');
beta   = 2%input('Enter beta: ');

%...generate Chebyshev differentiation matrices
[D0,D1,D2,D4] = ChebMat(N);

%...set up Orr-Sommerfeld matrices A and B
if (iflow == 1)
  [A,B] = PoiseuilleMatrix(N,alpha,beta,Re);
else
  [A,B] = CouetteMatrix(N,alpha,beta,Re);
end

%...generate energy weight matrix
k2 = alpha^2 + beta^2;
M  = EnergyMatrix(N+1,N+1,k2);

%...compute the Orr-Sommerfeld matrix (by inverting B)
OS = inv(B)*A;

%...compute the numerical range
EM    = GetMatrix(OS,M,k2);
e     = eig(EM);
theta = linspace(0,2*pi,100);
for i=1:100
  th          = theta(i);
  M           = exp(sqrt(-1)*th)*EM;
  Mhat        = (M + M')/2;
  [X,D]       = eig(Mhat);
  dd          = diag(D);
  [dmax,imax] = max(dd);
  xmax        = X(:,imax);
  numran(i)   = (xmax'*EM*xmax)/(xmax'*xmax);
end

%...graphics
figure(1)
fill(real(numran),imag(numran),[0.7 0.7 0.7])
hold on;
plot(real(numran),imag(numran),'b')
plot([-1 1.5],[0 0],'r','LineWidth',2)
plot(real(e),imag(e),'k.','MarkerSize',18);
xlabel('$\omega_{r}$','Interpreter','latex')
 ylabel('$\omega_{i}$','Interpreter','latex')
 ax=gca;
ax.FontSize=16;
ax.LabelFontSizeMultiplier=1.2;
box on
%grid on
hold off

output = {};

end