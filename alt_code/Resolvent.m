function output = Resolvent()

%  compute the resolvent norm for real and complex frequency omega
%
% INPUT
%
% Re        = Reynolds number
% alp       = alpha (streamwise wave number)
% beta      = beta  (spanwise wave number)
% iflow     = type of flow  (Poiseuille=1, Couette=2)
% N         = total number of modes for normal velocity

clear

global D0 D1 D2 D4

zi=sqrt(-1);
% input data

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


% generate Chebyshev differentiation matrices
[D0,D1,D2,D4] = ChebMat(N);

% set up Orr-Sommerfeld matrices A and B

if (iflow == 1)
  [A,B] = PoiseuilleMatrix(N,alpha,beta,Re);
else
  [A,B] = CouetteMatrix(N,alpha,beta,Re);
end

% generate energy weight matrix
k2 = alpha^2 + beta^2;
M  = EnergyMatrix(N+1,N+1,k2);

% compute the Orr-Sommerfeld matrix (by inverting B)
OS = inv(B)*A;

% compute the optimal
eee = eig(OS);
%figure(1);
%plot(real(eee),imag(eee),'o')
%axis([-2 2 -2 2]);
%axis image;

[F,e,invF] = GetMatrixParts(OS,M,k2);
nreso = 80;
for i=1:nreso
  for j=1:nreso
    zr = -0.5 + 2*(i-1)/(nreso-1);
    zi = -1 + 1.5*(j-1)/(nreso-1);
    zz = zr + sqrt(-1)*zi;
    dd = diag(1./(e-zz));
    Reso(i,j) = log(norm(F*dd*invF));
  end
end
for i=1:nreso
  zr = -0.5 + 2*(i-1)/(nreso-1);
  zz = zr;
  dd = diag(1./(e-zz));
  Reso_r(i) = (norm(F*dd*invF));
end

figure(1);subplot(1,1,1,'Fontsize',14)
semilogy(linspace(-0.5,1.5,nreso),Reso_r)
title('Resolvent norm')
ylabel('R');xlabel('\omega')
grid on
figure(2);subplot(1,1,1,'Fontsize',14)
contour(linspace(-0.5,1.5,nreso),linspace(-1,0.5,nreso),Reso')
hold on;
plot([-0.5 1.5],[0 0],'r','LineWidth',2)
plot(real(e),imag(e),'k.','MarkerSize',18);

title('Resolvent norm')
xlim([-0.5 0.8])
xlabel('$\omega_{r}$','Interpreter','latex')
 ylabel('$\omega_{i}$','Interpreter','latex')
 ax=gca;
ax.FontSize=16;
ax.LabelFontSizeMultiplier=1.2;
box on

%ylim([-1.4 0.4])
hold off

output = {};

end