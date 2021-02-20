A=[-5 4 4; 0 -2-2*i 4; 0 0 -.3+i];

e=eig(A);
plot(e,'*'),hold on

theta=0:pi/50:2*pi;
for j=1:length(theta)
  M=exp(i*theta(j))*A;
  Mt=.5*(M+M');
  [v,d]=iord2(Mt);
  x=v(:,1);
  z(j)=x'*A*x/(x'*x);
end
plot(z,'-r')