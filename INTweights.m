function IW=INTweights(N,boxsize)
%
% Computes the integration weigths for a function in the domain of length
% boxsize, interpolated to the [1,-1] interval.
%
% INPUTS:
% N  : the number of points (without the boundaries)
% boxsize  : the size of the domain
%
% OUTPUT:
% IW : Diagonal matrix of the integration weigths.
%
%

nW=0:1:N-1;
jW=0:1:N-1;

bW=ones(1,N);
bW(1)=0.5;
bW(N)=0.5;
cW=2*bW;
bW=bW/(N-1);
S=cos(nW(3:N)'*jW*(pi/(N-1)));
IW=boxsize/2*diag(bW.*[(2+(cW(3:N).*((1+(-1).^nW(3:N))./(1-nW(3:N).^2)))*S)]);

%IW=IW(2:end-1,2:end-1);


