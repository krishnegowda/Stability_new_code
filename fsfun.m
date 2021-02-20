function [fsfun]=fsfun(U)
%jerome

global DMfsc;
global Nfsc;
global Lfsc;
global BetaH;
global F;
global D1fsc
global D2fsc
global D3fsc


BC1=-2/Lfsc.*DMfsc(1,2:Nfsc,1)*1000;
BC2=-2/Lfsc.*DMfsc(Nfsc,2:Nfsc,1)*1000;


rhs=zeros(size(U));
rhs(Nfsc-1)=1000;
I=eye(Nfsc);
DU=D1fsc*U;



A=[BC1;D3fsc+diag(U(2:Nfsc-2))*D2fsc+BetaH*(-diag(DU)*D1fsc);BC2];
B=[0;BetaH*ones(Nfsc-3,1);0];

fsfun=A*U+B-rhs;
