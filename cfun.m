function cfun=cfun(U)
%
%
%
global DMfsc;
global Nfsc;
global Lfsc;
global BetaH;
global F;


BC=zeros(1,Nfsc-1);
BC(Nfsc-1)=1.;
rhs=zeros(size(U));
rhs(Nfsc-1)=1;

D1=DMfsc(2:Nfsc-1,2:Nfsc,1).*(-2/Lfsc);
D2=DMfsc(2:Nfsc-1,2:Nfsc,2).*(4/(Lfsc^2));




A=[D2+diag(F(1:Nfsc-2))*D1;BC];


cfun=A*U-rhs;
