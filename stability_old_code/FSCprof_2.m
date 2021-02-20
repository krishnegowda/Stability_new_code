function [yfsc,fprim,g]=FSCprof(Nprof,Lprof,m)
%  solve the Falkner-Skan-Cooke Boundary layer profiles
%			
%  Use chebyshev polynomials and explicit enforcement of boundary
%  conditions!
% 
%  The profiles in [0 Lprof] are then rescaled to [1,-1].
%
%  We use fsfun.m and cfun.m, two functions to specify the computation 
%  of the profiles.
%
% inputs:
% Nprof: number of grid nodes
% Lprof: boxsize
% m: pressure gradient parameter
%
% outputs: 
% yfsc: grid nodes
% fprim: u velocity profile
% g: w velocity profile
%
% code written by Markus Hogberg


% variable for fsfun och cfun
global DMfsc;                
global Nfsc;
global Lfsc;
global BetaH;
global F;
global D1fsc
global D2fsc
global D3fsc


Nfsc=Nprof;
Lfsc=Lprof;


% Specify problem parameters
BetaH=2*m/(m+1);


% other way to give the pressure gradiant          
%wedgeangle = 0;           
%BetaH=wedgeangle/180.;

% differenciation matrices and scaling
[x,DMfsc]=chebdif(Nfsc,3);            

%diferentiation matrices for fsfun
D1fsc=DMfsc(3:Nfsc-1,2:Nfsc,1).*(-2/Lfsc);      
D2fsc=DMfsc(3:Nfsc-1,2:Nfsc,2).*(4/(Lfsc^2));
D3fsc=DMfsc(3:Nfsc-1,2:Nfsc,3).*(-8/(Lfsc^3));

%  Initial guess can be very important
f0=-(2*x(2:Nfsc)-2)*Lfsc/4-1.21.*(1-exp((x(2:Nfsc)-1)*Lfsc/2));
c0=-(2*x(2:Nfsc)-2)*Lfsc/4-1.21.*(1-exp((x(2:Nfsc)-1)*Lfsc/2));

%%%%  Compute Falkner-Skan Profile
f=fsolve('fsfun',f0,[0,1e-14,1e-14]);
F=f;
f=[0;f];

% compute derivatives of f
fprim=-2/Lfsc*(DMfsc(:,:,1)*f);    
%fbis=-2/Lfsc*(DMfsc(:,:,1)*fprim);
%fbis0=-2/Lfsc*(DMfsc(1,:,1)*fprim);
y=(1-x(1:Nfsc)).*(Lfsc/2);

%%%% Compute Cooke profile
g=fsolve('cfun',c0,[0,1e-10,1e-10]);
g=[0;g];

%%%% Now normalize to give boundary layer thickness = 1
% compute integration weights
IW=INTweights(Nfsc,Lfsc);

% Now rescale to deltastar scale
deltastar=sum(IW*(ones(Nfsc,1)-fprim));
yfsc=y/deltastar;
Lfsc=Lfsc/deltastar;





