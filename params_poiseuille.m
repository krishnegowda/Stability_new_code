%%%  set the flow type:
%%%  iflow=1 is Poiseuille, iflow=2 is Couette
iflow=1;

%%% time parameters for transient growth
iflag=1; % find the time for maximal growth
T = 100;

%iflag=2; T=40;
%iflag=2 to have the initial and final optimal

%%% parameters for frequency response
omega=[-.1 1.2];


%%% blasius specific parameters
ylen=20;
m=0;
Nprof=99; Lprof=39;
%%% general variables
R=5772;
nosmod=120;
alpha=1.02;
beta=0;
