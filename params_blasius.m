%%%  set the flow type:
%%%  iflow=1 is Poiseuille, iflow=2 is Couette
iflow=3;

%%% time parameters for transient growth
iflag=1; % find the time for maximal growth
T = 300;
ntimes_opt = 60;

%iflag=2; T=40;
%iflag=2 to have the initial and final optimal

%%% parameters for frequency response
omega=[-.05 0.35];


%%% blasius specific parameters
ylen=30;
m=0;
Nprof=160; Lprof=70;
%%% general variables
R=550;
nosmod=150;
alpha=0.24;
beta=0.;
