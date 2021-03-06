 %%%%
% calculates transient growth for poiselle or blasius flow
%
%

% nosmod    = number of Orr-Sommerfeld modes
% R         = Reynolds number
% alpha     = alpha (streamwise wave number)
% beta      = beta  (spanwise wave number)
% iflow     = type of flow  (Poiseuille=1, Couette=2, Blasius=3)
% nosmod    = total number of modes for normal velocity
% iflag     = 
%             iflag = 1: compute the maximum growth and 
%                        initial condition in time 
%                        interval [0,T]
%             iflag = 2: compute the initial disturbance 
%                        yielding maximum growth at time T
%
%

close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% setting input parameters
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%params_poiseuille;
%params_blasius;
params_blasius_liftup;

k2=alpha^2+beta^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% generate differentiation matrices using matlab 
%%% differentiation suite
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[y,DM] = chebdif(nosmod+2,2);

%%% scale the differentiation matrices if area is [0, ylen]
if iflow == 3
  yphys = 0.5*ylen*(y+1);
else
  ylen = 2;
  yphys = y;
end
d1=2/ylen; d2=d1*d1; d3=d1*d2; d4=d1*d3;
D1=DM(:,:,1)*d1;    
D2=DM(:,:,2)*d2;    

%ylen=8;
%y_t=y*ylen/2
%u=exp(yt^2/sigma)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% generate profiles
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (iflow==1)
  u   = 1-y.^2;
  up  = -2*y;
  upp = -2*ones(size(y));
elseif (iflow==2)
  u = y;
  up = ones(size(y));
  upp = zeros(size(y));
else  
  %%% computing the similarity solution
  [yfsc,fprim,g]=FSCprof(Nprof,Lprof,m);     
  Lb=max(yfsc);
  if(Lb < ylen)
    error(['... box to large, max ' ,num2str(Lb),' or recompute profiles.'])
  end
  Nyfsc=length(fprim);

  %%% Interpolation from [0, ylen] to [-1,1]
  fprim=fprim(Nyfsc:-1:1);
  yto =-(1-(ylen/2*(1-cos((nosmod+1:-1:0)*pi/(nosmod+1))))*2/Lb);
  u=chebint(fprim,yto);
  %%% Differentiation
  up=D1*u;
  upp=D2*u;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% get the Orr-Sommerfeld-Squire matrix, boundary 
%%% conditions are implicitly set. Now on inner points 
%%% only
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=OSS(nosmod,D1,D2,d4,alpha,beta,R,u,up,upp);


%%% Get the energy matrix for inner points 
IWT=INTweights(nosmod+2,ylen);
M=ENER(D1,IWT,k2,nosmod);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% Compute the Transient Growth 
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imaglow=-0.5; 
%imaglow=-inf;
[flowin,flowout,gg]=optimal(A,T,M,k2,iflag,imaglow);
flowin = normalize(flowin, M);
flowout = normalize(flowout, M);
figure(1),clf; 
plot(gg(:,1),gg(:,2),'k-','Linewidth',1.2);
hold on;
plot(gg(:,1), gg(:,3), 'b-');
hold off;
%axis ([0 100 1 30])
set(gca,'YScale','log')
set(gca,'XScale','lin')
xlabel('t')
ylabel('G(t)')
ax=gca;
set (ax, "FontSize", 14);
%legend([a b],'Experiments','Simulations')
box on
grid on

% plot vorticity and velocity
figure(10);
plot(gg(:,1), gg(:,4), 'b-', 'Linewidth', 1.2);
hold on;
plot(gg(:,1), gg(:,5), 'r-', 'Linewidth', 1.2);
hold off;
grid on;
ylim([min(abs(gg(:, 4))) max(max(abs(gg(:, 4:5))))]);
set(gca,'YScale','log')
set(gca,'XScale','lin')
xlabel('t')
ylabel('G(t)')
ax=gca;
set (ax, "FontSize", 14);

% plot initial and final flows
figure(2);
plotflow(yphys, flowin);
figure(3);
plotflow(yphys, flowout);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% Compute the Optimal Response
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iflag=1;
[flowin,flowout,or]=opt_response(A,omega,M,k2,iflag,imaglow);
figure(6); semilogy(or(:,1),or(:,2));grid on
xlabel('f')
ylabel('R(f)')

% plot forcing and response
figure(7);
plotflow(yphys, flowin);
figure(8);
plotflow(yphys, flowout);
