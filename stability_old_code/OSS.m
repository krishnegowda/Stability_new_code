function OSS=OSSmat(nosmod,D1,D2,d4,alp,beta,R,u,up,upp)

% Compute the Orr-Somerfeld-Squire operator for the homogeneous problem.

k2=alp^2+beta^2;

% implement homogeneous boundary conditions 
D1=D1(2:nosmod+1,2:nosmod+1);
D2=D2(2:nosmod+1,2:nosmod+1);
u=u(2:nosmod+1);
up=up(2:nosmod+1);
upp=upp(2:nosmod+1);

% fourth derivative with clamped conditions
[x,D4]=cheb4c(nosmod+2);
D4=D4*d4;

% unity matrix
I=eye(size(D4));
Z=zeros(size(I));

%%%% laplacian
delta=(D2-k2*I);
delta2=(D4-2*k2*D2+k2*k2*I);


%%%% compute the matrices LOS, LC, LSQ
LOS  =                i*( alp*diag(u)  ) * delta  ...
                           -i*(  alp*diag(upp)  ) ...
                           - delta2 / R   ;
     
LC    = i*(     - beta*diag(up)    );

LSQ   =-i*( alp*diag(u) ) + delta/R;

OSS =    [-delta\LOS    Z         ; 
               LC                       LSQ              ];

OSS = i*OSS;

%in main ylen=1;and y=(y+1)/2 when computing base flow
%x=(x+1)/2;
%k2=alp^2*ones(size(x))+n2./x.^2;
%delta=diag(k2.*x.^2)
%T=diag(1./x.^2)-diag(1./x)*D1*diag(1./k2./x)*D1;
%S=diag(k2.^2.*x.^2)-diag(1./x)*D1*diag(k2.*x.^3)*D1;