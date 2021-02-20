function M=ENER(DM,IWT,k2,N)
% Compute energy measure matrix Q

MvT= 0.125 * (DM(:,:,1)' * IWT * DM(:,:,1) / k2 + IWT) ;
MetaT=0.125 * IWT / k2;
MT=      [MvT   zeros(N+2,N+2) ; zeros(N+2,N+2)   MetaT ];

Mv=MvT(2:N+1,2:N+1);
Meta= MetaT(2:N+1,2:N+1);
M=       [Mv   zeros(N,N) ;  zeros(N,N)   Meta];
