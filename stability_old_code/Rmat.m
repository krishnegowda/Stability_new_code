
    function [qb,invF]=Rmat(M,xu,e,omega)
%
% This function computes the matrix  Q=-i*F*(diag(e)-i*omega*I)*inv(F) 
% which is used to compute the maximum response
% (see Reddy and Henningson, "Energy Growth in Viscous 
% Channel Flows", JFM 252, page 209, 1993). 
% 
% INPUT
% M     =  energy weight matrix
% xu    =  matrix of eigenfunctions (expansion coefficients)
% e     =  vector of eigenvalues of the stability matrix
% omega =  frequency of the forcing 
%
% OUTPUT
% qb    =  output matrix Q
% invF  =  inverse of F  
%

% Phase 1: compute inner product of the eigenfunctions 
%          in energy norm

A=xu'*M*xu;

% Phase 2: compute decomposition A=F^*F

    [U,S,V]=svd(A);
    s=sqrt(diag(S));
    F=diag(s)*U';
    invF=U*diag(ones(size(s))./s);

% Phase 3: compute Q=-i*F*(diag(e)-i*omega*I)*inv(F)

    qb=-sqrt(-1)*F*inv(diag(e)-omega*eye(size(F)))*invF;



















