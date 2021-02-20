function [xs,es]=iord2(d);
%
% This function computes the eigenvalues of a matrix d and 
% orders the eigenvalues so that the imaginary parts are 
% decreasing.
%
% INPUT 
% d    = input matrix
%
% OUTPUT
% es   = ordered eigenvalues
% xs   = eigenvectors

    [v,e]=eig(d);
    e=diag(e);
    [eimag,is]=sort(-imag(e));
    xs=v(:,is); 
    es=e(is);


