


     function a=maxres(om,M,xu,eu);
%
% This function computes the norm of the matrix qb

     
     [qb,invF]=Rmat(M,xu,eu,om);
     a=-norm(qb);

