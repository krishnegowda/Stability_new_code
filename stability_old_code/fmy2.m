
    function t=fmy2(t1,t2,options,M,eu,xu);
% 
%
% This function uses the built-in function 'fminbnd' to find 
% the maximum value of a function on the interval [t1,t2]. 
% The function is in the file maxres.m
% 
% INPUT: 
% t1,t2   = lower and upper bounds of interval
% options = input parameters for minimization routine 
%
% OUTPUT
% t       = value at which function maxer(t) is minimized
%
% NOTE
% the function fmin has been replaced by fminbnd 2004.08.10
% which in turn requires a different setup of options, which is 
% now a struct array, Espen


    f1=maxres(t1,M,xu,eu);
    f2=maxres(t2,M,xu,eu);

    tt=fminbnd(@(om) maxres(om,M,xu,eu),t1,t2,optimset('TolX',options.tolx,'Display',options.disp));

    f3=maxres(tt,M,xu,eu);
    f=[f1 f2 f3]; 
    tm=[t1 t2 tt];
    [y,is]=sort(f);
    t=tm(is(1));
