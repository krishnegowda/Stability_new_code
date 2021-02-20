clc
clear all
iflag=input('Poiseuille (1) or Couette flow (2) ');

y1=-1:0.1:1;

if iflag==1
    U=1-y1.^2;
  
else
    U=y1;
end


plot(U,y1)

