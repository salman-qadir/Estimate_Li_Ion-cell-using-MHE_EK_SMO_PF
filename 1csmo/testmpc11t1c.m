%clearvars; close all; clc;
addpath('C:\Users\salma\MATLAB Drive\T11 MHE\1csmo');
p=param11t1c;
initial11t1c;
load('C:\Users\salma\MATLAB Drive\T11 MHE\onecdata')
g1=lpf(g1,7e-1);g2=lpf(g2,7e-1);%g1=[1,g1];g2=[g2(1),g2];
g3=1:g1(end); g4=interp1(g1,g2,g3);
q1=0.2e-3;qa=q1*randn(1,1e5);r1=0.2e-3;
ra=r1*randn(1,1e5);
save('C:\Users\salma\MATLAB Drive\T11 MHE\1aekfse\noise',"qa","ra","q1","r1")
 %%
p.gt=g3;p.gv=g4+ra(1:length(g4));
% tic
 %%
[dsfhm11t1csscec,dsfhm11t1cd2c]= nmpcss11t1c(p);     
% dsfhm11t1csscec.et=toc;
csm=mean(dsfhm11t1csscec.est(1,3:end))
save('C:\Users\salma\MATLAB Drive\T11 MHE\1csmo\dsfhm11t1csscec','dsfhm11t1csscec')
save('C:\Users\salma\MATLAB Drive\T11 MHE\1csmo\dsfhm11t1cd2c','dsfhm11t1cd2c')
figures11t1c
% rmse =(sqrt((sum((dsfhm11t1bd2.v(1,p.N:end)-dsfhm11t1bssce.y).^2))/length(dsfhm11t1bd2.v)))/mean(dsfhm11t1bd2.v);
%%
%end
%figures10t1a
%  format long
%  dsfhm10t1ss2c.soh(end)
%  dsfhm10t1ss2c.x(end,6)
%dsfhm4cc.rs=rs; dsfhm4cc.soh=soh;rs= [dsfhm4cv.rs,rs];soh=[dsfhm4cv.soh,soh];
%save('C:\Users\salma\MATLAB Drive\T10_SEI\sfhmmpc\dsfhm4cc','dsfhm4cc')
