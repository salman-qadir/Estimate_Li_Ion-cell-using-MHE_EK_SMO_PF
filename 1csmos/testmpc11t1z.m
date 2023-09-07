%clearvars; close all; clc;
addpath('C:\Users\salma\MATLAB Drive\T11 MHE\1csmos');
p=param11t1z;
initial11t1z;
load('C:\Users\salma\MATLAB Drive\T11 MHE\onecdata')
g1=lpf(g1,7e-1);g2=lpf(g2,7e-1);%g1=[1,g1];g2=[g2(1),g2];
g3=1:g1(end); g4=interp1(g1,g2,g3);
q1=2e-3;qa=q1*randn(1,1e5);r1=1e-3;
ra=r1*randn(1,1e5);
save('C:\Users\salma\MATLAB Drive\T11 MHE\1aekfse\noise',"qa","ra","q1","r1")
 %%
p.gt=g3;p.gv=g4+ra(1:length(g4));
tic
 %%
[dsfhm11t1zsscec,dsfhm11t1zd2c]= nmpcss11t1z(p);     
dsfhm11t1zsscec.et=toc;
save('C:\Users\salma\MATLAB Drive\T11 MHE\1csmos\dsfhm11t1zsscec','dsfhm11t1zsscec')
save('C:\Users\salma\MATLAB Drive\T11 MHE\1csmos\dsfhm11t1zd2c','dsfhm11t1zd2c')
figures11t1z
% rmse =(sqrt((sum((dsfhm11t1bd2.v(1,p.N:end)-dsfhm11t1bssce.y).^2))/length(dsfhm11t1bd2.v)))/mean(dsfhm11t1bd2.v);
%%
%end
%figures10t1a
%  format long
%  dsfhm10t1ss2c.soh(end)
%  dsfhm10t1ss2c.x(end,6)
%dsfhm4cc.rs=rs; dsfhm4cc.soh=soh;rs= [dsfhm4cv.rs,rs];soh=[dsfhm4cv.soh,soh];
%save('C:\Users\salma\MATLAB Drive\T10_SEI\sfhmmpc\dsfhm4cc','dsfhm4cc')
