function [d,ds] = nmpcss11t1a(w)
d=w.y0;d.x1=w.y0.x;d.qf=w.qb;
d.x0=d.x1;%repmat(d.x1,1,w.N);
d.t = [];d.t0 =w.tme;u0=w.y0.u0;d.x=w.y0.x;
w.opt=optt(1,d.x0);w.ua=u0(1);
d.y=w.y0.y;ds=d;w.u=[u0,u0];w.u0=u0;
ds.est=0;w.est=0;d.est=0;
for ij=1:w.N-1
    w.ij=ij;
    [ds,w] = plant(ds,w);
    w.tc=w.t1;
     da = fhmss11t1asslyte(d.qf(:,end),w);
     d.cel=[d.cel,da.cel];d.qf=[d.qf,fm(da.qf)'];
end
%% loop
for ij=w.N:3858-1
    %% 
      w.ij=ij;
      [ds,w] = plant(ds,w); 
      w.tc=w.t1;
         tic
      [d,p] = estimate(d,w);
       w.est=toc;
      [w] = control(w,d);
      g = loop(w,d,ds,ij); if g==0;   d.u=w.u;ds.u=w.u; break;  end
end
end
function [ds,w]=plant(ds,w)
      [d,w] = fhmss11t1assa(ds,w);
      ds=store(ds,d,w);
end
%%
function [d,p] = estimate(d,w)
   %  w.tc=w.t2/2;w.xm=w.xm2;w.xxm=w.xxm2;w.t=w.t2;w.mk=w.mk2;
 %    A   = [];b   = [];
%     Aeq = [];beq = []; % saving flag and  out in struct is good.
 %    [x,V,flag,out]=fmincon(@(x) costfe(w,d,x),d.x0,A,...
%     b,Aeq,beq,w.lb,w.ub,@(x) constre(w,d,x), w.opt);
 %    d.x1=x((w.N-1)*w.no+1:w.N*w.no);d.x0=repmat(d.x1,1,w.N);
  %   y = fhmss11t1assb1(d,w);
    [d,d1,p] = kalf(d,w);%d.x0=repmat(d.x1,1,w.N);
  % d1=fhmss11t1assd(d.x1,d.qf(:,end),w);    
    d = store(d,d1,w);        
end
function[d,d1,p] = kalf(d,p)
d.xf=d.x0;
yt= fhmss11t1assb1(d,p);
xp=fm(yt.x);
pp=p.smd.A*p.pq*p.smd.A'+p.q1*eye(p.bn);
p.smd.C=[p.hx(yt.x(1)/p.csp,yt.xnn(1)),0,0,0];
% predicted error covariance
yp=yt.v(1);
% predicted output
ye=p.y(:,p.ij)-yp;% predicted output error
%ye=p.gv(:,p.ij)-yp;% predicted output error
s=p.smd.C*pp*p.smd.C'+p.r1*eye(p.dn);
%correction term
kg=(pp*p.smd.C')/s; %kalman gain
da.xp=fm(xp)'+kg*ye; %updated state
p.pq=pp-kg*p.smd.C*pp'; %updated covariance
d.xf=da.xp;
d1= fhmss11t1assb1(d,p);  %updated output
dd.yp=d1.v(1);
d.ye(:,p.ij)=p.y(:,p.ij)-dd.yp; %updated output error
%d.ye(:,p.ij)=p.gv(:,p.ij)-dd.yp; %updated output error
%p.xe0=p.xp(:,p.ij);
d1.x=fm(da.xp);d.x1=d1.x;
d.x0=d.x1;
 end
function [w] = control(w, d)
   w.tc=w.t2/2;w.xm=w.xm2;w.xxm=w.xxm2;w.t=w.t2;w.mk=w.mk2;
%     A   = [];b   = [];
%     Aeq = [];beq = [];p=0;
%     lb  = -(((abs(w.cr)+p)*w.c)/w.a)*ones(1,w.M);
%     ub  = (((abs(w.cr))*w.c)/w.a)*ones(1,w.M)*0;
%     [u,~,~,~]=fmincon(@(u) costfc(w,d,u),w.u0,A,...
%     b,Aeq,beq,lb,ub,@(u) constrc(w,d,u), w.opt);
%     w.u0= [u(:,2:size(u,2)) u(:,size(u,2))];
%     w.u(w.ij+1)=u(1);w.ua=u(1);
   w.u(w.ij+1)=w.u0(1);w.ua=w.u0(1);
end
function cost = costfc(w,da,u)
      y = fhmss11t1assc(u,da,w);
 % c
%     q=1e9; k=1.5e2;j=0*1e-1; 
%                costr=sum( k*(fm(y.socn)-100).^2 +...
%                q*((fm(y.x(:,4)))) + j*u.^2);
%                costt =0*j*(100-y.socn(1)).^2;
   q=2e3; k=8e2; j=0e0; %q=1e9; k=1.5e1;
     if(y.socn>80)
     q=1e0;k=1e1;j=1;    
     end
               costr=sum( fm(k*(y.socn-w.socr(w.ij)).^2)' +...
               q*abs(fm(y.jsn1)') );
               costt =j*(100-y.socn(1)).^2;
% cccv
%     k=1.7e1;j=1e-1; 
%                costr=sum( k*(y.socn-100).^2);
%                costt =0*j*(100-y.socn(1)).^2;

   cost=costr+costt;
end
function [c,ceq] = constrc(w,da,u)
      y = fhmss11t1assc(u,da,w);
     cr   = [-fm(y.socn);fm(y.socn)-100;-fm(y.socp);fm(y.socp)-100;
     fm(y.v)-4.6; -fm(y.v)+1.125];% 4.125 for cc
     ceqr = [];
     ct   = [];
     ceqt = [];
     c = [cr ct];
     ceq = [ceqr ceqt];
end
% function cost = costfe(w,da,x)
%         y = fhmss11t1assb(x,da,w);
%         v=y.y;
%         zz1=w.y(w.ij-w.N+1:w.ij);
%         zz=zz1-v;
%         costt=sum(1e5*zz.^2);
%         costr=0;
%         cost=costr+costt;
% end
% function [c,ceq] = constre(w,da,x)
%      y = fhmss11t1assb(x,da,w);
%      cr   = [ -fm(y.socn);fm(y.socn)-100;-fm(y.socp);fm(y.socp)-100;
%      fm(y.v)-4.5; -fm(y.v)+1.3;fm(-min(y.cel))+1e-3;-0*fm(y.jsn1)];
%      ceqr = [];
%      ct   = [];
%      ceqt = [];
%      c = [cr ct];
%      ceq = [ceqr ceqt];
% end

%%
function o=fm(i)
    o=reshape(i,1,length(i));
end
function d = store(d,d1,w)
d.x=[d.x;fm(d1.x)];d.est=[d.est,w.est];%d.x0=d1.x0;
d.socp=[d.socp,fm(d1.socp)];d.opns1=[d.opns1,fm(d1.opns1)];
d.cel=[d.cel,(d1.cel)];d.qf=[d.qf,d1.qf];
%d.jsn=[d.jsn,fm(d1.jsn)];
%d.u = [ d.u; d.ua(:,1) ];
d.v=[d.v,d1.v];d.jsn1=[d.jsn1,fm(d1.jsn1)];d.y=[d.y,fm(d1.y)];
d.socn=[d.socn,fm(d1.socn)];d.opn=[d.opn,fm(d1.opn)];
d.soh=[d.soh,fm(d1.soh)];d.up=[d.up,fm(d1.up)];
d.rfilm=[d.rfilm,d1.rfilm'];d.un=[d.un,fm(d1.un)];
end
function g = loop(w,d,d1,ij)
    g=1;
     if rem(ij,100)==0
         fprintf('\rij\t=\t%5.1f',ij); 
         fprintf('\rSOCN\t=\t%5.2f',d.socn(ij-w.N));
        fprintf('\rSOCR\t=\t%5.2f',d.socn(ij-w.N));
         fprintf('\rinput\t=\t%5.2f  \r',w.u(end));
     end     
     if (max(d1.v)>4.5)||(min(d1.v)<1.4)||(min(min(d1.cel))<1e-3)||min(abs(d1.socp-100))<.9||min(abs(d1.socn-100))<.9||ij==length(w.iter)-w.N      
         fprintf('\rij\t=\t%1f',ij); fprintf('\rVolt\t=\t%5.2f',d.v(end));
         fprintf('\rSOCN\t=\t%5.2f',d.socn(end));g=0; 
     end
end
function opt= optt(i,x0)
tol_opt = 2.16e-10; %.00216e-5
switch i
   case 1
    opt = optimset('Display','off','TolFun', tol_opt,...
    'MaxIter', 100000,'Algorithm', 'active-set',...
    'FinDiffType','central','RelLineSrchBnd', [],...
    'RelLineSrchBndDuration', [],'TolConSQP', tol_opt,'MaxFunEvals',1e4);
   case 2
    opt = optimset('Display','off','TolFun', tol_opt,...
    'MaxIter', 20000,'Algorithm', 'interior-point',...
    'AlwaysHonorConstraints','bounds','FinDiffType','forward',...
    'HessFcn', [],'Hessian', 'bfgs','HessMult', [],...
    'InitBarrierParam', 0.1,'InitTrustRegionRadius', ...
    sqrt(size(x0,1)*size(x0,2)),'MaxProjCGIter',...
    2*size(x0,1)*size(x0,2),'ObjectiveLimit', -1e10,...
    'ScaleProblem', 'obj-and-constr','SubproblemAlgorithm', ...
    'cg','TolProjCG', 1e-10,'TolProjCGAbs', 1e-10);
   case 3
    opt = optimset('Display','off','TolFun', tol_opt,...
    'MaxIter', 20000,'Algorithm', 'trust-region-reflective',...
    'Hessian', 'off','MaxPCGIter', max(1,floor(size(x0,1)*size(x0,2)/2)),...
    'PrecondBandWidth', 0,'tolPCG', 1e-1);
   case 4
    opt = optimset('Display','off','TolFun', tol_opt,...
    'MaxIter', 10000,'Algorithm', 'sqp',...
    'FinDifftype', 'forward','RelLineSrchBnd', [],...
    'RelLineSrchBndDuration', 1,'TolConSQP', 1e-30);
 end
end
