function [d1,V, flag, out] = estimate11t1a(d,w)
    w.tc=w.t2/2;w.xm=w.xm2;w.xxm=w.xxm2;w.t=w.t2;w.mk=w.mk2;
    A   = [];b   = [];
    Aeq = [];beq = []; % saving flag and  out in struct is good.
    [x,V,flag,out]=fmincon(@(x) costfe(w,d,x),d.x0,A,...
    b,Aeq,beq,w.lb,w.ub,@(x) constre(w,d,x), w.opt);
    %%
    d.x0=repmat(x((w.N-1)*w.no+1:w.N*w.no),1,w.N);
    d1=fhmss11t1assd(x((w.N-1)*w.no+1:w.N*w.no),...
    d.qf(:,end),w);
end
function cost = costfe(w,da,x)
        y = fhmss11t1assb(x,da,w);
        v=y.y;
        zz1=w.y(w.ij-w.N+1:w.ij);
        zz=zz1-v;
        costt=sum(1e5*zz.^2);
        costr=0;
        cost=costr+costt;
end
function [c,ceq] = constre(w,da,x)
     y = fhmss11t1assb(x,da,w);
     cr   = [ -fm(y.socn);fm(y.socn)-100;-fm(y.socp);fm(y.socp)-100;
     fm(y.v)-4.5; -fm(y.v)+1.3;fm(-min(y.cel))+1e-3;-0+fm(y.jsn)];
     ceqr = [];
     ct   = [];
     ceqt = [];
     c = [cr ct];
     ceq = [ceqr ceqt];
end
%%
function o=fm(i)
    o=reshape(i,1,length(i));
end