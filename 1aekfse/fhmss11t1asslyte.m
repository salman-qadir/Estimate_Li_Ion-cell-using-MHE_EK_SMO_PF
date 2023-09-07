function[y] = fhmss11t1asslyte(qb,p)
% t=p.tc;
u=p.ua;
u=reshape(u,length(u),1);
qb=(fm(qb))';
g=p.ca*([0;0;0;0;qb]);
      qfn=(qb(1)+p.xm*((1e-2*p.tp1*u.*ones(p.M,1))/p.f+p.den*2*g(2)*p.ln^2));
       qfs=(qb(2)+p.xm*ones(p.M,1)*p.des*2*g(5)*p.ls^2);
      qfp=(qb(3)+p.xm*((-1e-2*p.tp1*u.*ones(p.M,1))/p.f+p.dep*2*g(7)*p.lp^2));
% %          qfn=qb(1)+(2*t*((p.tp1*u)/p.f+p.den*2*g(2)*p.ln^2));
% %          qfs=qb(2)+(2*t*p.des*2*g(5)*p.ls^2);
% %          qfp=qb(3)+(2*t*(-(p.tp1*u)/p.f+p.dep*2*g(7)*p.lp^2));
          qf=([qfn(end);qfs(end);qfp(end)]);
% %rt=size(qf);g=p.ca*([zeros(4,rt(1,2));qf]);
%g=p.ca*([0;0;0;0;qf]);
y.cen=(g(2)'.*((p.zn).^2)+g(1)'.*ones(length(p.zn),1)');
y.ces=(g(5)'.*((p.zs).^2)+g(4)'.*p.zs+g(3)'.*ones(length(p.zs),1)');
y.cep=fliplr(g(7)'.*((p.zp).^2)+g(6)'.*ones(length(p.zp),1)');
%  y.cen=p.ce*ones(length(p.zn),1)';
%  y.ces=p.ce*ones(length(p.zs),1)';
%  y.cep=p.ce*ones(length(p.zp),1)';
% qf=qb;
y.cel=[y.cen';y.ces';y.cep'];
y.qf=fm(qf);
end
function o=fm(i)
o=reshape(i,1,length(i));
end
