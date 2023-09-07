function [val,term,dir] = cutoff11t1(t,y,p)
% INPUTS
% t             The time at which the cut off condition has to be
%               calculated.
% y             The state at which the cut off condition has to be
%               calculated.
% OUTPUTS
% val           Value of the event. event(i) occurs if val(i) = 0.
% term          term(i) indicates if event(i) is terminal or not (i.e. if
%               time integration should be stopped if the event occurs or   
%               not).
% dir           dir(i) indicates the direction in which event (i) happens.
%%
qb=y(1:p.qbl);
cnn=y(p.qbl+1);cpp=y(p.qbl+2);
% phicn=y(p.qbl+3);
% phicp=y(p.qbl+4);
% jn=y(p.qbl+4+1:p.qbl+4+p.n);
% jp=y(p.qbl+4+p.n+1:p.qbl+4+p.n+p.p);
jn=p.u/(p.ln*p.f);
jp=-p.u/(p.lp*p.f);
%% electrode  model simplified
b=-(jn(end)*p.ln^2)/(p.f*p.dsn*6*p.nsn);
a= cnn-b/3;
csn=(a+b.*(p.zn'./p.ln).^2);
b1=-(jp(end)*p.lp^2)/(p.f*p.dsp*6*p.nsp);
a1= cpp-b1/3;
csp=fliplr(a1+b1.*((p.zp'/p.lp).^2));
%% electrolyte simplified
b=[0;0;0;0;qb];b1=b*2e9;
g=p.ca*b1;
cen=g(1)*((p.zn).^2)+g(2)*ones(1,length(p.zn));
ces=g(3)*((p.zs1).^2)+g(4)*p.zs1+g(5)*ones(1,length(p.zs1));
cep=g(6)*((p.l-p.zp1).^2)+g(7)*ones(1,length(p.zp1));
cel=[cep';ces';cen'];
% v=phicp-phicn-p.u*p.rc*p.a;
%%
[un,up,~,~] = ocp1(cnn/p.csn,cpp/p.csp);
%cep=cel(1:p.p);cen=cel(p.p+p.s+1:p.p+p.s+p.n);
ecdn=mean(p.kn*sqrt((cen.*csn').*(1-csn'/p.csn)));
ecdp=mean(p.kp*sqrt((cep.*csp').*(1-csp'/p.csp)));
opn=p.kb\asinh(jn/(2*ecdn));
opp=p.kb\asinh(jp/(2*ecdp));
phied=((p.ln+p.lp+2*p.ls)*p.u)/(2*p.ke) + p.kb*2*p.tp*p.ke*( log(cep(end))-log(cen(1))  );
v=opp-opn+phied+up-un-p.u*p.rc*p.a;
cnmin=min(csn)-p.snl(1);cnmax=p.snl(2)-max(csn);
cpmin=min(csp)-p.spl(1);cpmax=p.spl(2)-max(csp);
cemin=min(cel)-p.ecl(1);cemax=p.ecl(2)-max(cel);
vmin=v-p.vl(1);vmax=p.vl(2)-v;
val=[cnmin,cnmax,cpmin,cpmax,cemin,cemax,vmin,vmax];
%val=[cnmin,cnmax,cpmin,cpmax,cemin,cemax];
        term=ones(1,length(val));
      dir=-1*ones(1,length(val));     
%     % Calculate the voltage using the measurement function get_measurements
%     [V,~] = get_measurements(t,y,I,data,matrices_spm);
%     V = real(V);
% 
%     % Calculate the voltage 'window' left before the cut off voltage is
%     % reached. (When the window is <= 0, the V limits are exceeded)
%     Vspare_dis = V - V_limit(1); % V 'left' on discharge
%     Vspare_cha = V_limit(2)-V; % V 'left' on charge
% 
%     % If we reach the voltage limits, Vspare will become 0, and an event
%     % occurs
%     val = [Vspare_dis ; Vspare_cha];
%     % if term(i)=1, event(i) is terminal and time integration is halted
%     % when it happens. 
%     term = [1;1];
%     % event occurs only if 0 is reached from direction(i) (-1 means 'val'
%     % is decreasing and event occurs when val=0, 0 = from any direction, 1
%     % means 'val' is increasing and event when val=0 (val=negative))
%     dir = [-1;-1];
        % note: with current settings, an event occurs only if V is
        % decreasing and we cross the minimum voltage, or if V is
        % increasing and we cross the maximum voltage.
        % If you start charging at V < Vmin, no event will occur if
        % we cross Vmin. If you start discharging at V < Vmin, no event
        % will occur (eventually errors will be produced, e.g. if the
        % concentration becomes negative)
        % Similar for starting at V>Vmax (no event if you cross Vmax on
        % discharging, no event ever if you start charging)
end














