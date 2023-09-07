function [y] = lpf(u,e)

a=1-e; 

y=zeros(1,length(u));

y(1)=u(1);

for i=2:length(u)
    
    y(i)=e*y(i-1)+a*u(i);
    
end

end

