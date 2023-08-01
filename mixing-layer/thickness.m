function  yy = thickness( u,y,a,i )
%THICKNESS Summary of this function goes here
%   Detailed explanation goes here

jj=find(u(i,:)<=a);
if isempty(jj)
    yy=0;
else
    jj1=jj(end);
    if jj1==size(y,2)
        yy=y(size(y,2));
    else
        yy=y(jj1)+(y(jj1+1)-y(jj1))/(u(i,jj1+1)-u(i,jj1))*(a-u(i,jj1));
    end
end

end

