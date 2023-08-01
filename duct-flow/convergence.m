function error=convergence(Nmax_x,Nmax_y,A,Aold)
%convergence error

error=0;
norm=0;
for j=1:Nmax_y
    for i=1:Nmax_x       
        error=error+(A(i,j)-Aold(i,j))^2;
        norm=norm+Aold(i,j)^2;
     end
end
error=sqrt(error/(norm+1e-20));  %1e-20 added to avoid devision by zero