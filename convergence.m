function error=convergence(A,Aold)
%convergence error

error=sum((A(:)-Aold(:)).^2);
norm=sum(Aold(:).^2);
error=sqrt(error/(norm+1e-20));  %1e-20 added to avoid devision by zero