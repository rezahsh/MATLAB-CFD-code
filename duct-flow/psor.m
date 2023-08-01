function p1 = psor(Nx,Ny,beta,source)
%PSOR method

relax=1.9;       %relaxation factor
iter1max=5000;   %max number of iterations
tol=1e-4;        %tolerance

error1=1;
iter1=0;

p1=zeros(Nx+1,Ny+1);  %initial guess p1=0 everywhere

while  error1>tol & iter1<iter1max
    pold=p1;
    iter1=iter1+1;    
    
    for j=2:Ny
        for i=2:Nx        
            p1(i,j)=p1(i,j)*(1-relax)+relax/(2*(1+beta^2))*(p1(i+1,j)...
                   +p1(i-1,j)+beta^2*(p1(i,j+1)+p1(i,j-1))-source(i,j));
        end
    end
    
    %--- B.C. ---
    % bottom
    for i=1:Nx+1
        p1(i,1)=p1(i,2);  %wall
    end

    % top
    for i=1:Nx+1
        p1(i,Ny+1)=p1(i,Ny); %wall
    end

    % left
    for j=1:Ny+1
        p1(1,j)=p1(2,j);  %inlet
    end

    % right
    for j=1:Ny+1
        p1(Nx+1,j)=2*p1(Nx,j)-p1(Nx-1,j);  %fully developed
    end
    
    %--- convergence error ---
    error1=convergence(Nx+1,Ny+1,p1,pold);
    %disp(sprintf('                    PSOR iter=%d, error=%f',iter1,error1)) %display iter, error
end

if iter1 == iter1max
    disp(sprintf('warning: iter1max = iter1 in psor. iter1=%d, iter1max=%d, error1=%f', iter1,iter1max,error1))
end
    
    