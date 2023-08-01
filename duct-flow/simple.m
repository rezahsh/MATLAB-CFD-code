function [uNew,vNew,pNew]=simple (Nx,Ny,Re,dx,dy,dt,u,v,p)
%SIMPLE method

tol1=1e-6;
itermax=5000;

iter1=0;
error1=1;
source=zeros(Nx+1,Ny+1);
uNew=u;
vNew=v;
pNew=p;

while  error1>tol1 & iter1<itermax
    iter1=iter1+1;
    pold = p;
    
    %--- guessed u,v based on guessed pressure p ---
    uNew = x_momentum (Nx,Ny,Re,dx,dy,dt,u,v,p); 
    vNew = y_momentum (Nx,Ny,Re,dx,dy,dt,u,v,p);

    %--- pressure correction source term ---
    for j=2:Ny
        for i=2:Nx        
            source(i,j)=(uNew(i,j)-uNew(i-1,j))*dx/dt...
                    +(vNew(i,j)-vNew(i,j-1))*dx^2/(dt*dy);
        end
    end

    %--- solve for the pressure correction using PSOR ---
    p_correction = psor(Nx,Ny,dx/dy,source);

    %--- velocities corrections ---
    u_correction=zeros(Nx+1,Ny+1);
    v_correction=zeros(Nx+1,Ny+1);

    for j=2:Ny
        for i=2:Nx        
            u_correction(i,j)=(p_correction(i,j)-p_correction(i+1,j))*dt/dx;
            v_correction(i,j)=(p_correction(i,j)-p_correction(i,j+1))*dt/dy;
        end
    end

    %--- update guessed p,u,v ---
    for j=1:Ny+1
        for i=1:Nx+1        
            u(i,j)=uNew(i,j)+u_correction(i,j);
            v(i,j)=vNew(i,j)+v_correction(i,j);
            p(i,j)=p(i,j)+p_correction(i,j);
        end
    end
    
    %--- convergence error ---
    error1=convergence(Nx+1,Ny+1,p,pold);
    if mod(iter1,10)==0
        disp(sprintf('              SIMPLE: iter=%d, error=%f',iter1,error1)) %display iter, error
    end
end
uNew=u;
vNew=v;
pNew=p;

if iter1 == itermax
    disp(sprintf('warning: itermax = iter1 in simple. iter1=%d, iter1max=%d, error=%f',...
         iter1,itermax,error1))
end

        