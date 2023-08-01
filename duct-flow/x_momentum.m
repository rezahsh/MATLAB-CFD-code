function uNew = x_momentum (Nx,Ny,Re,dx,dy,dt,u,v,p)
%solves x-momentum

uNew=zeros(Nx+1,Ny+1);

%--- Average v velocity ---
for j=2:Ny
    for i=2:Nx-1
        v_ave(i,j)=0.25*(v(i,j)+v(i+1,j)+v(i,j-1)+v(i+1,j-1));
    end
end

%--- solve for uNew ---
for j=2:Ny
    for i=2:Nx-1
        conv=-u(i,j)*(u(i+1,j)-u(i-1,j))/(2*dx)...
             -v_ave(i,j)*(u(i,j+1)-u(i,j-1))/(2*dy);
         
        press=-(p(i+1,j)-p(i,j))/dx;
        
        diff=1/Re*((u(i+1,j)-2*u(i,j)+u(i-1,j))/(dx^2)...
                  +(u(i,j+1)-2*u(i,j)+u(i,j-1))/(dy^2));
                  
        uNew(i,j)= u(i,j)+(conv+diff+press)*dt;
    end
end

%--- B.C. ---
% bottom
for i=1:Nx+1
    uNew(i,1)=-uNew(i,2); %wall
end

% top
for i=1:Nx+1
    uNew(i,Ny+1)=-uNew(i,Ny); %wall
end

% left
for j=1:Ny+1
    uNew(1,j)=1; %inlet
end

% right
for j=1:Ny+1
    uNew(Nx,j)=uNew(Nx-1,j);  %fully developed
    uNew(Nx+1,j)=0;
end