function vNew = y_momentum (Nx,Ny,Re,dx,dy,dt,u,v,p)
%solves y-momentum

vNew=zeros(Nx+1,Ny+1);

%--- Average v velocity ---
for j=2:Ny-1
    for i=2:Nx
        u_ave(i,j)=0.25*(u(i,j)+u(i,j+1)+u(i-1,j)+u(i-1,j+1));
    end
end

%--- solve for vNew ---
for j=2:Ny-1
    for i=2:Nx
        conv=-u_ave(i,j)*(v(i+1,j)-v(i-1,j))/(2*dx)...
             -v(i,j)*(v(i,j+1)-v(i,j-1))/(2*dy);
         
        press=-(p(i,j+1)-p(i,j))/dy;
        
        diff=1/Re*((v(i+1,j)-2*v(i,j)+v(i-1,j))/(dx^2)...
                  +(v(i,j+1)-2*v(i,j)+v(i,j-1))/(dy^2));
              
        vNew(i,j)= v(i,j)+(conv+diff+press)*dt;     
    end
end

%--- B.C. ---
% bottom
for i=1:Nx+1
    vNew(i,1)=0;  %wall
end

% top
for i=1:Nx+1
    vNew(i,Ny)=0;
    vNew(i,Ny+1)=0;  %wall
end

% left
for j=1:Ny+1
    vNew(1,j)=-vNew(2,j);  %inlet
end

% right
for j=1:Ny+1
    vNew(Nx+1,j)=vNew(Nx,j); %fully developed
end