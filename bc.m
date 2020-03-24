%boundary conditions 

% bottom : bottom stream (zero gradient & v=0)
rho(:,1)=rho(:,2);
p(:,1)=rho(:,1).*RT;
u(:,1)=u(:,2);
ke(:,1)=ke(:,2);
diss(:,1)=diss(:,2);
nu1(:,1)=nu1(:,2);  %for spalart model
v(:,1)=0;

% top : top stream (zero gradient & v=0)
rho(:,Ny)=rho(:,Ny-1);
p(:,Ny)=rho(:,Ny).*RT;
u(:,Ny)=u(:,Ny-1);
ke(:,Ny)=ke(:,Ny-1);
diss(:,Ny)=diss(:,Ny-1);
nu1(:,Ny)=nu1(:,Ny-1);  %for spalart model
v(:,Ny)=0;

% left : inflow
rho(1,:)=rho(2,:);
p(2,:)=rho(2,:).*RT;
v(1,:)=0;
u(1,:)=0.5*(uh+ul)+0.5*(uh-ul).*tanh(2./vorThick.*(y-ycl));
ke(1,:)=((u(1,:)*intensity).^2)*3/2;
nu1(1,:)=ke(1,:).^(1/2)*vorThick;        %for spalart model
diss(1,:)=0.09*ke(1,:).^(3/2)/vorThick;  %for k-epsilon

% right : fully developed
rho(Nx,:)=rho(Nx-1,:);
p(Nx,:)=rho(Nx,:).*RT;
u(Nx,:)=u(Nx-1,:);
v(Nx,:)=v(Nx-1,:);
ke(Nx,:)=ke(Nx-1,:);
diss(Nx,:)=diss(Nx-1,:);
nu1(Nx,:)=nu1(Nx-1,:);  %for spalart model