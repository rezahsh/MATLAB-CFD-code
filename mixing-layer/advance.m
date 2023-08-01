%solves all transport equations

i1=2:Nx-1;
j1=2:Ny-1;
 
% continuity
src1=zeros(Nx,Ny);nuT1=zeros(Nx,Ny);
var1=ones(Nx,Ny);
rho1 = transport (var1,rho,u,v,nuT1,src,1,1);

% turbulence model 
turbulence;

% x-momentum
src(i1,:)=-(p(i1+1,:)-p(i1-1,:))./(2*dx);
u1 = transport (u,rho,u,v,nuT,src)./rho1;
% y-momentum
src(:,j1)=-(p(:,j1+1)-p(:,j1-1))./(2*dy);
v1 = transport (v,rho,u,v,nuT,src)./rho1;

rho=rho1;
u=u1;
v=v1;

%eq. of state
p=rho.*RT;     

%B.C.
bc;
