function fnew = transport (f,rho,u,v,nuT,src,c1,solve_continuity)
%solves a generic transport eq. using explicity second order central space
%and first order time accurate scheme
%for continuity set solve_continuity=1

global Nx Ny dx dy dt Ma gamma RT nu;

if nargin<8
    solve_continuity=0;
end
if nargin<7
    solve_continuity=0;
    c1=1;
end 

fnew=f;
convx=zeros(Nx,Ny);
convy=zeros(Nx,Ny);
diffxy=zeros(Nx,Ny);

j=2:Ny-1;
i=2:Nx-1;

%convection
convx(i,j)=(rho(i+1,j).*u(i+1,j).*f(i+1,j)...
    -rho(i-1,j).*u(i-1,j).*f(i-1,j))./(2*dx);
convy(i,j)=(rho(i,j+1).*v(i,j+1).*f(i,j+1)...
    -rho(i,j-1).*v(i,j-1).*f(i,j-1))./(2*dy);

%diffusion
if solve_continuity==0
    diffxy(i,j)=c1*(nu+nuT(i,j)).*((f(i+1,j)-2.*f(i,j)+f(i-1,j))./(dx^2)...
        +(f(i,j+1)-2.*f(i,j)+f(i,j-1))./(dy^2));
else
    diffxy(i,j)=0;
end

fnew(i,j)= f(i,j)+(-convx(i,j)-convy(i,j)+diffxy(i,j)+src(i,j)).*dt;
       
