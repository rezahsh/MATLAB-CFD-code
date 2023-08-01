%-----------------------------------------------
% Northeastern University
% MIE dept.
% ME 4565 - Matlab CFD code
% Solution of 2D incompressible Navier-Stokes in a 2D duct
% Developer Prof. Reza Sheikhi
% Fall 2010
%-----------------------------------------------
%The code is developed based on Chapter 6 of Patankar.
%It incorporates an explicit spatial second order, temporal first order
%accurate numerical scheme along with staggered grid arrangement.
%The boundary conditions are implemented using ghost cells.
%The pressure velocity link is handled using the Semi-Implicit Method for
%Pressure Linked Equations (SIMPLE) method.
%-----------------------------------------------

clear all;
clc;

%--- parameters ---
%Note: all variables are non-dimensionalized
%My computational grid is between [1,Nx+1],[1,Ny+1]. I solve for [2,Nx],[2,Ny].
% i=1,i=Nx+1,j=1,j=Ny+1 are my ghost cells for the B.C.

Nx=51;        %number of cells in x plus 1
Ny=21;        %number of cells in y plus 1
L=5;          %length of the pipe
h=1;          %height of the pipe
dx=1/(Nx-1);
dy=1/(Ny-1);
Re=10;        %Re number
tol=1e-6;      %steady state tolerance
itermax=5000;  %steady state max number of iterations

%--- calc. dt based on stability criteria ---
limit1=Re/2/(1/dx^2+1/dy^2); %diffusion limit
umax=1;                      %max(u(:));
limit2=dx/umax;              %convection limit
dt=0.9*min([limit1 limit2]);  %min of convection and diffusion criteria
                              %with a safety factor
%--- I.C. ---
u=zeros(Nx+1,Ny+1);
v=zeros(Nx+1,Ny+1);
p=zeros(Nx+1,Ny+1);
u(1,2:Ny)=1;         %inlet

iter=0;
error=1;
%===============================================================================
while  error>tol & iter<itermax   %time step loop
    iter=iter+1;

    %--- SIMPLE method ---
    [uNew,vNew,pNew]=simple(Nx,Ny,Re,dx,dy,dt,u,v,p);

    %--- steady state error ---
    error1=convergence(Nx+1,Ny+1,uNew,u);
    error2=convergence(Nx+1,Ny+1,vNew,v);
    error3=convergence(Nx+1,Ny+1,pNew,p);
    error=max([error1 error2 error3]);

    if mod(iter,10)==0
        disp(sprintf('> Steady-state: iter=%d, dt=%f, error=%f',iter,dt,error)) %display iter, error
    end

    %--- update old values ---
    u=uNew;
    v=vNew;
    p=pNew;

    %--- update dt: based on the stability criteria ---
    umax=max(u(:));
    limit2=dx/umax;              %convection limit in x
    vmax=max(v(:));
    limit3=dy/vmax;              %convection limit in y
    dt=0.9*min([limit1 limit2 limit3]);  %min of convection and diffusion criteria
                                         %with a saftety factor

end %time step loop end
%===============================================================================

if iter == itermax
    disp(sprintf('warning: end but no steady-state yet!!! iter=%d, itermax=%d, error=%f'...
               , iter,itermax,error))
else
    disp('   ')
    disp('>>> steady-state is reached <<<')
    disp(sprintf('  dt=%f, final time=%f',dt,iter*dt))
    disp('   ')
end

%--- interpolate face values to cell center location ---
[x_vec y_vec x_mat y_mat ucell vcell]=celldata(Nx,Ny,dx,dy,u,v);

%--- visualization ---
%figure;
%plot(ucell(???,:),y_vec);
%grid on; ylabel 'y'; xlabel 'u';title 'center';
%plot(x_vec,vcell(:,???));
%grid on; xlabel 'x'; ylabel 'v';title 'center';
%figure
%pcolor(x_mat,y_mat,ucell')
%shading interp;colormap(jet);colorbar;
%xlabel 'x'; ylabel 'y'; title 'u contours'
%figure;
%[sx,sy] = meshgrid(0.1:0.2:0.9,0.1:0.2:0.9);
%streamline(stream2(x_mat,y_mat,ucell',vcell',sx,sy));
%axis equal; axis tight;
%xlabel 'x'; ylabel 'y'; title 'streamlines'