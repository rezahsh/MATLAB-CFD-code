%-----------------------------------------------
% Solution of 2D turbulent spatial mixing layer using RANS
% Prof. Reza Sheikhi
% University of Connecticut
% Spring 2015
%-----------------------------------------------
% Experiment data for the mixing layer case:
% H. W. Liepmann and J. Laufer (1947) "Investigation of Free Turbulent
% Mixing," NACA Technical Note No. 1257.

clear all;
clc;

global Nx Ny dx dy dt Ma gamma RT nu;
global uh ul;

%--- parameters ---
Nx=51;                  %number of nodes in x
Ny=101;                 %number of nodes in y
L=2;                    %length of the pipe
h=2;                    %height of the pipe
tol=1e-6;               %steady state tolerance
itermax=50000;          %steady state max number of iterations
gamma=1.4;              %heat capacity ratio (cp/cv)
Ma=0.2;                 %Ref. Mach number (for eq. of state)
nu=0.001;               %kinematic viscosity
uh=2;                   %top stream velocity
ul=1;                   %bottom stream velocity
vorThick=0.1;           %vorticity thickness of the layer at the inlet
Re=vorThick*(uh-ul)/nu; %Re number
intensity=0.1;          %turbulent intensity at the inflow

%--- turbulence model ---
% Models:
% - none: no turbulence model
% - algebraic_const: algebraic model with constant eddy viscosity
% - algebraic_mixlength: algebraic with mixing length model
% - oneEq_k: one equation model solving transport equation for turbulent K.E.
% - spalart: one equation Spalart-Allmaras model
% - k_epsilon: two equation standard k-epsilon model

turbmodel='k_epsilon';

%--- Grid ---
dx=L/(Nx-1);   
dy=h/(Ny-1);
x=0:dx:L;
y=0:dy:h;

%--- Stability criteria ---
fac=0.01;
limit1=0.5/sqrt(1/dx^2+1/dy^2)/nu;   %diffusion limit
umax=uh;                             % max u   
limit2=dx/umax;                      %convection limit
dt=fac*min([limit1 limit2]);         %min of convection and diffusion criteria with a saftety factor
Rec=umax*dx/nu;                      %cell Re number


%--- I.C. ---
x=0:dx:L;
y=0:dy:h;
u=zeros(Nx,Ny);
ycl=y((Ny-1)/2+1);       %location of the splitter plate
u(2,1:Ny)=0.5*(uh+ul)+0.5*(uh-ul).*tanh(20./vorThick.*(y(1:Ny)-ycl));
for i=3:Nx
    u(i,:)=u(2,:);
end

%--- init ---
init;
Us=uh-ul;
Uc=0.5*(uh+ul); 
RT=Us^2./gamma/Ma^2;
bc;
    
iter=0;
error=1;
%===============================
while  error>tol & iter<itermax   %time step loop
    iter=iter+1; 
    
    rho0=rho;
    u0=u;
    v0=v;
    
    %--- advance one time step ---
    advance;
    
    %--- update the thicknesses ---
    un=(u-ul)./Us;   %normalized velocity
for i=1:Nx
    y095(i)=thickness(un,y,0.95,i);    %y where u=ul+0.95Us
    y01(i) =thickness(un,y,0.1,i);     %y where u=ul+0.1Us
    delta(i)=y095(i)-y01(i);           %thickness of the shear layer
end
    %--- steady state error ---
    error1=convergence(u,u0);
    error2=convergence(v,v0);
    error=max([error1 error2]);
    
    if (mod(iter,100)==0)
        disp(sprintf('> Steady-state: iter=%d, dt=%f, error=%f',iter,dt,error)) %display iter, error
    end
        
end %time step loop end
%===============================

if iter == itermax
    disp(sprintf('warning: end but no steady-state yet!!! iter=%d, itermax=%d, error=%f'...
               , iter,itermax,error))
else
    disp('   ')
    disp('>>> steady-state is reached <<<') 
    disp(sprintf('  dt=%f, final time=%f',dt,iter*dt))
    disp('   ')
end

