function [x_vector,y_vector,x_matrix,y_matrix,u_cell,v_cell]=celldata(Nx,Ny,dx,dy,u,v)
%interpolate face values to cell center for vizualization
%cell centers start at (dx/2,dy/2) and end at (L-dx/2,L-dy/2)

%--- cell center x coordinate as vectors ---
for i=1:Nx-1
   x_vector(i)=(i-1)*dx+dx/2;
end

%--- cell center y coordinate as vectors ---
for j=1:Ny-1
   y_vector(j)=(j-1)*dy+dy/2;
end
%--- cell center x,y coordinates as matices ---
[x_matrix y_matrix]=meshgrid(x_vector,y_vector);

%--- interpolate u ---
for j=2:Ny
    for i=1:Nx-1
        u_cell(i,j-1)=0.5*(u(i+1,j)+u(i,j));
    end
end

%--- interpolate v ---
for j=1:Ny-1
    for i=2:Nx
        v_cell(i-1,j)=0.5*(v(i,j+1)+v(i,j));
    end
end
