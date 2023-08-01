%include all turbulence models

i1=2:Nx-1;
j1=2:Ny-1;

switch turbmodel
    case 'none'
        nuT=zeros(Nx,Ny);
    case 'algebraic_const'
        %ReT: turbulent Re, ReT=60-110 for mixing layer
        ReT=60;
        for j=1:Ny
            nuT(:,j)=rho(:,j).*Us.*delta(:)./ReT;
        end
    case 'algebraic_mixlength'
        A=0.247;
        alpha=0.071;
        lmix=alpha*(A*x);
        dudy(:,j1)=(u(:,j1+1)-u(:,j1-1))/(2*dy);
        for i=1:Nx
            nuT(i,:)=rho(i,:).*(lmix(i).^2).*abs(dudy(i,:));
        end
    case 'oneEq_k'
        A=0.247;
        alpha=0.071;
        sigmak=1;
        lmix=alpha*(A*x);
        dudy(:,j1)=(u(:,j1+1)-u(:,j1-1))/(2*dy);  
        % we assume du/dy is the dominant gradient in the production
        for i=1:Nx
            nuT(i,:)=rho(i,:).*(lmix(i).*ke(i,:).^(1/2))./sigmak;
            src(i,:)=rho(i,:).*((lmix(i).^2).*(dudy(i,:).^2)...
                -ke(i,:).^(3/2)./lmix(i));
        end
        
        ke1 = transport (ke,rho,u,v,nuT,src)./rho1;
        ke=max(ke1,0); 
        for i=1:Nx
            nuT(i,:)=rho(i,:).*lmix(i).*ke(i,:).^(1/2);
        end
        
    case 'spalart'
        cb1=0.1355;
        cb2=0.622;
        cv1=7.1;
        kappa=0.41;
        sigma=2/3;
        dudy(:,j1)=(u(:,j1+1)-u(:,j1-1))/(2*dy);
        dnudx=zeros(Nx,Ny);
        dnudy=zeros(Nx,Ny);
        nu0=nu1;
        dnudx(i1,:)=(nu0(i1+1,:)-nu0(i1-1,:))/(2*dx);
        dnudy(:,j1)=(nu0(:,j1+1)-nu0(:,j1-1))/(2*dy);
        % we assume du/dy is the dominant gradient in the rotation tensor
        for i=1:Nx
            src(i,:)=rho(i,:).*(2*cb1.*nu0(i,:).*dudy(i,:)+cb2/sigma.*(dnudx(i,:).*dnudx(i,:)...
                +dnudy(i,:).*dnudy(i,:)));
        end
        nu1 = transport (nu0,rho,u,v,rho.*nu0,src,1/sigma)./rho1;
        xx=nu1./nu;
        fv1=xx.^3./(xx.^3+cv1^3);
        nuT=rho.*fv1.*nu1;
        
    case 'k_epsilon'
        ce1=1.44;
        ce2=1.92;
        cmu=0.09;
        sigmak=1;
        sigmae=1.3;
        dudy(:,j1)=(u(:,j1+1)-u(:,j1-1))/(2*dy);  
        nuT=cmu*ke.^2./diss;
        nuT(find(diss<1e-20))=0;
        
        % we assume du/dy is the dominant gradient in the production
        src=rho.*(nuT.*dudy.^2.-diss);
        nuT=rho.*nuT/sigmak;
        ke1 = transport (ke,rho,u,v,nuT,src)./rho1;

        % we assume du/dy is the dominant gradient in the production
        src=rho.*diss./ke.*(ce1*nuT.*dudy.^2-ce2*diss);
        nuT=rho.*nuT/sigmae;
        diss1 = transport (diss,rho,u,v,nuT,src)./rho1;
        
        ke=max(ke1,0);
        diss=max(diss1,0); 
        
        nuT=rho.*cmu.*ke.^2./diss;
        nuT(find(diss<1e-20))=0;
end
