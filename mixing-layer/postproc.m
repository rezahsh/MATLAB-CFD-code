%post-processing

% Experiment data for the mixing layer case:
% H. W. Liepmann and J. Laufer (1947) "Investigation of Free Turbulent
% Mixing," NACA Technical Note No. 1257.

%Obtained from Fig. 13 of Liepmann & Laufer
ye=[-2.4996643 -2.1687005 -1.942737 -1.7230953 -1.5001348 -1.276977 ...
    -1.0599434 -0.83631176 -0.6124827 -0.388476 -0.17399146 0.04693332...
    0.27753818 0.4953218 0.71917063 0.9395229 1.1562208 1.3758823 1.592205];


ue=[0.023191182 0.034588337 0.043536514 0.055230003 0.07646433...
    0.11134053 0.16260391 0.23022062 0.31147918 0.40501544 0.49994126...
    0.6003069 0.7101965 0.813299 0.8959217 0.95672596 0.98479813...
    0.99785584 1.00000];


% similarity variable
for i=1:Nx
    y05(i) =thickness(un,y,0.5,i);    
    yn(i,:)=(y(:)-y05(i));    
end

figure(1);
% clf

%xy plot the u vs. y form experiment
momthick=trapz(ye,ue.*(1.-ue));  %momentum thickness of the experiment
%the profile is adjusted to get u(y=0)=0.5 (this is to compensate for
%digitizing error)
plot(ue,ye*momthick+0.05,'ro');hold on;grid on


%xy plot the normalized u vs. y form CFD
 m=40;  %flow is self-similar at this point
 plot(un(m,:),yn(m,:),'k--')
 xlabel 'normalized u'
 ylabel 'normalized y'
 axis([0 1 -1 1])
 legend( 'Exp.', 'CFD','Location',"best")
%<u'v'>
 uv=nuT.*dudy;
 figure(2);
 plot(y,uv(40,:),'m-')
 xlabel 'y'
 ylabel(['$\overline{u \prime v \prime}$'],'interpreter','latex')


