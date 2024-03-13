function capmin = capfunmaxslitre(zvec)
% 
global hr Rd  et zet n cap
%
%%
zvec.'
m    = length(hr);
for k=1:m
    hz(k)  =  zvec(k);
end
%
n
t       =  (0:2*pi/n:2*pi-2*pi/n).';
% 
for k=1:m
    [Lc(k,1),crd(k,1)]   = HypDisk(hz(k),hr(k)/2);
end
Lk = 2*crd;
%
thetk = zeros(size(hz));
% 
map = PreImageStrSlit(Lc,Lk,thetk,1,n,1e-12,200);
%
zeto     =   exp(i.*t);
zetop    =   i.* exp(i.*t);
%
zeti    =   map.zet;
zetip   =   map.zetp;
zet     =  [zeto;zeti];
% 
eti     =   map.et;
etip    =   map.etp;
% 
eto     =  (zeto.'+fcau(zeti,zetip,eti-zeti,zeto.',n,0)).';
etop    =   derfft(real(eto))+i*derfft(imag(eto));
% 
et      =  [eto ;eti ];
etp     =  [etop;etip];
%
alphas  =  0.5i;
alpha   = (alphas+fcau(zeti,zetip,eti-zeti,alphas,n,0)).';
%
%    
deltv =[0];
for k=1:m
    deltv = [deltv; 1 ];
end
%
cent    =  map.cent;
alphav  =  cent;
%
[cap,~] =  capm(et,etp,alphav,deltv,m,alpha);
cap
capmin  = -cap;
%
% 
%%
figure(2);
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
k=1;
crv = et(1+(k-1)*n:k*n);crv(n)=crv(1);
plot(real(crv),imag(crv),'k-','LineWidth',1.5)
hold on; box on
for k=2:m+1
    crv = et(1+(k-1)*n:k*n);crv(n)=crv(1);
    plot(real(crv),imag(crv),'-b','LineWidth',1.5)
end
%
plot(real(alpha),imag(alpha),'or')
plot(real(alphav),imag(alphav),'dr')
%
plot(Rd*cos(t),Rd*sin(t),':k','LineWidth',1.5)
% plot(real(alpha),imag(alpha),'pk','LineWidth',1.5)
set(gca,'FontSize',12)
axis square
axis([-1.05 1.05 -1.05 1.05])
xticks([-1:0.5:1])
yticks([-1:0.5:1])
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.5; ax.MinorGridAlpha=0.5;
drawnow
hold off
%%
figure(1);
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
k=1;
crv = zet(1+(k-1)*n:k*n);crv(n)=crv(1);
plot(real(crv),imag(crv),'k-','LineWidth',1.5)
hold on; box on
for k=2:m+1
    crv = zet(1+(k-1)*n:k*n);crv(n)=crv(1);
    plot(real(crv),imag(crv),'-b','LineWidth',1.5)
end
%
plot(Rd*cos(t),Rd*sin(t),':k','LineWidth',1.5)
% plot(real(alpha),imag(alpha),'pk','LineWidth',1.5)
str = sprintf('cap$=$ %1.6f ', cap);
title(str)
set(gca,'FontSize',12)
axis square
axis([-1.05 1.05 -1.05 1.05])
xticks([-1:0.5:1])
yticks([-1:0.5:1])
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.5; ax.MinorGridAlpha=0.5;
drawnow
hold off
%
% FF = getframe(fig1);   % Video ========
% writeVideo(videow,FF); % Video ========
%
end