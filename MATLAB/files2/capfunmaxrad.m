function capmin = capfunmaxrad(zvec)
% 
global hr Rd  et zet n fig1 cap
%
%%
zvec.'
m    = length(hr);
hz(1)  =  zvec(1);
for k=2:m
    hz(k)  =  zvec(2*k-2)+i*zvec(2*k-1);
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
map    =  PreImageDiskRad ([0;Lc],[0;Lk],1,n,1e-12,200);
%
zet     =   map.zet;
zetp    =   map.zetp;
% 
et      =   map.et;
etp     =   map.etp;
%
%
%    
deltv =[0];
for k=1:m
    deltv = [deltv; 1 ];
end
%
cent   =  map.cent; 
alpha  =  cent(1);
alphav =  cent(2:end);
%
[cap,a] =  capm(et,etp,alphav,deltv,m,alpha);
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
fig1=figure(1);
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