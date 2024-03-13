clear all;
% 
addpath ../bie; addpath ../fmm; addpath ../files2; %addpath ../pcm
%%
% 
%%
hr    = 2.*[0.4; 0.4; 0.4; 0.4; 0.4]; 
%
n       =   2^10;
t       =  (0:2*pi/n:2*pi-2*pi/n).';
%
% 
%
%%
m    = length(hr);
hz   = [-0.9 -0.55 0 0.55 0.9];
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
%%
%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
k = 1;
crv = zet((k-1)*n+1:k*n); crv(n)=crv(1);
plot(real(crv),imag(crv),'k-','LineWidth',1.5);
for k=2:m+1
    crv = zet((k-1)*n+1:k*n); crv(n)=crv(1);
    plot(real(crv),imag(crv),'b-','LineWidth',2);
end
%
%
set(gca,'FontSize',14)
axis square
axis([-1.05 1.05 -1.05 1.05])
xticks([-1:0.5:1])
yticks([-1:0.5:1])
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
print -depsc FigSlit5figG
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
k = 1;
crv = et((k-1)*n+1:k*n); crv(n)=crv(1);
plot(real(crv),imag(crv),'k-','LineWidth',1.5);
for k=2:m+1
    crv = et((k-1)*n+1:k*n); crv(n)=crv(1);
    plot(real(crv),imag(crv),'b-','LineWidth',1.5);
end
%
%
set(gca,'FontSize',14)
axis square
axis([-1.05 1.05 -1.05 1.05])
xticks([-1:0.5:1])
yticks([-1:0.5:1])
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
print -depsc FigSlit5figO
%%
%
figure;
set(gca,'LooseInset',get(gca,'TightInset'))
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
for k=2:m+1
    crv = zet((k-1)*n+1:k*n); crv(n)=crv(1);
    plot(real(crv),imag(crv),'b-','LineWidth',2);
end
%
%
ax=gca; 
grid on; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',14)
xticks([-1:0.5:1])
yticks([-1:0.5:1])
axis equal
axis([-1.05 1.05 -0.55 0.55])
print -depsc FigSlit5figGu
%%
figure;
set(gca,'LooseInset',get(gca,'TightInset'))
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
for k=2:m+1
    crv = et((k-1)*n+1:k*n); crv(n)=crv(1);
    plot(real(crv),imag(crv),'b-','LineWidth',1.5);
end
%
%
ax=gca; 
grid on; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',14)
xticks([-1:0.5:1])
yticks([-1:0.5:1])
axis equal
axis([-1.05 1.05 -0.55 0.55])
print -depsc FigSlit5figOu
%%