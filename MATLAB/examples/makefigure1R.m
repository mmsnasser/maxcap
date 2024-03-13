clear;
% 
addpath ../bie; addpath ../fmm; addpath ../files2; %addpath ../pcm
%%
% 
%%
Rd   =  0.75;
%
hrad =  [0.2; 0.2; 0.2; 0.2 ; 0.2 ; 0.2];
m    = length(hrad);
%
zf = [   0.749999973835091 + 0.000000000000000i
  0.375000125816403 + 0.649518949985617i
 -0.374999690443153 + 0.649519201348318i
 -0.749999973834973 + 0.000000433271344i
 -0.375000113245405 - 0.649518957243438i
  0.374999867368982 - 0.649519099200461i];
%
n       =   2^10;
t       =  (0:2*pi/n:2*pi-2*pi/n).';
%%
etc     =   exp(i.*t);
etcp    =   i.* exp(i.*t);
%
for k=1:m
    [ecen(k),erad(k)]   =  HypDisk(zf(k),hrad(k));
    cr{k}   =  ecen(k)+erad(k).*exp(-i.*t);
    crp{k}  =      -i*erad(k).*exp(-i.*t);
end
%
alpha  =  0;
%
et  =  etc;  etp  =  etcp;  deltv =[0];
for k=1:m
    et    = [et  ; cr{k}];
    etp   = [etp ; crp{k}];
    deltv = [deltv; 1 ];
end
%
alphav =  ecen;
%
% 
%%
% 
[xh  , yh]  =  meshgrid([-0.999:0.001:0.999],[-0.999:0.001:0.999]);
zho         =  xh+i.*yh;
[mh,nh]     =  size(zho);
%
zho(abs(zho)>1-1e-3)=NaN+i*NaN;
%
for k=1:m
    zho(abs(zho-ecen(k))<erad(k)+1e-3)=NaN+i*NaN;
end
%
zhovo  =  zho(:);
zhov   =  zhovo(abs(zhovo)>=0).';
% 
%%
[cap,ak,uzov] =  capm(et,etp,alphav,deltv,m,alpha,zhov);
%%
% 
uzovo  =  NaN.*ones(size(zhovo));
uzovo(abs(zhovo)>=0) = uzov.';
uzo    =  NaN.*ones(size(zho));
uzo(:) =   uzovo;
%
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
surf(real(zho),imag(zho),uzo)
colormap jet
shading interp
caxis([0 1])
box on
axis equal
axis([-1.05 1.05 -1.05 1.05 -0.05 1.05])
view([135 25]) 
set(gca,'LooseInset',get(gca,'TightInset'))
set(gca,'XTick',[-1:1:8],'FontSize',18);
set(gca,'YTick',[0:1:6]);
set(gca,'LooseInset',get(gca,'TightInset'))
set(gcf,'Renderer','zbuffer')
print -depsc Fig3Dintf
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
surf(real(zho),imag(zho),uzo)
colormap jet
shading interp
caxis([0 1])
str = sprintf('cap$=$ %1.6f ', cap);
box on
axis equal
axis([-1.05 1.05 -1.05 1.05 -0.05 1.05])
view(2)
set(gca,'LooseInset',get(gca,'TightInset'))
set(gca,'XTick',[-1:1:8],'FontSize',18);
set(gca,'YTick',[0:1:6]);
title(str,'FontSize',14)
set(gca,'LooseInset',get(gca,'TightInset'))
set(gcf,'Renderer','zbuffer')
print -depsc Fig2Dintf
%%
