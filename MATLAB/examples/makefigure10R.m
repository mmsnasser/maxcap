clear;
% 
addpath ../bie; addpath ../fmm; addpath ../files2; %addpath ../pcm
%%
% 
%%
Rd   =  0.75;
%
hrad =  [1.5; 0.2; 0.2; 0.2 ; 0.2 ; 0.2];
m    = length(hrad);
%
hcen = [ 0.75 +          0i
         0.064304 +    0.74724i
        -0.51118 +    0.54882i
        -0.75 + 2.5252e-07i
        -0.51118 -    0.54882i
         0.064303 -    0.74724i];
%
n       =   2^10;
t       =  (0:2*pi/n:2*pi-2*pi/n).';
%%
etc     =   exp(i.*t);
etcp    =   i.* exp(i.*t);
%
for k=1:m
    [ecen(k),erad(k)]   =  HypDisk(hcen(k),hrad(k));
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
[zet,zetp] = hyppolygonp(hcen,2^8);
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
    if hrad(k-1)==0.2
        plot(real(crv),imag(crv),'b-','LineWidth',1.5);
    else 
        plot(real(crv),imag(crv),'r-','LineWidth',1.5);
    end
end
%
plot(real(zhov),imag(zhov),'.r','LineWidth',1.5);
% 
plot(real(zet),imag(zet),':m','LineWidth',1.25);
%
plot(Rd*cos(t),Rd*sin(t),':k','LineWidth',1.25)
% plot(real(alpha),imag(alpha),'pr','LineWidth',1.5)
str = sprintf('cap$=$ %1.6f ', cap);
title(str)
set(gca,'FontSize',14)
axis square
axis([-1.05 1.05 -1.05 1.05])
xticks([-1:0.5:1])
yticks([-1:0.5:1])
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.5; ax.MinorGridAlpha=0.5;
%%
vcont  = [0.0,0.05,0.25,0.5,0.7,0.7723,0.7978,0.85,0.95];
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
contourf(real(zho),imag(zho),uzo,vcont)
colormap jet
caxis([0 1])
colorbar 
for k=1:m+1
    crv = et((k-1)*n+1:k*n); crv(n)=crv(1);
    plot(real(crv),imag(crv),'-k','LineWidth',1.5);
end
box on
axis equal
axis([-1.05 1.05 -1.05 1.05])
xticks([-1:0.5:1])
yticks([-1:0.5:1])
set(gca,'LooseInset',get(gca,'TightInset'))
set(gca,'XTick',[-1:1:8],'FontSize',18);
set(gca,'YTick',[0:1:6]);
set(gca,'LooseInset',get(gca,'TightInset'))
set(gcf,'Renderer','zbuffer')
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
surf(real(zho),imag(zho),uzo)
colormap jet
shading interp
% map = 0.8.*colormap+0.2.*ones(size(colormap));
% brighten(map,0.4)
caxis([0 1])
% vcont  = [0.2,0.5,1,1.5,2,2.5,3,4,5,5.5]; %level values
% % vcont  = [0.2,0.5,1,1.5,2,3,4,5]; %level values
% [cnt1,cnt2] =  contour(real(z),imag(z),disz,vcont,'b');
% % clabel(cnt1, 'manual','FontSize',15,'Color','k')
% clabel(cnt1,'FontSize',15,'Color','k')
box on
axis equal
axis([-1.05 1.05 -1.05 1.05 -0.05 1.05])
view([135 25]) 
set(gca,'LooseInset',get(gca,'TightInset'))
set(gca,'XTick',[-1:1:8],'FontSize',18);
set(gca,'YTick',[0:1:6]);
set(gca,'LooseInset',get(gca,'TightInset'))
set(gcf,'Renderer','zbuffer')
print -depsc Figu3D1p5
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
surf(real(zho),imag(zho),uzo)
colormap jet
shading interp
caxis([0 1])
% vcont  = [0.2,0.5,1,1.5,2,2.5,3,4,5,5.5]; %level values
% % vcont  = [0.2,0.5,1,1.5,2,3,4,5]; %level values
% [cnt1,cnt2] =  contour(real(z),imag(z),disz,vcont,'b');
% % clabel(cnt1, 'manual','FontSize',15,'Color','k')
% clabel(cnt1,'FontSize',15,'Color','k')
box on
axis equal
axis([-1.05 1.05 -1.05 1.05 -0.05 1.05])
view(2)
set(gca,'LooseInset',get(gca,'TightInset'))
set(gca,'XTick',[-1:1:8],'FontSize',18);
set(gca,'YTick',[0:1:6]);
set(gca,'LooseInset',get(gca,'TightInset'))
set(gcf,'Renderer','zbuffer')
print -depsc Figu2D1p5
%%
% r1=0.5
% z = [         0.75 +          0i
%       0.29343 +    0.69022i
%      -0.41524 +    0.62456i
%         -0.75 +  5.412e-09i
%      -0.41524 -    0.62456i
%       0.29343 -    0.69022i]
%
% 
% r1=1.5
% z = [         0.75 +          0i
%      0.064304 +    0.74724i
%      -0.51118 +    0.54882i
%         -0.75 + 2.5252e-07i
%      -0.51118 -    0.54882i
%      0.064303 -    0.74724i]
% 