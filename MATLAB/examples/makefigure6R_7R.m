clear all;
% 
addpath ../bie; addpath ../fmm; addpath ../files2; %addpath ../pcm
%%
%
n       =   2^12;
t       =  (0:2*pi/n:2*pi-2*pi/n).';
%
% 
hr    =  [0.1 ; 0.1 ];
delt  =  [  1 ; 1   ];
deltv =  [0  ; delt ];
m     =  length(hr);
%
Rd    =  0.98;
dr    =  sum(hr)+0.02;
xr    =  tanh(dr/4);
hzv   =  [linspace(xr,0.5,26),linspace(0.51,Rd,26)];
%%
capd  = (-2*pi)./log(tanh(hr(1)./2));
UB    =  2*capd;
%%
% 
rho  =  @(x,y)(2*asinh(abs(x-y)/(sqrt(1-abs(x)^2)*sqrt(1-abs(y)^2))));
% 
etc     =   exp(i.*t);
etcp    =   i.* exp(i.*t);
%
for kk=1:length(hzv)
    hz    =  [ -hzv(kk)    hzv(kk) ];
    %
    for k=1:m
        [ecn(k),erd(k)]   = HypDisk(hz(k),hr(k));
        cr{k}   =  ecn(k)+erd(k).*exp(-i.*t);
        crp{k}  =      -i*erd(k).*exp(-i.*t);
    end
    %
    alpha  =  0.5i;
    %
    et  =  etc;  etp  =  etcp;  
    for k=1:m
        et    = [et  ; cr{k}];
        etp   = [etp ; crp{k}];
    end
    %
    alphav =  ecn;
    %
    [cap,~] = capm(et,etp,alphav,deltv,m,alpha);
    cap
    capv(kk,1) = cap;
%     
figure(1);
clf
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
k = 1;
crv = et((k-1)*n+1:k*n); crv(n)=crv(1);
plot(real(crv),imag(crv),'k-','LineWidth',1.5);
for k=2:m+1
    crv = et((k-1)*n+1:k*n); crv(n)=crv(1);
    if delt(k-1)==1
        plot(real(crv),imag(crv),'b-','LineWidth',1.5);
    else 
        plot(real(crv),imag(crv),'r-','LineWidth',1.5);
    end
end
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
if hzv(kk)==0.5
    print -depsc FigD2refig
end
str = sprintf('cap$=$ %1.6f ', cap);
title(str)
drawnow
% 
end
capv
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
plot([0 1],[UB UB],':k','LineWidth',2);
plot([0 1],[capd capd],'--k','LineWidth',2);
plot(hzv,capv,'-b','LineWidth',1.5);
set(gca,'FontSize',14)
axis square
axis([0 1 2.0 4.5])
% xticks([-1:0.5:1])
% yticks([-1:0.5:1])
xlabel('$x$','Interpreter','latex')
ylabel('Capacity','Interpreter','latex')
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
print -depsc FigD2Capre
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
plot(hzv,capv,'-b','LineWidth',1.5);
% plot(hzv,log(2-capv/capd),'-b','LineWidth',1.5);
% plot(hzv,capd+capd*tanh(hzv),':r','LineWidth',1.5);
% plot(hzv,2*capd-capd*exp(-4*hzv),'-r','LineWidth',1.5);
plot(hzv,capd+capd*(hzv./(hzv+0.1)).^(1/20),'-r','LineWidth',1.5);
plot([0 1],[2*capd 2*capd],':r','LineWidth',1.5);
plot([0 1],[capd capd],'--k','LineWidth',1.5);

set(gca,'FontSize',14)
axis square
xlabel('$x$','Interpreter','latex')
ylabel('Capacity','Interpreter','latex')
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.5; ax.MinorGridAlpha=0.5;
%%
