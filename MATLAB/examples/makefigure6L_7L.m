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
deltv =  [0   ; delt ];
m     =  length(hr);
%
%%
CapDisk = (-2*pi)./log(tanh(hr./2));
UB      = sum(CapDisk)
% 
rho  =  @(x,y)(2*asinh(abs(x-y)/(sqrt(1-abs(x)^2)*sqrt(1-abs(y)^2))));
% 
etc     =   exp(i.*t);
etcp    =   i.* exp(i.*t);
%
%%
thtvj = [];
Rdv   =  [0.25,0.5,0.75,0.9,0.95,0.98];
for jj=1:length(Rdv)
Rd    =  Rdv(jj);
dd    =  sum(hr)+0.02;
thts  =  asin(sinh(dd/2)*(1-Rd^2)/(2*Rd));
thtv  =  [linspace(thts,0.25,21),linspace(0.26,pi/4,21),linspace(0.8,pi/2,21)].';
thtvj =  [thtvj thtv];
%
%
for kk=1:length(thtv)
    hz    =  [ Rd*exp(-i*thtv(kk))    Rd*exp(i*thtv(kk)) ];
    %
    for k=1:m
        [ecn(k),erd(k)]   = HypDisk(hz(k),hr(k));
        cr{k}   =  ecn(k)+erd(k).*exp(-i.*t);
        crp{k}  =      -i*erd(k).*exp(-i.*t);
    end
    %
    alpha  =  0;
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
    capv(kk,jj) = cap;
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
%
%
cc = Rd.*exp(i.*t);
plot(real(cc),imag(cc),':k','LineWidth',1.0);
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
if thtv(kk)==pi/4 & Rd==0.5
    print -depsc FigD2cirfig
end
str = sprintf('cap$=$ %1.6f ', cap);
title(str)
drawnow
% 
end
end
capv
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
clo = ['m','c','g','k','r','b'];
for jj=1:length(Rdv)
    plot(thtvj(:,jj),capv(:,jj),clo(jj),'LineWidth',1.5);
end
Leg=legend({'$R=0.25$','$R=0.5$','$R=0.75$','$R=0.9$','$R=0.95$',...
        '$R=0.98$'},'Interpreter','LaTeX',...
        'location','southeast');
Leg.AutoUpdate = 'off';
plot([0 pi],[UB UB],':k','LineWidth',2);
plot([0 pi],CapDisk,'--k','LineWidth',2);
set(gca,'FontSize',14)
axis square
axis([0 pi/2 2.0 4.5])
xticks([0 pi/4 pi/2])
xticklabels({'$0$','$\pi/4$','$\pi/2$'})
% xticks([-1:0.5:1])
% yticks([-1:0.5:1])
xlabel('$\theta$','Interpreter','latex')
ylabel('Capacity','Interpreter','latex')
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
print -depsc FigD2Capv
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
plot((thtv),(capv),'k-','LineWidth',1.5);
set(gca,'FontSize',14)
axis square
% axis([0 pi 2.0 4.5])
% xticks([0 pi/4 pi/2 3*pi/4 pi])
% xticklabels({'$0$','$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'})
% xticks([-1:0.5:1])
% yticks([-1:0.5:1])
xlabel('$\theta$','Interpreter','latex')
ylabel('Capacity','Interpreter','latex')
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.5; ax.MinorGridAlpha=0.5;
%%