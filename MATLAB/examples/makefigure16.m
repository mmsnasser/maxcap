clear;
% 
addpath ../bie; addpath ../fmm; addpath ../files2; %addpath ../pcm
%%
% 
%%
alpha  =  0;
%
Rd    =  0.75;
m     =  6;
hz    =  0.75*exp([0:2/m:2-2/m]*pi*i);
%
hrv   = [0.01,0.02,0.03,0.05,0.07,0.1,0.13,0.16,linspace(0.2,1.2,21)];
%
%
for jj=1:length(hrv)
    hr    =  hrv(jj)*ones(m,1);
    if hrv(jj)<0.8
        n=2^10
    else
        n=2^12
    end
    %
% 
%
n     =   2^10;
t     =  (0:2*pi/n:2*pi-2*pi/n).';
% 
etc     =   exp(i.*t);
etcp    =   i.* exp(i.*t);
% 
for k=1:m
    [ecn(k),erd(k)]   = HypDisk(hz(k),hr(k));
    cr{k}   =  ecn(k)+erd(k).*exp(-i.*t);
    crp{k}  =      -i*erd(k).*exp(-i.*t);
end
%
%    
et  =  etc;  etp  =  etcp;  deltv =[0];
for k=1:m
    et    = [et  ; cr{k}];
    etp   = [etp ; crp{k}];
    deltv = [deltv; 1 ];
end
%
alphav =  ecn;
%
[cap,~] =  capm(et,etp,alphav,deltv,m,alpha);
capv(jj,1) = cap;
%
R  =  log(coth(pi/cap));
Rv(jj,1) = R;
%
hyparea6(jj,1) =  6*4*pi*sinh(hr(1)/2)^2;
hyparea1(jj,1) =    4*pi*sinh(R/2)^2;
% 
hypprem6(jj,1) =  6*2*pi*sinh(hr(1));
hypprem1(jj,1) =    2*pi*sinh(R);
% 2*pi*csch(2*pi/cap)
% 
%
figure(1),clf
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
plot(real(etc),imag(etc),'k-','LineWidth',1.5)
hold on; box on
for k=1:m
    plot(real(cr{k}),imag(cr{k}),'b-','LineWidth',1.5);
end
%
plot(Rd*cos(t),Rd*sin(t),':k','LineWidth',1.5)
plot(real(alpha),imag(alpha),'pr','LineWidth',1.5)
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
end
%%
figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
plot(hrv,Rv,'-b','LineWidth',1.5);
axis square
axis([0 1.2 0 3])
xticks([0:0.2:1.2])
% yticks([-1:0.5:1])
xlabel('$r$','Interpreter','latex')
ylabel('$R$','Interpreter','latex')
set(gca,'FontSize',14)
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
print -depsc FigmdR
%%
figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
plot(hrv,hyparea6,'-b','LineWidth',1.5);
plot(hrv,hyparea1,'-r','LineWidth',1.5);
axis square
legend({'$E$','$B_\rho(0,R)$'},'Interpreter','LaTeX',...
        'location','northwest');
xlabel('$r$','Interpreter','latex')
ylabel('Hyp area','Interpreter','latex')
set(gca,'FontSize',14)
axis([0 1.2 0 40])
xticks([0:0.2:1.2])
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
print -depsc FigmdA
%%
figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
plot(hrv,hypprem6,'-b','LineWidth',1.5);
plot(hrv,hypprem1,'-r','LineWidth',1.5);
axis square
legend({'$E$','$B_\rho(0,R)$'},'Interpreter','LaTeX',...
        'location','northwest');
xlabel('$r$','Interpreter','latex')
ylabel('Hyp perimeter','Interpreter','latex')
set(gca,'FontSize',14)
axis([0 1.2 0 60])
xticks([0:0.2:1.2])
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; 
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
print -depsc FigmdP
%%