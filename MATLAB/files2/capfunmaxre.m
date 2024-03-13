function cap = capfunmaxre(zvec)
% 
global hr Rd cap et n fig1 videow
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
etc     =   exp(i.*t);
etcp    =   i.* exp(i.*t);
%%
for k=1:m
    [ecn(k),erd(k)]   = HypDisk(hz(k),hr(k));
    cr{k}   =  ecn(k)+erd(k).*exp(-i.*t);
    crp{k}  =      -i*erd(k).*exp(-i.*t);
end
%
%
%
alpha  =  0.75i;
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
[cap,a] =  capm(et,etp,alphav,deltv,m,alpha);
cap
cap     = -cap;
%
% 
%%
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
% plot(real(alpha),imag(alpha),'pk','LineWidth',1.5)
str = sprintf('cap$=$ %1.6f ', -cap);
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