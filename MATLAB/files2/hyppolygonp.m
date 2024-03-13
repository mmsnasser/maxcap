function [et,etp] = hyppolygonp(ver,ns)
% Mohamed Nasser, 25-03-2019
%
n      =   length(ver)*ns;
t          =  (0:2*pi/n:2*pi-2*pi/n).';
% parametrization the polygon
m          = length(ver);ver(m+1)=ver(1); 
[s,sp]     =   deltw(t,m,3);
for k=1:m
    aa        = ver(k);  bb = ver(k+1);
    [cent(k),rd(k)] = my3Pts(aa,bb,bb/(abs(bb)^2));
    ang{k}    = carg([aa-cent(k),bb-cent(k)]);
end
for k=1:m
    alp(k)=ang{k}(1);
    bet(k)=ang{k}(2);
end
for k=1:m
    thet  =  (m*(bet(k)-alp(k))/(2*pi)).*s(1+(k-1)*n/m:k*n/m)+k*alp(k)-(k-1)*bet(k);
    thetp =  (m*(bet(k)-alp(k))/(2*pi)).*sp(1+(k-1)*n/m:k*n/m);
    et (1+(k-1)*n/m:k*n/m,1) =   cent(k)+rd(k).*exp(i.*thet);
    etp(1+(k-1)*n/m:k*n/m,1) =           rd(k).*exp(i.*thet).*(i.*thetp);
end
% 
end