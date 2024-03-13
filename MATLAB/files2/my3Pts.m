%FILE: my3Pts.m
function [cent,rad]=my3Pts(z1,z2,z3) 
%Gives center of circle through z1,z2,z3
D=(abs(z1).^2).*(z2-z3)+ (abs(z2).^2).*(z3-z1)+(abs(z3).^2).*(z1-z2);
d=z1.*conj(z3-z2)+z2.*conj(z1-z3)+z3.*conj(z2-z1);
cent=D./d;
rad=abs(cent-z1);
end
