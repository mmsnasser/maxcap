function [y,r] = HypDisk (x,M)
% HypDisk.m, 21-10-2020
% Based on: hypReuleaux20201021.nb (by Matti Vuorinen)
% y: is the Euclidean center of hyp disk B_{rho}(x,M),
% r: is its Euclidean radius
%  see HKVbook (4.20)
t = tanh(M./2);
y = x.*(1-t.^2)./(1-abs(x).^2.*t.^2);
r = (1-abs(x).^2).*t./(1-abs(x).^2.*t.^2);
end