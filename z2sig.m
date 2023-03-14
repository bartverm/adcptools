function sig = z2sig(z, zb, eta)
% Following the convention of Vermeulen: Sigma between 0 (bed ) and 1
% (surface)
% To be used with CTD data.

sig = (z - zb)./(eta - zb);
end