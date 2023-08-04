function x = symexp(y, c)
if nargin==1
    c = 1;
end
x = sign(y).*c.*(10.^abs(y)-1);
% np.sign(y)*l*(-1.0 + 10**np.abs(y))
end