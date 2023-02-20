function y = symlog(x,c)
if nargin==1
    c = 1;
end
y = sign(x).*log10(1+abs(x./c)) ;
end