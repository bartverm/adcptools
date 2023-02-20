function [lamwgs,phiwgs,LL]=rd2wgs(x,y)
% RD2WGS: Converts Dutch rd-coordinates to GE wgs lat(Easting) long(Northing) coordinates % % USAGE: % [long,lat,LL]=rd2wgs(x,y) % % fortran 90 routine received from Peter Vermeulen % converted to Matlab by TO 090916 % % SEE ALSO: wgs2rd, kmlpath kmlpath2rd getDinoXSec % % TO 090916
if nargin<2, [lamwgs,phiwgs,LL]=selftest; return 
end
[phibes, lambes]=rd2bessel(x, y); 
[phiwgs, lamwgs]=bessel2wgs84(phibes,lambes);
N=length(x(:)); 
LL=cell(N,1);

% allocate for i=1:length(phibes(:)) LL{i,1}=sprintf('N %.7g, E %.7g',phiwgs(i),lamwgs(i)); end

function [phi,lambda] = rd2bessel(x,y) %convert xy to Bessel
x0 = 1.55e5; y0 = 4.63e5; k=.9999079; bigr = 6382644.571; m = .003773953832; n = 1.00047585668; e = .08169683122;
lambda0 = pi * .029931327161111111; b0 = pi * .28956165138333334;
d__1 = x - x0; d__2 = y - y0;
r = sqrt(d__1 .* d__1 + d__2 .* d__2);
warning off 
sa = (x - x0) ./ r; ca = (y - y0) ./ r; 
warning on 
sa(r==0)=0.0; ca(r==0)=0.0;
psi = atan2(r, k * 2. * bigr) * 2.;
cpsi = cos(psi); spsi = sin(psi);
sb = ca * cos(b0) .* spsi + sin(b0) * cpsi; d__1 = sb; cb = sqrt(1. - d__1 .* d__1); b = acos(cb); sdl = sa .* spsi ./ cb; dl = asin(sdl); lambda = dl / n + lambda0; w = log(tan(b / 2.0 + pi / 4.0)); q = (w - m) / n; phiprime = atan(exp(q)) * 2. - pi / 2.;
for i = 1:4
    dq = e / 2.0 * log((e * sin(phiprime) + 1.) ./ (1. - e * sin(phiprime))); 
    phi = atan(exp(q + dq)) * 2. - pi / 2.; 
    phiprime = phi; 
end
lambda = lambda / pi * 180.; phi = phi / pi * 180.;

function [phiwgs,lamwgs]=bessel2wgs84(phibes, lambes) % convert Bessel2 WGS84
a = 52.0; b = 5.0; c = -96.862; d = 11.714; e = 0.125; f = 1e-5; g = 0.329; h = 37.902; i = 14.667;
dphi = phibes - a; dlam = lambes - b; phicor = (c - dphi * d - dlam * e) * f; 
lamcor = (dphi * g - h - dlam * i) * f; 
phiwgs = phibes + phicor; lamwgs = lambes + lamcor;

function [long,lat,LL]=selftest
%Amersfoort en andere punten x = [155000; 244000; 93000; 98000; 177000]; y = [463000; 601000; 464000; 471000; 439000];
[long,lat,LL]=rd2wgs(x,y);