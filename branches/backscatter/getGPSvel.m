function gpsvel=getGPSvel(adcp,varargin)
if nargin>1
    assert(ischar(varargin{1}) && numel(varargin{1})==1 && any(strcmpi(varargin{1},{'b','e','s','i'})));
    destcor=varargin{1};
else
    destcor='e';
end

[X,Y]=utmADCP(adcp);
misal=getExtMisalign(adcp); % Get beam 3 misalignment

time=datenum(adcp.timeV);
time=(time-time(1))*24*3600;

dxdt=gradient(X)./gradient(time);
dydt=gradient(Y)./gradient(time);

gpsvel=zeros(size(adcp.btvel));
gpsvel(:,[1 2])=-[dxdt,dydt];

if any(strcmpi(destcor,{'b','s','i'}))
    adcp2=adcp;
    adcp2.btvel=gpsvel;
    [~, gpsvel]=corADCP(adcp2,destcor,'forceOrigin','e','UseExtHeading',true,'Beam3Misalign',misal);
end
