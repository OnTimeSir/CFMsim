function img8 = UIS_Compress( data16, DRange, B)
% Compress 16 bit line data to 8 bit data,
% using compression curve specified by va/vb/vmax.
%
% Detailed explanation goes here:
% 

%% compression curve.
% 16 bit to 14 bit 
vmax = bitshift(1,14)-1;
dbase = 1.001 + B*.000215;
bs = 65535/(exp(255*log(dbase))-1);
%
x = 0:65535;
y1 = vmax*log(1+x/bs)/log(1+65535/bs);

%% apply dynamic range
% 14 bit to 8 bit
dr = DRange/100;
y = 255/dr*(y1/vmax+dr/2-1/2);
y(y<0) = 0;
y(y>255) = 255;


%% perform compression.
data16(data16<0) = 0;
data16(data16>65535) = 65535;

img8 = y(floor(data16)+1);
end
