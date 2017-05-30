function velocity8 = velocity_Compress(data16, DRange, B)
% Compress 16 bit line data to 8 bit data

%% compression curve.
% 16 bit to 14 bit 
vmax = bitshift(1,14)-1;
dbase = 1.001 + B*.000215;
bs = 32767/(exp(127*log(dbase))-1);
%
x = 0:32767;
y1 = vmax*log(1+x/bs)/log(1+32767/bs);

%% apply dynamic range
% 14 bit to 8 bit
dr = DRange/100;
y = 127/dr*(y1/vmax+dr/2-1/2);
y(y<0) = 0;
y(y>127) = 127;


%% perform compression.
velocity8 = int8(sign(data16).*y(floor(abs(data16))+1));
end
