clc;clear

%parameters set manual
c = 1540;
f0 = 5e6;
fprf = 3.0e3;
theta = 15/180*pi;
scale = 23.1;  %cm/s 
DB_range = 12;
line_density = 2;
WF_Order = 4;
threshold = 50;

load cfm_colormap.mat
cfm_clrmap=cfm_clrmap0./256;
load ../data/c_linedata7.mat
[sampleCnt,lineCnt,enssembleNum] = size(frmI);

%% Polynomial WallFilter 
% x=reshape(1:enssembleNum,1,1,enssembleNum);
% n = WF_Order;
% frmI_filtered = zeros(sampleCnt, lineCnt, enssembleNum);
% frmQ_filtered = zeros(sampleCnt, lineCnt, enssembleNum);
% for j=1:line_density:lineCnt
%     for i=1:sampleCnt
%         p = polyfit(x,frmI(i,j,:),n);
%         fI = polyval(p,x,n);
%         frmI_filtered(i,j,:) =  frmI(i,j,:) - fI;
%         
%         p = polyfit(x,frmQ(i,j,:),n);
%         fQ = polyval(p,x,n);
%         frmQ_filtered(i,j,:) =  frmQ(i,j,:) - fQ;
%     end
% end
% save tempdata.mat frmI_filtered frmQ_filtered
% 
load tempdata.mat
 %% Autocorrelation
re_auto = zeros(sampleCnt,lineCnt);
im_auto = zeros(sampleCnt,lineCnt);
for j=1:line_density:lineCnt
     for idx=1:enssembleNum-1
         re_auto(:,j) = re_auto(:,j) + frmI_filtered(:,j,idx) .* frmI_filtered(:,j,idx+1) + frmQ_filtered(:,j,idx) .* frmQ_filtered(:,j,idx+1);
         im_auto(:,j) = im_auto(:,j) + frmI_filtered(:,j,idx) .* frmQ_filtered(:,j,idx+1) - frmQ_filtered(:,j,idx) .* frmI_filtered(:,j,idx+1);
     end
end
v_est = c*fprf*cos(theta)/(4*pi*f0).* atan(im_auto./(re_auto+1e-6));

%% Compress to 8bit(-128 ~ +127)
frmV = floor(32768.*100.*v_est./scale);
frmV(abs(frmV)<threshold) = 0;
velocity_8bit = velocity_Compress(frmV, DB_range, 0.01);
%% interpolation
for j=2:line_density:lineCnt-1
    velocity_8bit(:,j) = (velocity_8bit(:,j-1) + velocity_8bit(:,j+1))./2;
end
%% smooth
velocity_8bit = medfilt2(velocity_8bit,[5,5]);
%% scan conversion
for i=1:sampleCnt
    for j=1:lineCnt
        CoordX = j+floor((sampleCnt-i)*sin(theta));
        velocity_8bitSC(i,CoordX) = velocity_8bit(i,j);
    end
end

figure;image(floor(velocity_8bitSC)+128);
colormap(cfm_clrmap)
colorbar
     
        