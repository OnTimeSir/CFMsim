clc;clear

%parameters set manual
c = 1540;
f0 = 5e6;
fprf = 3.0e3;
theta = 15/180*pi;
scale = 23.1;%cm/s
DB_range = 12;
WF_Order = 4;

threshold = 0.02;

load cfm_colormap.mat
cfm_clrmap=cfm_clrmap0./256;
load ../data/c_linedata2.mat
[sampleCnt,lineCnt,enssembleNum] = size(frmI);

%% Polynomial WallFilter 
% x=reshape(1:enssembleNum,1,1,enssembleNum);
% n = WF_Order;
% frmI_filtered = zeros(sampleCnt, lineCnt, enssembleNum);
% frmQ_filtered = zeros(sampleCnt, lineCnt, enssembleNum);
% for j=1:lineCnt
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

load tempdata.mat
 %% Autocorrelation
re_auto = zeros(sampleCnt,lineCnt);
im_auto = zeros(sampleCnt,lineCnt);
for j=1:2:lineCnt
     for idx=1:enssembleNum-1
         re_auto(:,j) = re_auto(:,j) + frmI_filtered(:,j,idx) .* frmI_filtered(:,j,idx+1) + frmQ_filtered(:,j,idx) .* frmQ_filtered(:,j,idx+1);
         im_auto(:,j) = im_auto(:,j) + frmI_filtered(:,j,idx) .* frmQ_filtered(:,j,idx+1) - frmQ_filtered(:,j,idx) .* frmI_filtered(:,j,idx+1);
     end
end
v_est = c*fprf*cos(theta)/(4*pi*f0).* atan(im_auto./(re_auto+1e-6));

%% Compress to 8bit(-128 ~ +127)


%% interpolation
for j=2:2:lineCnt-1
    v_est(:,j) = (v_est(:,j-1) + v_est(:,j+1))./2;
end
%% scan conversion
for i=1:sampleCnt
    for j=1:lineCnt
        CoordX = j+floor((sampleCnt-i)*sin(theta));
        v_estSC(i,CoordX) = v_est(i,j);
    end
end

figure;imagesc(100.*v_estSC./scale);
colormap(cfm_clrmap)
colorbar
     
        