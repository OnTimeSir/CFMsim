clc;clear

%parameters set manual
c = 1540;
f0 = 5e6;
fprf = 3.0e3;
theta = 15/180*pi;
scale = 23.1;  %cm/s 
DB_range = 60;
line_density = 2;
WF_Order = 4;
ROthre = 50;

FRthre = 32;
load cfm_colormap.mat
cfm_clrmap=cfm_clrmap0./256;

pathname = '../data/';
for idx=1:10
    filename = ['c_linedata',num2str(idx),'.mat'];
    load ([pathname filename]);
    [sampleCnt,lineCnt,enssembleNum] = size(frmI);

%% Polynomial WallFilter 
%     x=reshape(1:enssembleNum,1,1,enssembleNum);
%     n = WF_Order;
%     frmI_filtered = zeros(sampleCnt, lineCnt, enssembleNum);
%     frmQ_filtered = zeros(sampleCnt, lineCnt, enssembleNum);
%     for j=1:line_density:lineCnt
%         for i=1:sampleCnt
%             p = polyfit(x,frmI(i,j,:),n);
%             fI = polyval(p,x,n);
%             frmI_filtered(i,j,:) =  frmI(i,j,:) - fI;
% 
%             p = polyfit(x,frmQ(i,j,:),n);
%             fQ = polyval(p,x,n);
%             frmQ_filtered(i,j,:) =  frmQ(i,j,:) - fQ;
%         end
%     end
%     save ([pathname 'tempdata',num2str(idx),'.mat'], 'frmI_filtered', 'frmQ_filtered')
    
    load ([pathname 'tempdata',num2str(idx),'.mat'])
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

    if (line_density~=1)
        for j=2:line_density:lineCnt-1
            v_est(:,j) = (v_est(:,j-1)+v_est(:,j+1))./2;
        end
    end

%% Compress to 8bit(-128 ~ +127)
    frmV = floor(32768.*100.*v_est./scale);
%     frmV(abs(frmV)<ROthre) = 0;
    velocity_frms(:,:,idx) = velocity_Compress(frmV, DB_range, 0.01);
%% FlashReject
    velocity_8bit = FlashReject_UIS( velocity_frms, idx, FRthre);
%     a = 0.2;
%     if idx>1
%         R0_Prev = velocity_frms(:,:,idx-1);
%         velocity_8bit = velocity_frms(:,:,idx) *a + R0_Prev*(1-a);
%     else
%         velocity_8bit = velocity_frms(:,:,idx);
%     end

%% smooth
    velocity_8bit = medfilt2(velocity_8bit,[5,5]);
    
%% scan conversion
    for i=1:sampleCnt
        for j=1:lineCnt
            CoordX = j+floor((sampleCnt-i)*sin(theta));
            velocity_8bitSC(i,CoordX) = velocity_8bit(i,j);
        end
    end
%% display images
    figure(2);
    image(double(velocity_8bitSC)+128);
    colormap(cfm_clrmap)
    colorbar('Ticks',[1,256],...
         'TickLabels',{['-',num2str(scale),' cm/s'],['+',num2str(scale),' cm/s']})
end     
        