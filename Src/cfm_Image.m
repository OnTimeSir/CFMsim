clc;clear
%% parameters set manual
theta = 15/180*pi; %deflection angle
WF_Order = 4;      %wall filter order
line_density = 2;  %1-high density;2-middle density;4-low density
scale = 23.1;      %the upper limit velocity (cm/s) 
Reverse = 1;       %1-left deflection;0-right deflection 
FenceRemove = 1;   %if remove fence line

DB_range_min = -12;
DB_range_max = 12;
RO_thre = 50;

Modulus_thre = 32000;  %modulus threshold

FR_thre = 32;    %flash reject threshold
%% 
load cfm_colormap.mat
cfm_clrmap=cfm_clrmap0./256;
if (Reverse==1)
    cfm_clrmap = flipud(cfm_clrmap);
end

if (FenceRemove==1)
    ens_start = 2;
else
    ens_start = 1;
end

pathname = '../data/';
for frm_idx=1:10
    filename = ['c_linedata',num2str(frm_idx),'.mat'];
    load ([pathname filename]);
    [sampleCnt,lineCnt,enssembleNum] = size(frmI);
%% PreEstimate
     %remove DC
     SumI = zeros(sampleCnt,lineCnt);
     SumQ = zeros(sampleCnt,lineCnt);
     for idx=ens_start:enssembleNum
         SumI = SumI + frmI(:,:,idx);
         SumQ = SumQ + frmQ(:,:,idx);
     end
     SumI = SumI./(enssembleNum-ens_start+1);
     SumQ = SumQ./(enssembleNum-ens_start+1);
     for idx = 1:enssembleNum
         frmI_rmvDC(:,:,idx) = frmI(:,:,idx)-SumI;
         frmQ_rmvDC(:,:,idx) = frmQ(:,:,idx)-SumQ;
     end
     %Calculate Modulus
    PreEst_re = zeros(sampleCnt,lineCnt);
    PreEst_im = zeros(sampleCnt,lineCnt);
    for j=1:line_density:lineCnt
        for idx=ens_start:enssembleNum-1
            PreEst_re(:,j) =  PreEst_re(:,j) + frmI_rmvDC(:,j,idx) .* frmI_rmvDC(:,j,idx+1) + frmQ_rmvDC(:,j,idx) .* frmQ_rmvDC(:,j,idx+1);
            PreEst_im(:,j) =  PreEst_im(:,j) + frmI_rmvDC(:,j,idx) .* frmQ_rmvDC(:,j,idx+1) - frmQ_rmvDC(:,j,idx) .* frmI_rmvDC(:,j,idx+1);
        end
    end
    Modulus = ;
%% Polynomial WallFilter 
    x=reshape(1:enssembleNum,1,1,enssembleNum);
    n = WF_Order;
    frmI_filtered = zeros(sampleCnt, lineCnt, enssembleNum);
    frmQ_filtered = zeros(sampleCnt, lineCnt, enssembleNum);
    for j=1:line_density:lineCnt
        for i=1:sampleCnt
            if (Modulus(i,j)<Modulus_thre)
                p = polyfit(x,frmI(i,j,:),n);
                fI = polyval(p,x,n);
                frmI_filtered(i,j,:) =  frmI(i,j,:) - fI;
    
                p = polyfit(x,frmQ(i,j,:),n);
                fQ = polyval(p,x,n);
                frmQ_filtered(i,j,:) =  frmQ(i,j,:) - fQ;
            end
        end
    end
    save ([pathname 'tempdata',num2str(frm_idx),'.mat'], 'frmI_filtered', 'frmQ_filtered')
    
    load ([pathname 'tempdata',num2str(frm_idx),'.mat'])

%% Autocorrelation and Normalized to 8bit(-128 ~ +127)
    re_auto = zeros(sampleCnt,lineCnt);
    im_auto = zeros(sampleCnt,lineCnt);
    for j=1:line_density:lineCnt
         for idx=ens_start:enssembleNum-1
             re_auto(:,j) = re_auto(:,j) + frmI_filtered(:,j,idx) .* frmI_filtered(:,j,idx+1) + frmQ_filtered(:,j,idx) .* frmQ_filtered(:,j,idx+1);
             im_auto(:,j) = im_auto(:,j) + frmI_filtered(:,j,idx) .* frmQ_filtered(:,j,idx+1) - frmQ_filtered(:,j,idx) .* frmI_filtered(:,j,idx+1);
         end
    end
    v_est = atan2(im_auto(:,1:line_density:end),re_auto(:,1:line_density:end)).*127/pi;
    velocity_frms(:,:,frm_idx) = int8(round(v_est));
%% FlashReject
    velocity_8bit = FlashReject_UIS( velocity_frms, frm_idx, FR_thre);

%% smooth
%     velocity_8bit = medfilt2(velocity_8bit,[5,5]);
    
%% scan conversion
    for i=1:sampleCnt
        for j=1:floor(lineCnt/line_density)
            CoordX = j+floor((sampleCnt-i)*sin(theta));
            velocity_8bitSC(i,CoordX) = velocity_8bit(i,j);
        end
    end
%% display images
    figure(1);
    image(double(velocity_8bitSC)+128);
    colormap(cfm_clrmap)
    colorbar('Ticks',[1,255],...
         'TickLabels',{['-',num2str(scale),' cm/s'],['+',num2str(scale),' cm/s']})
end     
        