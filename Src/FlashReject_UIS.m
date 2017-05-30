function img_out = FlashReject_UIS( cfmimg, idx, thre )
%FLASHREJECT_UIS Summary of this function goes here
%   cfmimg  - the matrix of cfm data, in signed char.
%              [depth][width][frames]...
%   thre - the thre value of the algorithm, ranged 0~127
%   idx - the coming order of cfmimg.
%

% for the first 2 frames, let it go.
if idx <= 2
    img_out = cfmimg(:,:,idx);
    return;
end

imgcalc(:,:,1) = cfmimg(:,:,idx);
imgcalc(:,:,2) = cfmimg(:,:,idx-1);
imgcalc(:,:,3) = cfmimg(:,:,idx-2);

% get max velocity value
img_maxabs = max(abs(imgcalc),[],3);
lowidx = img_maxabs<thre;

%
[maxval,maxval_idx] = max(imgcalc,[],3);
% get 
for i=1:size(imgcalc,1)
    for j=1:size(imgcalc,2)
        i2min = 1:3;
        i2min(maxval_idx(i,j)) = [];
        img2min(i,j,:) = imgcalc(i,j,i2min);
    end
end
imgmin1 = int8(img2min(:,:,1));
imgmin2 = int8(img2min(:,:,2));

%% out put
img_out = cfmimg(:,:,idx);
img_out(lowidx) = Velocity8(imgmin1(lowidx)) + Velocity8(imgmin2(lowidx)-imgmin1(lowidx))/2;

end

