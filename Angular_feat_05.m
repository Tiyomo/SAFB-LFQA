function wlbp = Angular_feat_05(lf)

% Version 1.0: For Zerui Yu "Spatial-Angular Features Based No-Reference Light
% Field Quality Assessment" submitted to ESWA.

% inputs:
% lf - the light field images with size [u, v, s, t, Ch]

% outputs:
% wlbp - the angular feats of downsampled distorted LF, size: (108, 1)

% Downsampling the LF by resizing SAIs
[U, V, S, T, Ch] = size(lf);
for uu = 1:U
    for vv= 1:V
        sai = squeeze(lf(uu, vv, :, :, :));
        sai = imresize(sai,0.5);
        newLF(uu, vv, :, :, :) = sai;
    end
end
lf = newLF;
[U, V, S, T, Ch] = size(newLF);

R_P         = {[1,8],[2,16],[3,24]};
C           = 1e-06;
mapp     = {load('mapping1.mat'),load('mapping2.mat'),load('mapping3.mat')};
Rlbp_feats  = [];    
Clbp_feats  = [];
for s = 1:length(R_P)
    radius      = R_P{s}(1); 
    Neighbors   = R_P{s}(2);
    mapping     = mapp{s}.mapping;
    
    Rdata       = 0;
    for row = 1:U
        for width = 1:S
            dis_Slice_epi = rgb2gray(squeeze(lf(row, :, width, :, :)));
            lbp_feat            = LBPV(dis_Slice_epi,radius,Neighbors,mapping);
            lbp_feat(isnan(lbp_feat)) = 0;
            tmp_data            = lbp_feat + C;
            Rweight(row,width)  = -sum(tmp_data.*log(tmp_data));
            Rdata               = Rdata + Rweight(row,width)*tmp_data;
        end
    end
    Rlbp_feats = [Rlbp_feats; Rdata'/(sum(Rweight(:)))];
    
    Cdata = 0;
    for col = 1:V
        for height = 1:T
            dis_Slice_epi = rgb2gray(squeeze(lf(:, col, :, height, :)));
            lbp_feat            = LBPV(dis_Slice_epi,radius,Neighbors,mapping);
            lbp_feat(isnan(lbp_feat)) = 0;
            tmp_data            = lbp_feat + C;
            Cweight(col,height)  = -sum(tmp_data.*log(tmp_data));
            Cdata               = Cdata + Cweight(col,height)*tmp_data;
        end
    end
    Clbp_feats = [Clbp_feats; Cdata'/(sum(Rweight(:)))];

end
wlbp = [Rlbp_feats; Clbp_feats];