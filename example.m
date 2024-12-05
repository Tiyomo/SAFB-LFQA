LFpath = 'HEVC_Bikes_29.bmp';
LF = WIN5ReadLF(LFpath);
feat0 = Angular_feat(LF);
feat1 = Angular_feat_05(LF);
feat2 = Spatial_feat(LF,0.2,10);
feat = cat(2,feat0',feat1',feat2);