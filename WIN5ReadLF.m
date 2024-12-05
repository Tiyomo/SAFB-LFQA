function LF = WIN5ReadLF(LF_path)
LF = imread(LF_path);
if size(LF)==[3906, 5625, 3]
   LF = permute(reshape(LF,[[9, 434, 9, 625, 3]]),[1,3,2,4,5]);
else
    LF = permute(reshape(LF,[[9, 512, 9, 512, 3]]),[1,3,2,4,5]);
end
