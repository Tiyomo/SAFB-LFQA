function feat2 = Spatial_feat(LF,step,k)

% Version 1.0: For Zerui Yu "Spatial-Angular Features Based No-Reference Light
% Field Quality Assessment" submitted to ESWA.

% inputs:
% LF - the light field images with size [u, v, s, t, Ch]
% step - determine the number of refocused images in the focus stack, default: 0.2
% k - determine the number of bins of joint statistics, default: 10

% outputs:
% feat2 - the spatial feats of distorted LF, size: (1, 4*k)

slope=-1:step:1;
[u, v, s, t, Ch] = size(LF);
imstack = zeros(s,t,Ch,length(slope));
stack2 = zeros(length(slope),4*k);
for j=1:length(slope)
    img = LFFiltShiftSum(LF,slope(j));
    imstack(:,:,:,j) = double(img(:,:,1:3));
    img = rgb2gray(img(:,:,1:3)); 
    img = (img-min(img(:)))/(max(img(:))-min(img(:)))*255;
    img = double(img);
    out = Joint_statistics_Phase_MSCN(img,k);
    stack2(j,:) = out;
end
feat2 = mean(stack2);
end

function out =  Joint_statistics_Phase_MSCN(imd,k)

% inputs:
% imd - the distorted refocus image of LF (grayscale image, double type, 0~255)
% k -  the length of the marginal probability distribution vector (default: k=10)

% outputs:
% out(1:20): the marginal distributions for phasecong and MSCN coefficients. 
% out(21:40): the probability distributions between phasecong and MSCN coefficients. 

sigma = 0.5;
phase_im = phasecong2(imd,4,6);
phase_im = (phase_im-min(phase_im(:)))/(max(phase_im(:))-min(phase_im(:)))*255;

window2 = fspecial('gaussian',7,7/6);
window2 = window2/sum(window2(:));
mu            = filter2(window2, imd, 'same');
mu_sq         = mu.*mu;
sigma0         = sqrt(abs(filter2(window2, imd.*imd, 'same') - mu_sq));
mscn_im     = (imd-mu)./(sigma0+1);
mscn_im = (mscn_im-min(mscn_im(:)))/(max(mscn_im(:))-min(mscn_im(:)))*255;

%Normalization
c0 = 4*0.05;
sigmaN = 2*sigma;
window1 = fspecial('gaussian',2*ceil(3*sigmaN)+1, sigmaN);
window1 = window1/sum(window1(:));
Nmap = sqrt(filter2(window1,mean(cat(3,phase_im,mscn_im).^2,3),'same'))+c0;
phase_im = (phase_im)./Nmap;
mscn_im = (mscn_im)./Nmap;
% remove the borders, which may be the wrong results of a convolution
% operation
h = ceil(3*sigmaN);
phase_im = abs(phase_im(h:end-h+1,h:end-h+1,:));
mscn_im = abs(mscn_im(h:end-h+1,h:end-h+1));

ctrs{1} = 1:k;ctrs{2} = 1:k;
% histogram computation
step1 = 0.20;
step2 = 0.20;
phase_qun = ceil(phase_im/step1);
MSCN_qun = ceil(mscn_im/step2);

N1 = hist3([phase_qun(:),MSCN_qun(:)],ctrs);
N1 = N1/sum(N1(:));
NP = sum(N1,2); NM = sum(N1,1);

alpha1 = 0.0001;
cp_PM = N1./(repmat(NM,size(N1,1),1)+alpha1);
cp_PM_H=  sum(cp_PM,2)';
cp_PM_H = cp_PM_H/sum(cp_PM_H);
cp_MP = N1./(repmat(NP,1,size(N1,2))+alpha1);
cp_MP_H = sum(cp_MP,1);
cp_MP_H = cp_MP_H/(sum(cp_MP_H));

out = [NP', NM, cp_PM_H, cp_MP_H];
end