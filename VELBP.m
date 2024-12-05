%  J = LBPV(I,R,N,MAPPING) returns the feat of local binary pattern variance
%  of an intensity image I. The LBP codes are computed using P sampling 
%  points on a circle of radius R and using mapping table defined by MAPPING. 
%  The VAR values are computed for the P sampling points. Instead of
%  computing the joint histogram of LBP and VAR globally, the LBPV computes the 
%  VAR from a local region and accumulates it into the LBP bin. This can be regarded 
%  as the integral projection along the VAR coordinate. 
%  Thanks for the origin code from Marko Heikkil?, Timo Ahonen, and Zhenhua Guo.

%  Examples
%  --------
%       I=imread('rice.png');
%       mapp=load('mapping1.mat'); 
%       mapping=mapp.mapping
%       lbp_feat=LBPV(I,1,8,mapping);


function result = VELBP(I,R,P,MAPPING)
% Version 1.0
% Authors: Zhenhua Guo, Lei Zhang and David Zhang
% Copyright @ Biometrics Research Centre, the Hong Kong Polytechnic University

if nargin<1
    disp('No input image')
    return
end
if nargin<2
    R = 1;
end
if nargin<3
    P = 8;
end   
if nargin<4
    MAPPING = getmapping(P,'riu2'); 
end

% Get LBP value for each pixel of the input image
LBPMap = multi_threshold_lbp(I,R,P,MAPPING,'x',R/2);  

% Get VAR value for each pixel of the input image
spoints=zeros(P,2);

% Angle step.
a = 2*pi/P;

for i = 1:P
    spoints(i,1) = -R*sin((i-1)*a);
    spoints(i,2) = R*cos((i-1)*a);
end

[ysize xsize] = size(I);

miny=min(spoints(:,1));
maxy=max(spoints(:,1));
minx=min(spoints(:,2));
maxx=max(spoints(:,2));

% Block size, each VAR value is computed within a block of size bsizey*bsizex
bsizey=ceil(max(maxy,0))-floor(min(miny,0))+1;
bsizex=ceil(max(maxx,0))-floor(min(minx,0))+1;

% Coordinates of origin (0,0) in the block
origy=1-floor(min(miny,0));
origx=1-floor(min(minx,0));

% Minimum allowed size for the input image depends
% on the radius of the used VAR operator.
if(xsize < bsizex || ysize < bsizey)
  error('Too small input image. Should be at least (2*radius+1) x (2*radius+1)');
end

% Calculate dx and dy;
dx = xsize - bsizex;
dy = ysize - bsizey;

% convert the input image to "double" type.
d_image = double(I);

% Reorganize the image to compute VAR value
for i = 1:P
  y = spoints(i,1)+origy;
  x = spoints(i,2)+origx;
  % Calculate floors, ceils and rounds for the x and y.
  fy = floor(y); cy = ceil(y); ry = round(y);
  fx = floor(x); cx = ceil(x); rx = round(x);
  % Check if interpolation is needed.
  if (abs(x - rx) < 1e-6) && (abs(y - ry) < 1e-6)
    % Interpolation is not needed, use original datatypes
    T = d_image(ry:ry+dy,rx:rx+dx);
    
    % convert the matrix to a column vector
    NArray(:,i) = reshape(T,prod(size(T)),1); 
  else
    % Interpolation needed, use double type images 
    ty = y - fy;
    tx = x - fx;

    % Calculate the interpolation weights.
    w1 = roundn((1 - tx) * (1 - ty),-6);
    w2 = roundn(tx * (1 - ty),-6);
    w3 = roundn((1 - tx) * ty,-6) ;
    % w4 = roundn(tx * ty,-6) ;
    w4 = roundn(1 - w1 - w2 - w3, -6);

    % Compute interpolated pixel values
    T = w1*d_image(fy:fy+dy,fx:fx+dx) + w2*d_image(fy:fy+dy,cx:cx+dx) + ...
        w3*d_image(cy:cy+dy,fx:fx+dx) + w4*d_image(cy:cy+dy,cx:cx+dx);
    T = roundn(T,-4);
    % convert the matrix to a column vector
    NArray(:,i) = reshape(T,prod(size(T)),1);
  end  
end

% after get all neighborhood for each pixel, compute variance for the image
VAR = var(NArray,1,2);
VAR = reshape(VAR,size(T));

MapNum = max(MAPPING.table(:));

% Initialize the result matrix with zeros.
result=zeros(1,MapNum+1);
for k=0:MapNum;
    index = find(LBPMap==k);
    result(k+1) = sum(VAR(index));
end
result = result/(sqrt(sum(result.^2)));
end

function [result,wo_mapping_img,mapping_img] = multi_threshold_lbp(varargin) % image,radius,neighbors,mapping,mode)
% Version 0.3.3
% Authors: Marko Heikkil?and Timo Ahonen

% Changelog
% Version 0.3.2: A bug fix to enable using mappings together with a
% predefined spoints array
% Version 0.3.1: Changed MAPPING input to be a struct containing the mapping
% table and the number of bins to make the function run faster with high number
% of sampling points. Lauge Sorensen is acknowledged for spotting this problem.


% Check number of input arguments.
narginchk(1,6);

image=varargin{1};
d_image=double(image);

if nargin==1
    spoints=[-1 -1; -1 0; -1 1; 0 -1; -0 1; 1 -1; 1 0; 1 1];
    neighbors=8;
    mapping=0;
    mode='h';
end

if (nargin == 2) && (length(varargin{2}) == 1)
    error('Input arguments');
end

if (nargin > 2) && (length(varargin{2}) == 1)
    radius=varargin{2};
    neighbors=varargin{3};
    
    spoints=zeros(neighbors,2);
    threshold = varargin{6};
    % Angle step.
    a = 2*pi/neighbors;
    
    for i = 1:neighbors
        spoints(i,1) = -radius*sin((i-1)*a);
        spoints(i,2) = radius*cos((i-1)*a);
    end
    
    if(nargin >= 4)
        mapping=varargin{4};
        if(isstruct(mapping) && mapping.samples ~= neighbors)
            error('Incompatible mapping');
        end
    else
        mapping=0;
    end
    
    if(nargin >= 5)
        mode=varargin{5};
    else
        mode='h';
    end
end

if (nargin > 1) && (length(varargin{2}) > 1)
    spoints=varargin{2};
    neighbors=size(spoints,1);
    
    if(nargin >= 3)
        mapping=varargin{3};
        if(isstruct(mapping) && mapping.samples ~= neighbors)
            error('Incompatible mapping');
        end
    else
        mapping=0;
    end
    
    if(nargin >= 4)
        mode=varargin{4};
    else
        mode='h';
    end   
end

% Determine the dimensions of the input image.
[ysize xsize] = size(image);



miny=min(spoints(:,1));
maxy=max(spoints(:,1));
minx=min(spoints(:,2));
maxx=max(spoints(:,2));

% Block size, each LBP code is computed within a block of size bsizey*bsizex
bsizey=ceil(max(maxy,0))-floor(min(miny,0))+1;
bsizex=ceil(max(maxx,0))-floor(min(minx,0))+1;

% Coordinates of origin (0,0) in the block
origy=1-floor(min(miny,0));
origx=1-floor(min(minx,0));

% Minimum allowed size for the input image depends
% on the radius of the used LBP operator.
if(xsize < bsizex || ysize < bsizey)
  error('Too small input image. Should be at least (2*radius+1) x (2*radius+1)');
end

% Calculate dx and dy;
dx = xsize - bsizex;
dy = ysize - bsizey;

% Fill the center pixel matrix C.
C = image(origy:origy+dy,origx:origx+dx);
d_C = double(C) + threshold;

bins = 2^neighbors;

% Initialize the result matrix with zeros.
result=zeros(dy+1,dx+1);

%Compute the LBP code image

for i = 1:neighbors
  y = spoints(i,1)+origy;
  x = spoints(i,2)+origx;
  % Calculate floors, ceils and rounds for the x and y.
  fy = floor(y); cy = ceil(y); ry = round(y);
  fx = floor(x); cx = ceil(x); rx = round(x);
  % Check if interpolation is needed.
  if (abs(x - rx) < 1e-6) && (abs(y - ry) < 1e-6)
    % Interpolation is not needed, use original datatypes
    N = image(ry:ry+dy,rx:rx+dx);
    D = N >= C; 
  else
    % Interpolation needed, use double type images 
    ty = y - fy;
    tx = x - fx;

    % Calculate the interpolation weights.
    w1 = roundn((1 - tx) * (1 - ty),-6);
    w2 = roundn(tx * (1 - ty),-6);
    w3 = roundn((1 - tx) * ty,-6) ;
    % w4 = roundn(tx * ty,-6) ;
    w4 = roundn(1 - w1 - w2 - w3, -6);
            
    % Compute interpolated pixel values
    N = w1*d_image(fy:fy+dy,fx:fx+dx) + w2*d_image(fy:fy+dy,cx:cx+dx) + ...
        w3*d_image(cy:cy+dy,fx:fx+dx) + w4*d_image(cy:cy+dy,cx:cx+dx);
    N = roundn(N,-4);
    D = N >= d_C; 
  end  
  % Update the result matrix.
  v = 2^(i-1);
  result = result + v*D;
end

wo_mapping_img = result;
%Apply mapping if it is defined
if isstruct(mapping)
    bins = mapping.num;
    for i = 1:size(result,1)
        for j = 1:size(result,2)
            result(i,j) = mapping.table(result(i,j)+1);
        end
    end
end
mapping_img = result;

if (strcmp(mode,'h') || strcmp(mode,'hist') || strcmp(mode,'nh') || strcmp(mode,'l2'))
    % Return with LBP histogram if mode equals 'hist'.
    result=hist(result(:),0:(bins-1));
    if (strcmp(mode,'nh'))
        result=result/sum(result);
    end
    if (strcmp(mode,'l2'))
        result = result/(sqrt(sum(result.^2)));
    end
else
    %Otherwise return a matrix of unsigned integers
    if ((bins-1)<=intmax('uint8'))
        result=uint8(result);
    elseif ((bins-1)<=intmax('uint16'))
        result=uint16(result);
    else
        result=uint32(result);
    end
end

end

function x = roundn(x, n)

narginchk(2, 2)
validateattributes(x, {'single', 'double'}, {}, 'ROUNDN', 'X')
validateattributes(n, ...
    {'numeric'}, {'scalar', 'real', 'integer'}, 'ROUNDN', 'N')

if n < 0
    p = 10 ^ -n;
    x = round(p * x) / p;
elseif n > 0
    p = 10 ^ n;
    x = p * round(x / p);
else
    x = round(x);
end

end
