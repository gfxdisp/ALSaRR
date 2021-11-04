% This script gives an example on how to calculate optimal VRS map given a
% budget and data from g-buffer

% Some sample data has been extracted from a Unity scene of N'th frame
% for this demo.
% sample_data/frame_n        : N-1 frame to approximate luminance
% sample_data/materials_n    : (material_id, mipmap_level) of nth frame
% sample_data/motion_n       : (x velocity, y velocity) in pixels/frame of nth
%                              frame
% sample_data/all_polynomials: All fitted polynomial for this scene
%                              (extract all_polynomials.zip)

if ~exist('exrread','file')
    addpath('../utils/MatlabEXR/')
end

color_image = exrread('sample_data/frame_n.exr');
motion_image = exrread('sample_data/motion_n.exr');
material_image = exrread('sample_data/materials_n.exr');
coeffs_all = loadPolys('sample_data/all_polynomials');

% display configuration
DISPLAY_PEAK_LUMINANCE = 300; %cd/m2
DISPLAY_REFRESH_RATE = 60; %Hz
DISPLAY_PERSISTENCE = 1;
DISPLAY_RESOLUTION = [1920, 1080];
PPD = 60;

% Prepare data for knapsack
% Correct the data units
color_image = color_image * DISPLAY_PEAK_LUMINANCE; % to absolute luminance
motion_image(:,:,1) = motion_image(:,:,1) * DISPLAY_RESOLUTION(1) * DISPLAY_REFRESH_RATE / PPD / 100; % to degrees/sec
motion_image(:,:,2) = motion_image(:,:,2) * DISPLAY_RESOLUTION(2) * DISPLAY_REFRESH_RATE / PPD / 100; % to degrees/sec
motion_image(:,:,3) = 0;
[imageHeight,imageWidth,~] = size(color_image);
% Calculate mean luminance
luminance_image = color_image(:,:,1) * 0.2126 + color_image(:,:,2) * 0.7152 + color_image(:,:,3) * 0.0722;
luminance_image = imresize(luminance_image,1/16,'Method','bilinear');
numBlocksX = floor(imageWidth/16);
numBlocksY = floor(imageHeight/16);
% weights and values for knapsack
N = numBlocksX * numBlocksY;
sr = [[4,4]; [2,4]; [4,2]; [2,2]; [1,2]; [2,1]; [1,1]];
b = zeros([N,length(sr)]);
q = zeros([N,length(sr)]);

% For every tile
nn = 1;
for j = 0:(numBlocksX-1)
    for i = 0:(numBlocksY-1)
        mean_lum = luminance_image(i+1,j+1) + 0.01;
        min_vel_x = motion_image(i+1,j+1,1);
        min_vel_y = motion_image(i+1,j+1,2);
        material_id = material_image(i+1,j+1,1) + 1;
        mip_id = clamp(material_image(i+1,j+1,2)+1,1,12);
        coeff = coeffs_all{material_id};
        
        % For every shading rate of this tile
        for k = 1:length(sr)
            % Calculate quality
            quality_half = @(x,y,z) (coeff.coeff_half(mip_id,1) + coeff.coeff_half(mip_id,2)*x + coeff.coeff_half(mip_id,3)*y + coeff.coeff_half(mip_id,4)*x^2 + coeff.coeff_half(mip_id,5)*x*y + coeff.coeff_half(mip_id,6)*y^2)*(coeff.coeff_half(mip_id,7)*z^3 + coeff.coeff_half(mip_id,8)*z^2 + coeff.coeff_half(mip_id,9)*z + coeff.coeff_half(mip_id,10));
            quality_full = @(x,y,z) (coeff.coeff_full(mip_id,1) + coeff.coeff_full(mip_id,2)*x + coeff.coeff_full(mip_id,3)*y + coeff.coeff_full(mip_id,4)*x^2 + coeff.coeff_full(mip_id,5)*x*y + coeff.coeff_full(mip_id,6)*y^2)*(coeff.coeff_full(mip_id,7)*z^3 + coeff.coeff_full(mip_id,8)*z^2 + coeff.coeff_full(mip_id,9)*z + coeff.coeff_full(mip_id,10));
            if sr(k,1) == 4
                qx = 0;
            elseif sr(k,1) == 2               
                qx = quality_half(min_vel_x,DISPLAY_REFRESH_RATE,log10(mean_lum));
            else
                qx = quality_full(min_vel_x,DISPLAY_REFRESH_RATE,log10(mean_lum));
            end
            if sr(k,2) == 4
                qy = 0;
            elseif sr(k,2) == 2               
                qy = quality_half(min_vel_y,DISPLAY_REFRESH_RATE,log10(mean_lum));
            else
                qy = quality_full(min_vel_y,DISPLAY_REFRESH_RATE,log10(mean_lum));
            end
            q(nn,k) = mean([qx qy]) + 0.1; % adding offset because greedy method doesn't like 0 ratios
            b(nn,k) = 16/sr(k,1)/sr(k,2);
        end        
        nn = nn +1;
    end
end

% Run approximate Greedy and optimal Dynamic Programming solution
gtic = tic;
[vrs_img_g,quality_greedy] = vrs_greedy(.25*N*16, b, q, N, numBlocksX, numBlocksY); % budget = 25% 
time_greedy = toc(gtic);
fprintf("Time (greedy) = %f sec\n",time_greedy);
dtic = tic;
[vrs_img_dp,quality_dp] = vrs_dp(.25*N*16, b, q, N, numBlocksX, numBlocksY);  % budget = 25% 
time_dp = toc(dtic);
fprintf("Time (dp) = %f sec (grows with budget)\n",time_dp);

imshow([color_image/DISPLAY_PEAK_LUMINANCE, imresize((motion_image./(max(motion_image(:)))).^(1/2.2),size(color_image,1:2));imresize(vrs_img_g,size(color_image,1:2)), imresize(vrs_img_dp,size(color_image,1:2))])

function coeffs_all = loadPolys(dir_name)
    files = dir(dir_name);
    coeffs_all = {};
    ii = 1;
    for jj = 3:length(files)
        coeff = load([dir_name,'/',files(jj).name]);
        coeffs_all{ii} = coeff.poly;
        ii = ii + 1;
    end
end
