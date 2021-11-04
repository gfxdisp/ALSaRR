function [Es,Et,freq_bands] = model_Dst(ref_img, mean_lum, refresh_rate, velocity, persistence, disp_ppd, shading_rates_x, shading_rates_y, mipmap_level_ref, beta, ablation_mask)
%model_Dst Calculate the spatial and temporal distortion in a moving
%texture
%   Takes the reference texture and relevant display parameters as input
%   and returns the amount of spatial and temporal distortions. For
%   predictions scaled in JND units, refer to camojab.m
% ---Output---
% Es        : spatial distortion
% Et        : temporal distortion
% freq_bands: DFT of reference and distorted signal (used for training)
% ---Input----
% ref_img           : grayscale 0-1 normalized image
% mean_lum          : mean (adaptation) luminance in cd/m^2
% refresh_rate      : display refresh rate in Hz
% velocity          : image velocity in deg/s
% persistence       : display persistence (0 to 1)
% disp_ppd          : display pixels-per-degree (assumes square pixel)
% shading_rates_x   : VRS shading rate in x-direction
% shading_rates_y   : VRS shading rate in y-direction [optional]
% mipmap_level_ref  : texture's mipmap level at 1x1 shading rate [optional]
% beta              : power parameter of the model [optional]
% ablation_mask     : on/off mask for VRS distortions, eye motion blur and
%                     temporal distortion [optional]


% Cache the spatio-chromatic CSF class instance
persistent sccsf;

% Which function to approximate sinc function
USE_LORENTZIAN = true;
USE_GAUSSIAN = false;

% load trained beta if not provided
if ~exist('beta', 'var') || isempty(beta)   
    fits = load('camojab_fits.mat');
    beta = fits.beta;    
end

% assume finest mipmap level if not provided
if ~exist('mipmap_level_ref', 'var') || isempty(mipmap_level_ref)    
    mipmap_level_ref = 2;
end

% assume uniform shading rate if not provided
if ~exist('shading_rates_y', 'var') || isempty(shading_rates_y)    
    shading_rates_y = shading_rates_x; % Assume uniform scaling
end

% assumes all components are enabled if not provided
if ~exist('ablation_mask', 'var') || isempty(ablation_mask)    
    ablation_mask = [1 1 1];    
end

% Create an instance of spatio-chromatic CSF
if isempty(sccsf)  
    sccsf = SCCSF_ConeContrastMat();    
end

% Calculate Frequency Spectrum of the texture
F_img = fftshift(abs(fft2(ref_img,512,512)));

% Determine spatial frequency of each entry of the spectrum
[rho_x, rho_y] = create_rho_2D( [512 512], disp_ppd/2, disp_ppd/2 );
rho_x = fftshift(rho_x);
rho_y = fftshift(rho_y);
rho = max(abs(rho_x),abs(rho_y));%sqrt(rho_x.^2 + rho_y.^2);%

% Calculate CSF filter
L_bkg = ones(size(rho)) * mean_lum;
S = stcsf_cam_daly(rho, 0, L_bkg, sccsf);
S(isnan(S)) = 1;

% Initialize Es and Et. The function can take an array of shading rates to
% speed up computation
Es = zeros(size(shading_rates_x));
Et = zeros(size(shading_rates_x));

% build hold-type blur and eye tracking blur filter
b_d = persistence * velocity / refresh_rate;
b_e = persistence * (0.001648*velocity + 0.079818); % coeffiencets from [Denes 2020 pred]
if USE_LORENTZIAN
    B_eye = abs(lorentzian(b_d,rho_x).*lorentzian(b_e,rho_x));
elseif USE_GAUSSIAN
    sigma = sqrt(b_e.^2 + b_d.^2)/pi;
    B_eye = exp(-2.*(pi.*sigma.*rho_x).^2);
else
    B_eye = abs(sinc(b_d*rho_x).*sinc(b_e*rho_x));
end

ref_spectrum = F_img; % reference spectrum
freq_bands = cell([2 length(shading_rates_x)]); % distorted spectrum
% for each shading rate
for ii = 1:length(shading_rates_x)
    rx = shading_rates_x(ii)*disp_ppd/2; % effective x-resolution
    ry = shading_rates_y(ii)*disp_ppd/2; % effective y-resolution
    F_img = ref_spectrum;
%     freq_bands{1,ii} = ref_spectrum  .* S; % used for optimisation
    % build mipmap filter
    b_mip = (2^mipmap_level_ref)/(disp_ppd * min([shading_rates_x(ii), shading_rates_y(ii)]));
    if USE_LORENTZIAN
        B_mip = abs(lorentzian(b_mip,rho_x).*lorentzian(b_mip,rho_y));
    elseif USE_GAUSSIAN
        sigma = b_mip/pi;
        B_mip = exp(-2.*(pi.*sigma.*rho).^2);
    else
        B_mip = abs(sinc(b_mip*rho_x).*sinc(b_mip*rho_y));
    end
    
    % Apply mipmap/VRS distortions
    if ablation_mask(1) == 1
        % Mipmap Filtering
        F_img = F_img .* B_mip;  
        % Mipmap Downsampling
        F_img = convDiracComb(F_img,rx,ry,rho_x,rho_y);
        F_img(abs(rho_x) > rx | abs(rho_y) > ry) = 0; 
    end
    
    % Calculate temporal distortions
    [Et(ii),tfreqb] = model_Et(F_img, beta,refresh_rate,mean_lum,velocity,rho_x,rho, ablation_mask);
    freq_bands{2,ii} = tfreqb; % used for optimisation
    
    % Apply eye motion blur distortion
    if ablation_mask(2) == 1
        % Hold-type and eye blurring
        F_img = F_img .* B_eye; 
    end
    freq_bands{1,ii} = F_img  .* S; % used for optimisation
    
    % Calculate disotrted energy w.r.t to reference
    % CSF normalized visual energy
    Es(ii) = (sum((ref_spectrum(:).*S(:)) .^ beta)/numel(ref_spectrum)) - (sum((F_img(:) .* S(:)) .^ beta)/numel(F_img));
%     Es(ii) = -(sum((F_img(:) .* S(:)) .^ beta)/numel(F_img));
%     Es(ii) = (sum((abs(ref_spectrum(:)-F_img(:)).*S(:)) .^ beta)/numel(ref_spectrum));% - (sum((F_img(:) .* S(:)) .^ beta)/numel(F_img));
end

% Convolution with a Dirac comb to simulate downsampling
function fs = convDiracComb(F,rsx,rsy,rho_x,rho_y)
    xstep = rho_x(1,2) - rho_x(1,1);
    ystep = rho_y(2,1) - rho_y(1,1);
    %Only consider first replica in dirac comb(Accuracy vs Smoothness)
    dcx = zeros([1 round(rsx/xstep*2+1)]);
    dcy = zeros([round(rsy/ystep*2+1) 1]);
    dcx([1, round(length(dcx)/2)]) = 1; 
    dcy([1, round(length(dcy)/2)]) = 1; 
    
%     dcx = dcx ./ sum(dcx(:)); % Normalize energy
%     dcy = dcy ./ sum(dcy(:));
    fs = conv2(F,dcx,'same');
    fs = conv2(fs,dcy,'same');
end

% Lorenztian function (to approximate sinc())
function L = lorentzian(b,rho)
    wl = pi/(2*b);
    L = (wl^2./(wl^2 + 4*rho.^2));
    L(isnan(L)) = 1;
end

% Calculate energy of temporal distortion
function [et,tfreqb] = model_Et(F_img, beta,refresh_rate,mean_lum,velocity,rho_x,rho, ablation_mask)
    if ablation_mask(3) == 1
        % Calculate the first signal alias due to temporal sampling
        F_img_alias = zeros(size(F_img));
        stepsize = abs(rho_x(1,2)-rho_x(1,1));
        velocity = max([velocity,0.1]);
        shift = floor(refresh_rate/2/velocity/stepsize);
        F_img_alias(:,shift+1:end) = F_img(:,1:end-shift);
    else
        et = 0;
        tfreqb = zeros(size(F_img));
        return
    end
    % Calculate energy CSF-normalized energy    
    S_t = stcsf_cam_daly(rho, ones(size(rho))*refresh_rate, ones(size(rho))*mean_lum, sccsf);
    S_t(isnan(S_t)) = 1;    
    F_img_alias = F_img_alias .*S_t;
    tfreqb = F_img_alias;
    et = (sum(F_img_alias(:) .^ beta)/numel(F_img_alias));
end

end