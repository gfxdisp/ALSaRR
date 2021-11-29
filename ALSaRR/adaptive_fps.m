% This script provides an example of how to use CaMoJAB to generate a
% refresh rate vs velocity LUT for a given bandwidth and a scene snapshot


if ~exist('camojab','file')
    addpath('../CaMoJAB/');
end
if ~exist('imread2double','file')
    addpath('../utils/')
end

% Input
max_resolution = [1920 1080];
max_rr = 120; % Hz
min_rr = 30; % Hz
ppd = 45;
bandwidths_percent = [0.0625 .25];% Bandwidths of interest
bandwidths = bandwidths_percent * max_resolution(1) * max_resolution(2) * max_rr;
mean_lum = 100;
texture_path = './sample_data/scene_snapshot.png';
texture = rgb2gray(imread2double(texture_path).^(1));
% Adjust sampling rate of LUT
rr = max_rr:-5:min_rr;
velocities = [linspace(0,50,5),50];

% build display configuration
disp_test_params.ppd = ppd;
disp_test_params.mean_lum = mean_lum;
disp_test_params.persistence = 1;
disp_test_params.mipmap_level = 0;
disp_test_params.srx = 1;
disp_test_params.sry = 1;

[B,V] = meshgrid(bandwidths,velocities);

LUT = zeros([size(B) length(rr)]);
LUT_Max = zeros(size(B));

ppds = zeros(size(rr));
% For each bandwidth
for iB = 1:size(B,1)-1
    % For each velocity
    for iV = 1:size(V,2)
        disp_test_params.velocity = V(iB,iV);
        % For each possible refresh rate
        for iR = 1:length(rr)
           % calculate resolution from refresh rate
           resolution_reduction = sqrt((max_resolution(1)/max_resolution(2)) * B(iB,iV) / rr(iR)) / max_resolution(1);
           if resolution_reduction > 1.001
               LUT(iB,iV,iR) = -100; % resolution cannot be larger than maximum resolution
           else
               disp_test_params.fps = rr(iR);
               disp_test_params.srx = resolution_reduction; % assume uniform shading rate
               disp_test_params.sry = resolution_reduction;
               LUT(iB,iV,iR) =  camojab(texture, disp_test_params); % calculate quality for refresh rate-resolution pair
           end
        end
        tmp = squeeze(LUT(iB,iV,:)); % for this bandwidth and velocity
        tmp = tmp - min(tmp); % shift origin to lowest quality
        tmp(abs(tmp)<1e-2)=0; % to avoid precision errors
        [~,max_q_iR] = max(tmp); % Pick the refresh rate that maximizes quality
        LUT_Max(iB,iV) = rr(max_q_iR);
    end
end
LUT_Max(end,:) = [];
velocities(end) = [];

% Plot LUT
plot(velocities(1:end),LUT_Max(1:end,1),'DisplayName',""+bandwidths_percent(1)*100+"%");hold on;
plot(velocities(1:end),LUT_Max(1:end,2),'DisplayName',""+bandwidths_percent(2)*100+"%");hold on;
xlabel('Velocity (deg/s)')
ylabel('Refresh rate (Hz)')
legend;
hold off
