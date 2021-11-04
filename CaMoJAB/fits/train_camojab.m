function [ws,wt, beta] = train_camojab(ablation_mask)
% Fit CaMoJAB to VRS Experiment and MARRR Experiemnt data.
% ws : fitted weight of spatial artifacts
% wt : fitted weight of temporal artifacts
% beta : fitted power parameter of the model
% ablation_mask : on/off mask for different model components [optional]

if ~exist('model_Dst','file')
    addpath('../')
end
if ~exist('ablation_mask', 'var') || isempty(ablation_mask)    
    ablation_mask = [1 1 1];
end

% Load the stimulus textures
imgs = load('../../Experiment_Data/vrs_exp_stimulus.mat');

%-------------------Model predictions for VRS Experiment-------------------

% Experiment configurations
configs = ["Config1","Config2","Config3"];
velocities = [0,10,30;3,10,20;10,45,75];
resolutions = [1,1/2,1/4;1,1/2,1/4;1,1/2,1/4];
stimuli = ["checkerbox","gradient","grass","noise"];

freq_bands_vrs_st = cell([length(stimuli) size(resolutions,2) size(velocities,2) length(configs) 2]);

for cc=1:length(configs)    
    if strcmp(configs(cc), 'Config1')==1    % Block PC
        ppd = 50;
        fps = 144;
        persistence = 1;
        mean_lum = 150;
        mipmap_level_ref = 3; % extracted from Unity
    elseif strcmp(configs(cc), 'Config2')==1   % Block Mobile
        ppd = 90;
        fps = 60;
        persistence = 1;
        mean_lum = 75;
        mipmap_level_ref = 1;
    elseif strcmp(configs(cc), 'Config3')==1   % Block VR
        ppd = 26; % 2x SSAA
        fps = 144;
        persistence = 0.05;
        mean_lum = 2.5;
        mipmap_level_ref = 5;
    end
    srx = [1,1/2,1/4];
    sry = ones(size(srx));
    for ii=1:length(stimuli)
        switch stimuli(ii)
            case "checkerbox"
                ref_img = (imgs.checkerboard);
            case "gradient"
                ref_img = (imgs.gradient_linear);
            case "noise"
                ref_img = (imgs.noise);
            case "grass"
                ref_img = (imgs.grass);
        end
        for vv=1:size(velocities,2)
            vel = abs(velocities(cc,vv) - 2); % -2 to account for persistent rotational motion in the experiment
            for rr = 1:length(srx)
                % Only calculate frequency bands and calculate energy
                % during model fitting to speed up optimisation
                [~,~,fband] = model_Dst(ref_img, mean_lum, fps, vel, persistence, ppd, srx(rr), sry(rr),mipmap_level_ref,1,ablation_mask);
                freq_bands_vrs_st{ii,rr,vv,cc,2} = fband{2}; % temporal distortion band
                freq_bands_vrs_st{ii,rr,vv,cc,1} = fband{1}; % spatial distortion band
            end
        end
    end
end

%---------------------Model Predictions Denes Exp 1------------------------
% Experiment Configuration
fov = 30; %deg
screen_width = 2560; %px
ppd = screen_width / fov; % pix per deg
persistence = 1;
mean_lum = 100;

freq_bands_marrr_st = cell([3 24 2]);
marrr_refresh_rates = 50:5:165;
marrr_velocities = [15,30,45];

% Model predictions
for vv = 1:length(marrr_velocities)
    for rr=1:length(marrr_refresh_rates)     
        [~,~,fband] = model_Dst(imgs.marrr_checkerboard, mean_lum, marrr_refresh_rates(rr), marrr_velocities(vv), persistence, ppd, 1, 1, 1,1,ablation_mask);
        freq_bands_marrr_st{vv,rr,2} = fband{2};
        freq_bands_marrr_st{vv,rr,1} = fband{1};                 
    end
end

%---------------------Model Predictions Denes Exp 3------------------------
data = load('../../Experiment_Data/external_data/marrr_exp3_jnds.mat');
% Experiment configuration
refresh_rates = data.f;

fbands_marrr_exp3 = cell(size(refresh_rates));

% Model predictions
for rr = 1:length(refresh_rates)
    [~,~,fband] = model_Dst(imgs.marrr_checkerboard, 100, refresh_rates(rr), 25, 1, ppd, 1, 1, 1,1,ablation_mask);    
    fbands_marrr_exp3{rr} = fband{2};        
end
%----------------------------------------------------
disp('INFO: Fitting data (may take >30mins')
[ws,wt,beta,a1,a2,a3] = fit_to_data_together(freq_bands_vrs_st,freq_bands_marrr_st, fbands_marrr_exp3);

camojab_vrs = freq_bands2energy_vrs(freq_bands_vrs_st, beta, ws, wt) .* a1;
camojab_marrr_exp1 = freq_bands2energy_marrr(freq_bands_marrr_st, beta, ws, wt) .* a2;
camojab_marrr_exp3 = freq_bands2energy_exp3(fbands_marrr_exp3, beta, a3);
fitting_err(camojab_vrs, camojab_marrr_exp1, camojab_marrr_exp3);

disp('INFO: Plotting results (may take >5mins')
plot_training_results(camojab_vrs,camojab_marrr_exp1,ws,wt, beta)
% save 'camojab_fits.mat' 'ws' 'wt' 'beta' 'a1' 'a2' 'a3'

% Vectorized method to turn freq bands to energy for VRS experiment
function E = freq_bands2energy_vrs(freq_bands, beta, ws, wt)
    energyFun = @(fb) sum( fb(:).^beta ) / numel(fb);
    E_s = -cellfun(energyFun,freq_bands(:,:,:,:,1));
    E_t = cellfun(energyFun,freq_bands(:,:,:,:,2));
    E = ws*E_s + wt*E_t;
    E = E - E(:,1,:,:);    
end

% Vectorized method to turn freq bands to energy for MARRR experiment 1
function E = freq_bands2energy_marrr(freq_bands, beta, ws, wt)
    energyFun = @(fb) sum( fb(:).^beta ) / numel(fb);
    E_s = -cellfun(energyFun,freq_bands(:,:,1));
    E_t = cellfun(energyFun,freq_bands(:,:,2));
    E = ws*E_s + wt*E_t;
    E = E - E(:,1);    
end

% Vectorized method to turn freq bands to energy for MARRR experiment 3
function E = freq_bands2energy_exp3(freq_bands, beta, w)
    energyFun = @(fb) sum( fb(:).^beta ) / numel(fb);
    E = cellfun(energyFun,freq_bands);
    E = w*E;
    E = E - E(end);
end

% Linearly fit model parameters to experiment data
function [ws,wt,beta,a1,a2,a3] = fit_to_data_together(freq_bands_vrs_st,freq_bands_marrr_st, fbands_marrr_exp3)
    % Scaled experiment results for VRS experiment
    GT = load('../../Experiment_Data/vrs_exp_data_jnd.mat');
    GT = GT.GT(:,:,:,:);
    % Scaled experiment results for MARRR experiment 1
    marrr_data = load('../../Experiment_Data/external_data/marrr_exp1_jnds.mat');
    marrr_data_exp1 = [marrr_data.res.JODs(1,:); marrr_data.res.JODs(2,:); marrr_data.res.JODs(3,:)];
    marrr_data_exp1 = marrr_data_exp1 - marrr_data_exp1(:,1);
    % Scaled experiment results for MARRR experiment 3
    marrr_data_exp3 = load('../../Experiment_Data/external_data/marrr_exp3_jnds.mat');
    marrr_jnds_exp3 = marrr_data_exp3.Q; marrr_jnds_exp3 = marrr_jnds_exp3(:,1);
    
    fun = @(vars) sqrt(mean([reshape((freq_bands2energy_vrs(freq_bands_vrs_st,vars(3),vars(1),vars(2))*vars(4) - GT) .^ 2,...
        [1 , numel(GT)]), reshape(((freq_bands2energy_marrr(freq_bands_marrr_st,vars(3),vars(1), vars(2))...
        ).*vars(5) - marrr_data_exp1) .^ 2,[1 , numel(marrr_data_exp1), ])...
        ,reshape((freq_bands2energy_exp3(fbands_marrr_exp3, vars(3),vars(6)) - marrr_jnds_exp3).^2,[1,numel(marrr_jnds_exp3)])])); 
    options = optimset('Display','iter','PlotFcns',@optimplotfval);
    [vars,rmse] = fminsearch(fun,[-1.4538, -2.6915, 0.348, 3, 7.8395, 5.4194]);
    ws = vars(1);
    wt = vars(2);
    beta = vars(3);
    a1 = vars(4);
    a2 = vars(5);
    a3 = vars(6);
    disp("rmse = " + rmse);
end

% Calculate RMSE of fitted model
function [rmse_c1,rmse_c2,rmse_c3,rmse_marrr_exp1,rmse_marrr_exp3,total] = fitting_err(camojab_vrs, camojab_marrr_exp1, camojab_marrr_exp3)
    GT = load('../../Experiment_Data/vrs_exp_data_jnd.mat');GT = GT.GT;
    marrr_data = load('../../Experiment_Data/external_data/marrr_exp1_jnds.mat');
    marrr_data_exp1 = [marrr_data.res.JODs(1,:); marrr_data.res.JODs(2,:); marrr_data.res.JODs(3,:)];
    marrr_data_exp1 = marrr_data_exp1 - marrr_data_exp1(:,1);
    marrr_data_exp3 = load('../../Experiment_Data/external_data/marrr_exp3_jnds.mat');
    marrr_jnds_exp3 = marrr_data_exp3.Q; marrr_jnds_exp3 = marrr_jnds_exp3(:,1);

    rmse_c1 = sqrt(mean(reshape((camojab_vrs(:,:,:,1) - GT(:,:,:,1)) .^ 2,[1 , numel(camojab_vrs(:,:,:,1))])));
    rmse_c2 = sqrt(mean(reshape((camojab_vrs(:,:,:,2) - GT(:,:,:,2)) .^ 2,[1 , numel(camojab_vrs(:,:,:,2))])));
    rmse_c3 = sqrt(mean(reshape((camojab_vrs(:,:,:,3) - GT(:,:,:,3)) .^ 2,[1 , numel(camojab_vrs(:,:,:,3))])));
    rmse_marrr_exp1 = sqrt(mean(reshape((camojab_marrr_exp1 - marrr_data_exp1) .^ 2,[1 , numel(marrr_data_exp1)])));
    rmse_marrr_exp3 = sqrt(mean(reshape((camojab_marrr_exp3 - marrr_jnds_exp3) .^ 2,[1 , numel(marrr_jnds_exp3)])));
    total = sqrt(mean([reshape((camojab_vrs - GT) .^ 2,[1 , numel(camojab_vrs)]), reshape((camojab_marrr_exp1 - marrr_data_exp1) .^ 2,[1 , numel(marrr_data_exp1)]), reshape((camojab_marrr_exp3 - marrr_jnds_exp3) .^ 2,[1 , numel(marrr_jnds_exp3)]) ]));
    fprintf("RMSE | C1 = %f | C2 = %f | C3 = %f | MARRR_E1 = %f | MARRR_E3 = %f | Total = %f\n", rmse_c1,rmse_c2,rmse_c3,rmse_marrr_exp1,rmse_marrr_exp3,total);
end

end