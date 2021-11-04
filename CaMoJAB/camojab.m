function [Q,Dt,Ds] = camojab(texture, display_test_params)
%CAMOJAB Content-aware model of judder, aliasing and blur
%   The function takes as input the observed texture and the display
%   configuration. It returns the total perceived quality in JND and also
%   the individual spatial and temporal distortions.
% texture                           : [0-1] normalized grayscale diffuse texture image
% display_test_params               : struct with following fields:
% display_test_params.mean_lum      : mean (adaptation) luminance  in cd/m^2
% display_test_params.fps           : display refresh rate in Hz
% display_test_params.velocity      : image velocity in deg/s
% display_test_params.persistence   : display persistence (0 to 1)
% display_test_params.ppd           : display pixels-per-degree (assumes square pixel)
% display_test_params.srx           : VRS shading rate in x-direction
% display_test_params.sry           : VRS shading rate in y-direction [optional]
% display_test_params.mipmap_level  : mipmap level of the texture at 1x1
%                                     shading rate  [optional]
    
    % load pre-fitted model paramters
    params = load('camojab_fits.mat');       
    % Calculate spatial and temporal distorion energy
    [Ds,Dt,~] = model_Dst(texture, display_test_params.mean_lum, display_test_params.fps, display_test_params.velocity, display_test_params.persistence, display_test_params.ppd, display_test_params.srx, display_test_params.sry, display_test_params.mipmap_level, params.beta);
    % Calculate perceived quality in JND units
    Q = Ds*params.ws + Dt*params.wt;
end

