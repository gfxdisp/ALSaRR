function [Q,Dt,Ds] = camojab(texture, display_test_params)
%CAMOJAB Summary of this function goes here
%   Detailed explanation goes here
% texture : [0-1] normalized grayscale diffuse texture image
% display_test_params : struct with following fields:
% display_test_params.mean_lum : mean (adaptation) luminance
% display_test_params.fps : display refresh rate
% display_test_params.velocity : image velocity
% display_test_params.persistence : display persistence
% display_test_params.ppd : display pixels-per-degree
% display_test_params.srx : VRS shading rate in x-direction
% display_test_params.sry : VRS shading rate in y-direction
% display_test_params.mipmap_level : mipmap level of the texture
    params = load('camojab_fits.mat');        
    [Ds,Dt,~] = model_Dst(texture, display_test_params.mean_lum, display_test_params.fps, display_test_params.velocity, display_test_params.persistence, display_test_params.ppd, display_test_params.srx, display_test_params.sry, display_test_params.mipmap_level, params.beta);
    Q = Ds*params.ws + Dt*params.wt;
end

