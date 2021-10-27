texture = rgb2gray(imread2double('wood_texture.png'));

%-----display configuration-----
disp_params.fps = 60; %Hz
disp_params.ppd = 60; %ppd
disp_params.persistence = 1;
disp_params.velocity = 15; %deg/s
disp_params.mean_lum = 150; %cd/m2
disp_params.persistence = 1;
disp_params.mipmap_level = 2;
disp_params.srx = 1;
disp_params.sry = 1;
%--------------------------------

[Q_1x1,~,~] = camojab(texture, disp_params);
disp_params.srx = 0.5;
[Q_2x1,~,~] = camojab(texture, disp_params);
disp_params.srx = 0.25;
[Q_4x1,~,~] = camojab(texture, disp_params);

fprintf("Quality loss for 2x1 shading rate = %f JND\n",Q_2x1 - Q_1x1)
fprintf("Quality loss for 4x1 shading rate = %f JND\n",Q_4x1 - Q_1x1)