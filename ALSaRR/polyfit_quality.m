% This script provides an example of how to sample camojab for a given 
% texture and fit a polynomial to it

if ~exist('camojab','file')
    addpath('../CaMoJAB/');
end
if ~exist('imread2double','file')
    addpath('../utils/')
end

% Sampling range
refresh_rates = 30:15:120;
velocities = [0:10:50];
log_luminance = [-1 1 2 2.5];
shading_resolutions = [1,0.5,0.25];

texture = rgb2gray(imread2double('../CaMoJAB/sample_texture.png'));
if(size(texture,1)>512)
    % crop the image (fft2 is slow for high resolutions)
    texture = texture(1:512,1:512,:); 
end

% We are only interested in quarter, half and full shading rates
delta_Q_quarter_res = zeros([length(refresh_rates) length(velocities) length(log_luminance)]);
delta_Q_half_res = zeros([length(refresh_rates) length(velocities) length(log_luminance)]);
delta_Q_full_res = zeros([length(refresh_rates) length(velocities) length(log_luminance)]);

% setup display configuration
disp_params.ppd = ppd;
disp_params.persistence = 1;
disp_params.mipmap_level = 0; 

% sample camojab
for rr = 1:length(refresh_rates)
    for vv = 1:length(velocities)
        for ll = 1:length(log_luminance)
            q = zeros([1 length(shading_resolutions)]);
            for ss = 1:length(shading_resolutions)
                disp_params.velocity = velocities(vv);
                disp_params.mean_lum = 10^log_luminance(ll);
                disp_params.fps = refresh_rates(rr);
                disp_params.srx = shading_resolutions(ss);
                disp_params.sry = shading_resolutions(ss);
                [q(ss),~,~] = camojab(texture,disp_params);
    %             [q,~] = model_Dst(texture, texture_ppd, disp_params.mean_lum, disp_params.fps, disp_params.velocity, disp_params.persistence, disp_params.ppd, [1 0.5 .25], [1 0.5 .25], disp_params.mipmap_level, params.beta);
    %             q = q*params.ws;
            end
            % set the quarter resolution quality to 0 (origin) and
            % calculate half and full quality w.r.t. to quarter
            delta_Q_quarter_res(rr,vv,ll) = 0;
            delta_Q_half_res(rr,vv,ll) = max(q(2) - q(3),0);
            delta_Q_full_res(rr,vv,ll) = max(q(1) - q(3),0);
        end
    end
end

% take care of floating point precision errors
delta_Q_half_res(abs(delta_Q_half_res)<1e-4)=0;
delta_Q_full_res(abs(delta_Q_full_res)<1e-4)=0;

% Fit a seprate polynomial for both half and full resolution
try 
    [X,Y] = meshgrid(velocities,refresh_rates);
    [xData, yData, zData] = prepareSurfaceData( X, Y, delta_Q_half_res(:,:,length(log_luminance)) );
    [fitresult_half, gof] = fit( [xData, yData], zData, fittype( 'poly22' ) );
    [xData, yData] = prepareCurveData( log_luminance,  squeeze(delta_Q_half_res(6,1,:)./delta_Q_half_res(6,1,end))' );
    [fitresult_half_lum, gof] = fit( xData, yData, fittype( 'poly3' ) );
    quality_half = @(x,y,z) (fitresult_half.p00 + fitresult_half.p10*x + fitresult_half.p01*y + fitresult_half.p20*x^2 + fitresult_half.p11*x*y + fitresult_half.p02*y^2)*(fitresult_half_lum.p1*z^3 + fitresult_half_lum.p2*z^2 + fitresult_half_lum.p3*z + fitresult_half_lum.p4);
    coeff_half = [fitresult_half.p00, fitresult_half.p10, fitresult_half.p01 , fitresult_half.p20 , fitresult_half.p11 , fitresult_half.p02 ,fitresult_half_lum.p1 , fitresult_half_lum.p2 , fitresult_half_lum.p3 ,fitresult_half_lum.p4];

    [xData, yData, zData] = prepareSurfaceData( X, Y, delta_Q_full_res(:,:,length(log_luminance)) );
    [fitresult_full, gof] = fit( [xData, yData], zData, fittype( 'poly22' ) );
    [xData, yData] = prepareCurveData( log_luminance,  squeeze(delta_Q_full_res(6,1,:)./delta_Q_full_res(6,1,end))' );
    [fitresult_full_lum, gof] = fit( xData, yData, fittype( 'poly3' ) );
    quality_full = @(x,y,z) (fitresult_full.p00 + fitresult_full.p10*x + fitresult_full.p01*y + fitresult_full.p20*x^2 + fitresult_full.p11*x*y + fitresult_full.p02*y^2)*(fitresult_full_lum.p1*z^3 + fitresult_full_lum.p2*z^2 + fitresult_full_lum.p3*z + fitresult_full_lum.p4);
    coeff_full = [fitresult_full.p00, fitresult_full.p10, fitresult_full.p01 , fitresult_full.p20 , fitresult_full.p11 , fitresult_full.p02 ,fitresult_full_lum.p1 , fitresult_full_lum.p2 , fitresult_full_lum.p3 ,fitresult_full_lum.p4];
    
    clf
    %visualize
    figure(1)
    for ll = 1:length(log_luminance)
        subplot(2,3,ll);
        [X,Y] = meshgrid(velocities,refresh_rates);
        surf(X,Y,delta_Q_half_res(:,:,ll),'EdgeColor','none','FaceColor','interp'); hold on;
        surf(X,Y, arrayfun(quality_half,X,Y,ones(size(X))*log_luminance(ll)),'EdgeColor','none','FaceAlpha',0.5,'FaceColor','red');hold off;
        xlabel('Velocity (deg/s)')
        ylabel('Refresh Rate (Hz)')
        zlabel('\DeltaQ_{0.5}')
        title( 10^log_luminance(ll) + "cd/m^2")
    %     zlim([0 6])
    end
    %saveas(gcf,"./polyfits_figs/"+name+"_half.fig")

    figure(2)
    for ll = 1:length(log_luminance)
        subplot(2,3,ll);
        [X,Y] = meshgrid(velocities,refresh_rates);
        surf(X,Y,delta_Q_full_res(:,:,ll),'EdgeColor','none','FaceColor','interp'); hold on;
        surf(X,Y, arrayfun(quality_full,X,Y,ones(size(X))*log_luminance(ll)),'EdgeColor','none','FaceAlpha',0.5,'FaceColor','red');hold off;
        xlabel('Velocity (deg/s)')
        ylabel('Refresh Rate (Hz)')
        zlabel('\DeltaQ_{1}')
        title( 10^log_luminance(ll) + "cd/m^2")
    %     zlim([0 6])
    end
    %saveas(gcf,"./polyfits_figs/"+name+"_full.fig")
catch
    coeff_half = zeros([1,10]);
    coeff_full = zeros([1,10]);
    disp("Warning: Fitting error. Outputting 0 poly : " + name);
end

disp('Polynomial coefficients for half rate shading and mipmap level = 0 w.r.t. to quarter rate:')
disp(coeff_half)
disp('Polynomial coefficients for full rate shading and mipmap level = 0 w.r.t. to quarter rate:')
disp(coeff_full)