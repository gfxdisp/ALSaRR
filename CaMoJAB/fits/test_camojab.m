
if ~exist('camojab','file')
    addpath('../')
end

%------------------------------------------------------------------------------------------------------------------------------
% Fit to Series, B. T. "The present state of ultra-high definition television." (2012).
%Figure 20 data (extracted using WebPlotDigitizer)
disp('INFO: Testing on ITU FPS dataset')
x = [[30,1];[60,1];[60,.5];[120,1];[120,.5];[240,1]]; % [fps,shutter angle]
err_low = [-3.1,-2.2,-1.8,-1.2,-0.92,-0.066];
err_high = [-2.8,-1.9,-1.4,-0.90,-0.66,0.099];
quality_fps = [-2.9,-2.0,-1.6,-1.0,-0.79,0.011];

mask = [3,5]; % We do not model shutter angle
x(mask,:) = [];
err_low(mask) = [];
err_high(mask) = [];
quality_fps(mask) = [];

% Reference content
% Due to lack of original videos, we assume checkerboard
imgs_config1 = load('../../Experiment_Data/vrs_exp_stimulus.mat');
texture = imgs_config1.checkerboard;

% display configuration
disp_test_params.ppd = round(1080 / 2 / rad2deg(atan(1/6)));
disp_test_params.mean_lum = 300;
disp_test_params.persistence = 1;
disp_test_params.velocity = 20;
disp_test_params.srx = 1;
disp_test_params.sry = 1;
disp_test_params.mipmap_level = 1;

% Linear fit and report error
camojab_q = zeros(size(quality_fps));
for ii=1:length(x)       
    disp_test_params.fps = x(ii,1);
    camojab_q(ii) = camojab(texture, disp_test_params); 
end
plcc = corr([camojab_q',quality_fps'],'type','Pearson');plcc = plcc(1,2);
srocc = corr([camojab_q',quality_fps'],'type','Spearman');srocc = srocc(1,2);
disp("PLCC = " + plcc + " | SROCC = " + srocc);
camojab_q = camojab_q - camojab_q(end);
[c1,c2] = fit_camojab_to_itu_data(camojab_q,quality_fps);

% Densly sample for continuous plots
fps = 30:10:240;
camojab_q = zeros(size(fps));
for ff = 1:length(fps)
    disp_test_params.fps = fps(ff);
    camojab_q(ff) = camojab(texture, disp_test_params)*c1+c2; 
end
camojab_q = camojab_q - camojab_q(end) + c2;

% Plot
figure(1)
plot(fps,camojab_q,'--k');hold on
errorbar(x(:,1),quality_fps,quality_fps-err_low,err_high-quality_fps,'-k'); hold off   
ylim([-3.5 0.2])
xlim([25,250])
xticks(x(:,1))
yticks([-3 -2 -1 0])
set(gca,'TickLabelInterpreter','latex','FontSize', 16)
ylabel("Quality Scale",'FontSize', 20,'Interpreter','latex')
xlabel("Refresh Rate (Hz)",'FontSize', 20,'Interpreter','latex')
title(["Average of 5 sequences, 3H viewing distance", "(with $95\%$ Confidence Intervals shown)"],'FontSize', 18,'Interpreter','latex')
legend(["CaMoJAB", "ITU-R BT"],'Location','best','FontSize', 20,'Interpreter','latex')
grid on
% exportgraphics(gcf,"figs/stdm_ITU_fps.pdf",'ContentType','vector')

%------------------------------------------------------------------------------------------------------------------------------
% Fit to Series, B. T. "The present state of ultra-high definition television." (2012).
%Figure 7 data (extracted using WebPlotDigitizer)
disp('INFO: Testing on ITU PPD dataset')
max_ppd = round(2160/ 2 / rad2deg(atan(1/3)));
ppds = [.125, .25, .25, .5, .5, 1] * max_ppd;
x = ["480x270 Lanczos", "960x540 NN", "960x540 Lanczos", "1920x1080 NN", "1920x1080 Lanczos", "3840x2160"];
quality_ppd =  [12.7,39.0,41.8,73.0,79.5,89.0];
err_low =  [11.0,36.3,39.0,70.2,78.0,87.5];
err_high = [13.5,41.3,44.3,75.0,81.0,90.2];
mask = [1,3,5]; % We do not model Lanczos filter

ppds(mask) = [];
x(mask) = [];
quality_ppd(mask) = [];
err_low(mask) = [];
err_high(mask) = [];

% Reference content
% Due to lack of original videos, we assume checkerboard
texture = imgs_config1.checkerboard;

% display configuration
disp_test_params.ppd = max_ppd;
disp_test_params.mean_lum = 300;
disp_test_params.persistence = 1;
disp_test_params.fps = 60;
disp_test_params.velocity = 6;
disp_test_params.mipmap_level = 0;

% Linear fit and report error
camojab_q = zeros(size(quality_ppd));
for ii=1:length(ppds)       
    disp_test_params.srx = ppds(ii)/disp_test_params.ppd;
    camojab_q(ii) = camojab(texture, disp_test_params); 
end
plcc = corr([camojab_q',quality_ppd'],'type','Pearson');plcc = plcc(1,2);
srocc = corr([camojab_q',quality_ppd'],'type','Spearman');srocc = srocc(1,2);
disp("PLCC = " + plcc + " | SROCC = " + srocc);
[c1,c2] = fit_camojab_to_itu_data(camojab_q,quality_ppd);

% Densly sample for continuous plots
cppds = 10:3:65;
camojab_q = zeros(size(ppds));
for pp = 1:length(cppds)    
    disp_test_params.srx = cppds(pp)/max(cppds);
    camojab_q(pp) = camojab(texture,  disp_test_params)*c1+c2; 
end

% Plot
figure(2)
plot(cppds,camojab_q,'LineWidth',1,'Color',[0,0,0],'LineStyle','--','DisplayName','CaMoJAB');hold on
bar(ppds,quality_ppd,'FaceColor',[0.5,0.5,0.5], 'BarWidth', .25, 'FaceAlpha', 0.5, 'DisplayName', 'ITU-R BT');hold on
er = errorbar(ppds,quality_ppd,quality_ppd - err_low,err_high - quality_ppd,'HandleVisibility','off');    
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
ylim([0 100])
xlim([10,65])
% xticks(ppds)
set(gca,'TickLabelInterpreter','latex','FontSize', 16)
set(gca,'xtick',ppds,'xticklabel',{'960x540','1920x1080','3840x2160'})
yticks(0:20:100)
ylabel("Perceived Video Quality",'FontSize', 20,'Interpreter','latex')
xlabel("Resolution(PPD)",'FontSize', 20,'Interpreter','latex')
title(["Impact of image resolution on preceived video quality", "using HD/UHD displays at 1.5H viewing distance"],'FontSize', 18,'Interpreter','latex')
legend('Location','best','FontSize', 20,'Interpreter','latex')
grid on
hold off;
% exportgraphics(gcf,"figs/stdm_ITU_ppd.pdf",'ContentType','vector')

%------------------------------------------------------------------------------------------------------------------------------
% Fit to Mackins et al. [2016] data
% Data extracted using WebPlotDigitizer

disp('INFO: Testing on Mackins dataset')
velocities = [10,30,50,70];
mis = zeros([7,2,4]);
mis_err_high = zeros([7,1,4]);
mis_err_low = zeros([7,1,4]);

mis(:,:,1) = [[60, 3.724137931034483];[100, 4.565517241379311];[150, 4.855172413793103];[300, 4.910344827586207];[600, 4.855172413793103];[1000, 4.958620689655173];[2000, 4.758620689655173]];
mis(:,:,2) = [[60, 2.406896551724138];[100, 3.2068965517241383];[150, 3.8965517241379315];[300, 4.551724137931035];[600, 4.862068965517242];[1000, 4.9655172413793105];[2000, 4.841379310344828]];
mis(:,:,3) = [[60, 2.103448275862069];[100, 2.606896551724138];[150, 3.158620689655173];[300, 4.255172413793104];[600, 4.710344827586207];[1000, 4.606896551724138];[2000, 4.710344827586207]];
mis(:,:,4) = [[60, 1.9034482758620692];[100, 2.344827586206897];[150, 2.462068965517242];[300, 2.710344827586207];[600, 4.662068965517242];[1000, 4.855172413793103];[2000, 4.5103448275862075]];

mis_err_high(:,:,1) = [3.91;4.68;4.93;4.97;4.95;4.99;4.87];
mis_err_high(:,:,2) = [2.60;3.40;4.07;4.68;4.94;4.97;4.93];
mis_err_high(:,:,3) = [2.34;2.80;3.38;4.40;4.79;4.70;4.82];
mis_err_high(:,:,4) = [2.11;2.59;2.64;2.86;4.76;4.92;4.63];

mis_err_low(:,:,1) = [3.49;4.43;4.78;4.83;4.76;4.90;4.64];
mis_err_low(:,:,2) = [2.21;3.01;3.74;4.43;4.78;4.83;4.78];
mis_err_low(:,:,3) = [1.87;2.41;2.92;4.11;4.61;4.50;4.58];
mis_err_low(:,:,4) = [1.69;2.11;2.26;2.54;4.54;4.78;4.38];

figure(3)
for vv = 1:length(velocities)
    plot(mis(:,1,vv),mis(:,2,vv));hold on;
end
hold off;
xlabel('Frame Rate(Hz)')
ylabel('Mean Impairment Score')
legend({'10 deg/s','30 deg/s','50 deg/s','70 deg/s'},'Location','southeast')

camojab_q = mis;
% Reference content (checkerboard is a good approximation of printed lines)
texture = imgs_config1.checkerboard;

% display configuration
disp_test_params.ppd = 60;
disp_test_params.mean_lum = 150;
disp_test_params.persistence = 0.1;
disp_test_params.srx = 1;
disp_test_params.sry = 1;
disp_test_params.mipmap_level = 1;

% Linearly fit and report error
for vv=1:length(velocities)
    vel = velocities(vv);
    disp_test_params.velocity = vel;
    for rr = 1:size(mis,1)
        fps = mis(rr,1,vv);
        disp_test_params.fps = fps;
        camojab_q(rr,2,vv) = camojab(texture, disp_test_params); 
    end
end
camojab_q(:,2,:) = camojab_q(:,2,:) - camojab_q(end,2,end);

X = reshape(camojab_q(:,2,:),[7,4]);
Y = reshape(mis(:,2,:),[7,4]);
plcc = corr(X(:),Y(:),'type','Pearson');
srocc = corr(X(:),Y(:),'type','Spearman');
disp("PLCC = " + plcc + " | SROCC = " + srocc);
[c1,c2] = fit_camojab_to_mackins_data(camojab_q,mis);
camojab_q(:,2,:) = camojab_q(:,2,:) * c1 + c2;
disp("Error (CaMoJAB) = " + sqrt(mean(reshape((camojab_q(:,2,:) - mis(:,2,:)) .^ 2,[1 , numel(mis(:,2,:))]))));

% Densly sample
all_rr = logspace(1.6,3.301);%[50:5:100,120:20:200, 250:50:500,600:100:1000,1200:200:2000];
all_camojab_q = zeros([length(velocities),length(all_rr)]);
for vv=1:length(velocities)
    vel = velocities(vv);
    disp_test_params.velocity = vel;
    for rr = 1:length(all_rr)
        fps = all_rr(rr);
        disp_test_params.fps = fps;
        all_camojab_q(vv,rr) = camojab(texture, disp_test_params); 
    end
end
all_camojab_q = all_camojab_q - all_camojab_q(end,end);
all_camojab_q = all_camojab_q* c1 + c2;

% Plot
clf
lines = ['r','k','b','m'];
lineColors = ["#4DBEEE"	,"#D95319",  "#EDB120", "#7E2F8E"];
for vv = 1:length(velocities)
    semilogx(all_rr,all_camojab_q(vv,:),'-','Color',lineColors(vv),'LineWidth',1);hold on;
    semilogx((mis(:,1,vv)),mis(:,2,vv),'x','Color',lineColors(vv),'LineWidth',1);hold on;
    xticks(mis(:,1,vv))
    axis([0 2000 0 5])
    errorbar((mis(:,1,vv)),mis(:,2,vv),mis(:,2,vv) - mis_err_low(:,1,vv),mis_err_high(:,1,vv) - mis(:,2,vv),'x','Color',lineColors(vv),'HandleVisibility','off','LineWidth',1);hold on;    
end
hold off;
set(gca,'TickLabelInterpreter','latex','FontSize', 16)
xlabel('Frame Rate(Hz)','FontSize', 20,'Interpreter','latex')
ylabel('Mean Impairment Score','FontSize', 20,'Interpreter','latex')
title('Checkerboard $\vert$ 60 ppd $\vert$ 150 $cd/m^2$ $\vert$ Persistence = 0.1','FontSize', 18,'Interpreter','latex')
grid on

legend({'10 deg/s(CaMoJAB)','10 deg/s(Mackin et al. 2016)','30 deg/s(CaMoJAB)','30 deg/s(Mackin et al. 2016)','50 deg/s(CaMoJAB)','50 deg/s(Mackin et al. 2016)','70 deg/s(CaMoJAB)','70 deg/s(Mackin et al. 2016)'},'Location','southeast','FontSize', 13,'Interpreter','latex')

function [c1,c2] = fit_camojab_to_itu_data(cmarrrq,mis)
    fun = @(vars) sqrt(mean(reshape(((cmarrrq.*vars(1) + vars(2)) - mis) .^ 2,[1 , numel(mis)])));
    [vars,rmse] = fminsearch(fun,[-1,0]);
    c1 = vars(1);c2=vars(2);
    disp("rmse = " + rmse);
end

function [c1,c2] = fit_camojab_to_mackins_data(marrr,mis)
%     fun = @(vars) sqrt(mean(reshape(((marrr(:,2,1:3).*vars(1) + vars(2)) - mis(:,2,1:3)) .^ 2,[1 , numel(mis(:,2,1:3))])));
    fun = @(vars) sqrt(mean(reshape(((marrr(:,2,:).*vars(1) + vars(2)) - mis(:,2,:)) .^ 2,[1 , numel(mis(:,2,:))])));
    [vars,rmse] = fminsearch(fun,[-1,0]);
    c1 = vars(1);c2=vars(2);
    disp("rmse = " + rmse);
end