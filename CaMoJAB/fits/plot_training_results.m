function plot_training_results(camojab_vrs, camojab_marrr, ws,wt, beta)
%plot_training_results Plot the model fits to training data
%   Detailed explanation goes here

%--------- Plot CaMoJAB predictions for VRS Experiment Results-------------
figure(1)
GT = load('../../Experiment_Data/vrs_exp_data_jnd.mat');ErrorLow = GT.ErrorLow;ErrorHigh=GT.ErrorHigh; GT = GT.GT;
resolutions = [1,1/2,1/4;1,1/2,1/4;1,1/2,1/4];
velocities = [0,10,30;3,10,20;10,45,75];
pltNdx = 1;
lines = ['r','k','b','m'];
configs = ["Block-PC";"Block-Mobile";"Block-VR"];
lineColors = ["#4DBEEE"	,"#D95319",  "#EDB120", "#7E2F8E"];
for cc=1:size(camojab_vrs,4)
    for vv=1:size(camojab_vrs,3)
        subplot(3,3,pltNdx);
        for tt = 1:4
            plot([1,1/2,1/4],camojab_vrs(tt,:,vv,cc),'--','Color',lineColors(tt));hold on;
            e = errorbar(resolutions(cc,:),GT(tt,:,vv,cc),ErrorLow(tt,:,vv,cc),ErrorHigh(tt,:,vv,cc),'-x','Color',lineColors(tt));hold on; 
        end
        hold off;
        set(gca,'TickLabelInterpreter','latex','FontSize',12)
        xlabel("Shading Resolution",'FontSize', 15,'Interpreter','latex')
        ylabel("$\Delta Q$ (JND)",'FontSize', 15,'Interpreter','latex')
        titleStr = sprintf("%s $\\vert$ %d deg/s",configs(cc),velocities(cc,vv));
        title(titleStr,'FontSize', 15,'Interpreter','latex')        
        axis([0.1 1 -9 2])        
        grid on
    pltNdx = pltNdx + 1;
    end
end
legend({'Checkerboard (CaMoJAB)','Checkerboard (Exp.)', 'Gradient (CaMoJAB)','Gradient (Exp.)', 'Grass (CaMoJAB)','Grass (Exp.)', 'Noise (CaMoJAB)', 'Noise (Exp.)'},'Location','southeast','FontSize', 13,'Interpreter','latex')

%--------- Plot CaMoJAB predictions for MARRR Experiment 1 Results-----------
figure(2)
marrr_data = load('../../Experiment_Data/external_data/marrr_exp1_jnds.mat');
marrr_refresh_rates = 50:5:165;
marrr_velocities = [15,30,45];
jnds = [marrr_data.res.JODs(1,:); marrr_data.res.JODs(2,:); marrr_data.res.JODs(3,:)];
jnds = jnds - jnds(:,1);
errors = load('../../Experiment_Data/external_data/marrr_exp1_confidence75.mat');
error_low(:,:) = errors.errors(1,1:3,:);
error_high(:,:) = errors.errors(2,1:3,:);
pltNdx =1;
for vv=1:length(marrr_velocities)
    subplot(1,3,pltNdx);
    errorbar(marrr_refresh_rates,jnds(vv,:),-error_low(vv,:),error_high(vv,:),'-r','DisplayName','Denes et al. 2020 (Exp.)'); hold on;
    plot(marrr_refresh_rates, camojab_marrr(vv,:),'--r');hold off;
    if vv==1
        xlabel('Refresh Rate (Hz)','FontSize', 15,'Interpreter','latex')
        ylabel('$\Delta Q$ (JND)','FontSize', 15,'Interpreter','latex')
    else
        xlabel('Refresh Rate (Hz)','FontSize', 15,'Interpreter','latex','Color',[1 1 1])
        ylabel('$\Delta Q$ (JND)','FontSize', 15,'Interpreter','latex','Color',[1 1 1])
    end
    title(marrr_velocities(vv)+" deg/s",'FontSize', 15,'Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex','FontSize',15)
    grid on
    if vv == 3
        legend({'MARRR (Exp. 1)', 'CaMoJAB'},'Location', 'southeast','FontSize', 15,'Interpreter','latex');
    end
    pltNdx = pltNdx + 1;
end

%------------ Plot CaMoJAB predictions for VRS Experiment Results (Continuous plots)-----------
figure(3)
plot_vrs_continous(ws,wt, beta)

end

function plot_vrs_continous(ws,wt, beta)
%PLOT_VRS_CONTINOUS Summary of this function goes here
imgs = load('../../Experiment_Data/vrs_exp_stimulus.mat');

% Model predictions for VRS Experiment
configs = ["Config1","Config2","Config3"];
velocities = [0,10,30;3,10,20;10,45,75];
resolutions = [1:-0.05:0.2;1:-0.05:0.2;1:-0.05:0.2];
stimuli = ["checkerbox","gradient","grass","noise"];

camojab_es_vrs = zeros([length(stimuli) size(resolutions,2) size(velocities,2) length(configs)]);
camojab_et_vrs = zeros([length(stimuli) size(resolutions,2) size(velocities,2) length(configs)]);

for cc=1:length(configs)    
    if strcmp(configs(cc), 'Config1')==1    
        ppd = 50;
        fps = 144;
        persistence = 1;
        mean_lum = 150;
        mipmap_level_ref = 3;
    elseif strcmp(configs(cc), 'Config2')==1   
        ppd = 90;
        fps = 60;
        persistence = 1;
        mean_lum = 75;
        mipmap_level_ref = 1;
    elseif strcmp(configs(cc), 'Config3')==1    
        ppd = 26; % 2x SSAA
        fps = 144;
        persistence = 0.05;
        mean_lum = 2.5;
        mipmap_level_ref = 5;
    end
    srx = resolutions(cc,:);
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
            [camojab_es_vrs(ii,:,vv,cc),camojab_et_vrs(ii,:,vv,cc),~] = model_Dst(ref_img, mean_lum, fps, vel, persistence, ppd, srx, sry,mipmap_level_ref,beta,[]);
        end
    end
end
camojab_vrs = ws*camojab_es_vrs + wt * camojab_et_vrs;
camojab_vrs = camojab_vrs - camojab_vrs(:,1,:,:);

GT = load('../../Experiment_Data/vrs_exp_data_jnd.mat');ErrorLow = GT.ErrorLow;ErrorHigh=GT.ErrorHigh; GT = GT.GT;
pltNdx = 1;
lines = ['r','k','b','m'];
lineColors = ["#4DBEEE"	,"#D95319",  "#EDB120", "#7E2F8E"];
configNames = ["PC","Mobile","VR"];
configRR = [144,60,144];
configPPD = [50,90,15];
configLum = ["150","75","2.5"];
for cc=1:size(camojab_vrs,4)
    for vv=1:size(camojab_vrs,3)
        subplot(3,3,pltNdx);
        for tt = 1:4
            plot(resolutions(cc,:),camojab_vrs(tt,:,vv,cc),'--','Color',lineColors(tt));hold on;
            e = errorbar([1 1/2 1/4],GT(tt,:,vv,cc),ErrorLow(tt,:,vv,cc),ErrorHigh(tt,:,vv,cc),'-x','Color',lineColors(tt));hold on; 
        end
        hold off;
        set(gca,'TickLabelInterpreter','latex','FontSize',12)
        if cc==3
            xlabel("Shading Resolution",'FontSize', 15,'Interpreter','latex')
        else
            xlabel("Shading Resolution",'FontSize', 15,'Interpreter','latex','Color',[1 1 1])
            set(gca, 'XTickLabel', [])
        end
        if vv==1
            ylabel("$\Delta Q$ (JND)",'FontSize', 15,'Interpreter','latex')
        else
            ylabel("$\Delta Q$ (JND)",'FontSize', 15,'Interpreter','latex','Color',[1 1 1])
            set(gca, 'YTickLabel', [])
        end
        titleStr = sprintf("Block-%s $\\vert$ %d Hz$\\vert$ %d PPD $\\vert$ %s cd/m$^2$ $\\vert$ %d deg/s",configNames(cc),configRR(cc),configPPD(cc),configLum(cc),velocities(cc,vv));%Config2 (Mobile) | 60Hz | 90 PPD | 75 cd/m2 | 3 deg/s
        title(titleStr,'FontSize', 15,'Interpreter','latex')
        axis([0.1 1 -10 2])        
        grid on
    pltNdx = pltNdx + 1;
    end
end
legend({'Checkerboard (CaMoJAB)','Checkerboard (Exp.)', 'Gradient (CaMoJAB)','Gradient (Exp.)', 'Grass (CaMoJAB)','Grass (Exp.)', 'Noise (CaMoJAB)', 'Noise (Exp.)'},'Location','southeast')

end

