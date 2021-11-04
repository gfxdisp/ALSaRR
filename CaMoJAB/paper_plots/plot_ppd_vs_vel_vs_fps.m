if ~exist('model_Dst','file')
    addpath('../')
end

mean_lum = 150; % cd/m^2
persistence = 1;
fps = [30 60 120];
ppds = 10:2:60;
velocities = [0:1:9 10:2:60];
ref_img = load('../../Experiment_Data/vrs_exp_stimulus.mat');ref_img = ref_img.grass;

qt = zeros([length(velocities) length(ppds) length(fps)]);
qs = zeros([length(velocities) length(ppds) length(fps)]);
load('../camojab_fits.mat')

for ff = 1:length(fps)
    for vv = 1:length(velocities)
        [qs(vv,:,ff),qt(vv,:,ff),~] = model_Dst(ref_img, mean_lum, fps(ff), velocities(vv), persistence, ppds(end), ppds./ppds(end), ppds./ppds(end), 1, beta,[]);
    end
end

camojab_preds = qt * wt + qs * ws;
camojab_preds = camojab_preds - camojab_preds(end,1,1);

figure(1)
[X,Y] = meshgrid(ppds,velocities);
col = ['#0072BD';'#EDB120';'#77AC30'];
for ff = 1:length(fps)
    surf(X./ppds(end),Y,camojab_preds(:,:,ff),'DisplayName',fps(ff)+" Hz",'EdgeColor','none','FaceColor',col(ff,:),'FaceAlpha',0.5); hold on;    
end
hold off;
set(gca,'TickLabelInterpreter','latex','FontSize',12)
xlabel('Shading Resolution','Interpreter','latex','FontSize',14,'Rotation',15)
ylabel('Velocity (deg/s)','Interpreter','latex','FontSize',14,'Rotation',-25)
zlabel('$\Delta Q$ (JND)','Interpreter','latex','FontSize',14)
legend('Location','best','Interpreter','latex','FontSize',11)
zlim([-2 Inf])
xlim([ppds(1)/ppds(end) Inf])
% exportgraphics(gcf,"figs/res_vs_vel_vs_rr_3d.pdf",'ContentType','vector')

figure(2)
res = [10:10:60];
for rr = 1:length(res)
    ndx = find(ppds==res(rr));
    plot(velocities,camojab_preds(:,ndx,2),'DisplayName',res(rr)+" PPD"); hold on;
end
hold off
xlabel('Velocity (deg/s)','FontSize', 15,'Interpreter','latex')
ylabel('$\Delta Q$ (JND)','FontSize', 15,'Interpreter','latex')
title(fps(2)+" Hz",'FontSize', 15,'Interpreter','latex')
legend('Location','best','FontSize', 15,'Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','FontSize',15)
grid on
% exportgraphics(gcf,"figs/res_vs_vel_vs_rr_2d.pdf",'ContentType','vector')

figure(3)
vel = [0:10:50];
for vv = 1:length(vel)
    ndx = find(velocities==vel(vv));
    plot(ppds./ppds(end),camojab_preds(ndx,:,2),'DisplayName',vel(vv)+" deg/s",'LineWidth',1); hold on;
end
hold off
set(gca,'TickLabelInterpreter','latex','FontSize',15)
xlabel('Shading Resolution','FontSize', 20,'Interpreter','latex')
ylabel('$\Delta Q$ (JND)','FontSize', 20,'Interpreter','latex')
title(fps(2)+" Hz",'FontSize', 20,'Interpreter','latex')
legend('Location','best','FontSize', 18,'Interpreter','latex')
grid on
xlim([ppds(1)/ppds(end) Inf])
% exportgraphics(gcf,"figs/res_vs_vel_vs_rr_2d.pdf",'ContentType','vector')