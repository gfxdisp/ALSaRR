mean_lum = 150; % cd/m^2
persistence = 1;
fps = 30:3:120;
ppds = 10:5:60;
velocities = [1 15 45];
ref_img = load('../../Experiment_Data/vrs_exp_stimulus.mat');ref_img = ref_img.grass;

qt = zeros([length(velocities) length(ppds) length(fps)]);
qs = zeros([length(velocities) length(ppds) length(fps)]);
load('../camojab_fits.mat')

for ff = 1:length(fps)
    for vv = 1:length(velocities)
        [qs(vv,:,ff),qt(vv,:,ff),~] = model_Dst(ref_img, mean_lum, fps(ff), velocities(vv), persistence, ppds(end), ppds./ppds(end), ppds./ppds(end), 3, beta);
    end
end

camojab_preds = qt * wt + qs * ws;
camojab_preds = camojab_preds - camojab_preds(end,1,1);

figure(1)

[X,Y] = meshgrid(fps,ppds);
col = ['#0072BD';'#EDB120';'#77AC30'];
for vv = 1:length(velocities)
    surf(X,Y./ppds(end),squeeze(camojab_preds(vv,:,:)),'DisplayName',velocities(vv)+" deg/s",'EdgeColor','none','FaceColor',col(vv,:),'FaceAlpha',0.5); hold on;    
end
hold off;
xlabel('Refresh Rate (Hz)','Interpreter','latex','FontSize',12,'Rotation',15)
ylabel('Shading Resolution','Interpreter','latex','FontSize',12,'Rotation',-25)
zlabel('$\Delta Q$ (JND)','Interpreter','latex','FontSize',12)
legend('Location','best','Interpreter','latex','FontSize',11)
set(gca,'TickLabelInterpreter','latex','FontSize',12)
ylim([ppds(1)/ppds(end) Inf])
zlim([-1 Inf])
% exportgraphics(gcf,"figs/res_vs_vel_vs_rr_3d_2.pdf",'ContentType','vector')

figure(2)
res = [15 20 30 45 60];%;10:10:60;
for rr = 1:length(res)
    ndx = find(ppds==res(rr));
    plot(fps,squeeze(camojab_preds(2,ndx,:)),'DisplayName',"$r_x,r_y = $ " + round(res(rr)./ppds(end),2)); hold on;
end
hold off
xlabel('Refresh Rate (Hz)','FontSize', 15,'Interpreter','latex')
ylabel('$\Delta Q$ (JND)','FontSize', 15,'Interpreter','latex')
title(velocities(2)+" deg/s",'FontSize', 15,'Interpreter','latex')
legend('Location','best','FontSize', 12,'Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','FontSize',15)
xlim([30,120])
grid on
% exportgraphics(gcf,"figs/res_vs_vel_vs_rr_2d_2.pdf",'ContentType','vector')

figure(3)
refrate = [30,45,60,75,90,120];
for rr = 1:length(refrate)
    ndx = find(fps==refrate(rr));
    plot(ppds./ppds(end),squeeze(camojab_preds(2,:,ndx)),'DisplayName',refrate(rr)+" Hz"); hold on;
end
hold off
xlabel('Shading Resolution','FontSize', 15,'Interpreter','latex')
ylabel('$\Delta Q$ (JND)','FontSize', 15,'Interpreter','latex')
title(velocities(2)+" deg/s",'FontSize', 15,'Interpreter','latex')
legend('Location','best','FontSize', 14,'Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','FontSize',15)
grid on
xlim([ppds(1)/ppds(end) Inf])
% exportgraphics(gcf,"figs/res_vs_vel_vs_rr_2d2_2.pdf",'ContentType','vector')