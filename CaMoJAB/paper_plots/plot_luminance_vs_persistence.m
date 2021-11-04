if ~exist('model_Dst','file')
    addpath('../')
end

mean_lum = logspace(-2,3,20); % cd/m^2
persistence = 0.05:0.05:1;
fps = 60;
ppd = 60;
vel = 10;
ref_img = load('../../Experiment_Data/vrs_exp_stimulus.mat');ref_img = ref_img.grass;

qt = zeros([length(mean_lum) length(persistence)]);
qs = zeros([length(mean_lum) length(persistence)]);
load('../camojab_fits.mat')

for ll = 1:length(mean_lum)
    for pp = 1:length(persistence)
        [qs(ll,pp),qt(ll,pp),~] = model_Dst(ref_img, mean_lum(ll), fps, vel, persistence(pp), ppd, 1, 1, 2, beta);
    end
end

camojab_preds =  qt * wt + qs * ws;
camojab_preds = camojab_preds - camojab_preds(end,end);

figure(1)
[X,Y] = meshgrid(persistence,mean_lum);
surf(X,Y,camojab_preds);
set(gca,'YScale','log')
xlabel('Persistence')
ylabel('Luminance (cd/m^2)')
zlabel('\DeltaQ (JND)')

figure(2)
pers = [0.05, 0.1,0.25,0.5,0.75,1];
for pp = 1:length(pers)
   pNdx = find(persistence==pers(pp));
   plot(mean_lum,camojab_preds(:,pNdx),'DisplayName',"$p^{[d]} = $"+pers(pp),'LineWidth',1); hold on; 
end
set(gca,'XScale','log')
hold off;
set(gca,'TickLabelInterpreter','latex','FontSize', 18)
xlabel('Luminance ($cd/m^2$)','FontSize', 20,'Interpreter','latex')
ylabel('$\Delta Q$ (JND)','FontSize', 20,'Interpreter','latex')
titleStr = sprintf('Grass $$\\vert$$ %d deg/s $$\\vert$$ %d Hz $$\\vert$$ %d ppd', vel,fps,ppd);
title(titleStr,'FontSize', 20,'Interpreter','latex')
legend('Location','best','FontSize', 18,'Interpreter','latex')
grid on
xlim([0.01,1000])
% exportgraphics(gcf,"figs/lum_vs_per.pdf",'ContentType','vector')

figure(3)
for ll = [1:4:length(mean_lum),length(mean_lum)]
   lum = mean_lum(ll);
   if(mod(lum,1) == 0)
       lumStr = num2str(lum);
   else
       lumStr = sprintf('%.2f',lum);
   end
   plot(persistence,camojab_preds(ll,:),'DisplayName',lumStr+" $cd/m^2$",'LineWidth',1); hold on; 
end
hold off;
set(gca,'TickLabelInterpreter','latex','FontSize', 18)
xlabel('Persistence','FontSize', 20,'Interpreter','latex')
ylabel('$\Delta Q$ (JND)','FontSize', 20,'Interpreter','latex')
titleStr = sprintf('Grass $$\\vert$$ %d deg/s $$\\vert$$ %d Hz $$\\vert$$ %d ppd', vel,fps,ppd);
title(titleStr,'FontSize', 20,'Interpreter','latex')
legend('Location','best','FontSize', 18,'Interpreter','latex')
grid on
xlim([0.05 1])
% exportgraphics(gcf,"figs/per_vs_lum.pdf",'ContentType','vector')