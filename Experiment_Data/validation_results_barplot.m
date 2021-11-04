% This script visualizes validation experiment results as a bar plot

if ~exist('myBinomTest','file')
    addpath('../utils/')
end

observer = 'all';
t = readtable('alsarr_validation_experiment.csv');
if strcmp(observer,'all') ~= 1
    t = t(string(t.id)==observer ,:);
end

% aggregate data
num_trials = zeros([2 3]);
num_success = zeros([2 3]);

Conditions = {{"MARRR" , "4x4 120FPS", "2x2 30FPS"}; {"MARRR" , "4x2 120FPS", "2x2 60FPS"}};

for ii = 1:height(t)
    conditionID = condition2Id(t(ii,:).budget, t(ii,:).against, t(ii,:).uniformSetting);
    num_trials(conditionID(1),conditionID(2)) = num_trials(conditionID(1),conditionID(2)) + 1;
    if t(ii,:).answer == 1
        num_success(conditionID(1),conditionID(2)) = num_success(conditionID(1),conditionID(2)) + 1;
    end
end

figure(1)

[expected_prob, conf_ints] = binofit(num_success(:),num_trials(:));
expected_prob = expected_prob * 100;
conf_ints = conf_ints * 100;
expected_prob = reshape(expected_prob,[2,3]);
conf_low = expected_prob - reshape(conf_ints(:,1),[2,3]);
conf_high = reshape(conf_ints(:,2),[2,3]) - expected_prob;

b = bar(expected_prob, 'grouped', 'FaceColor','none','EdgeColor',[1 0 0],'LineWidth',1.5);
set(gca,'xticklabel',{'6.25\% (15 MP/s)', '12.5\% (30 MP/s)'},'TickLabelInterpreter','latex','FontSize',12)
yticks(0:25:100)
ylabel('\% Choosing our method','Interpreter','latex','FontSize',15);
xlabel('Bandwidth and reference methods','Interpreter','latex','FontSize',15)
% set(gca,'TickLabelInterpreter','latex','FontSize',12)
grid on
ylim([0,105])

hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(expected_prob);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',expected_prob,conf_low,conf_high,'k','linestyle','none');

for i = 1:nbars
%     myBinomTest(num_success(1,i),num_trials(1,i),0.5,'two')
%     myBinomTest(num_success(2,i),num_trials(2,i),0.5,'two')
    if myBinomTest(num_success(1,i),num_trials(1,i),0.5,'two') <= 0.05
        text(x(i,1)-.015, expected_prob(1,i)+conf_high(1,i) + 3, '*');
    end
    if myBinomTest(num_success(2,i),num_trials(2,i),0.5,'two') <= 0.05
        text(x(i,2)-.015, expected_prob(2,i)+conf_high(2,i) + 3, '*');
    end
end

yline(50,'--k')

for i = 1:nbars    
    text(x(i,1)-.015, 5, Conditions{1}{i}, 'Rotation', 90,'Interpreter','latex','FontSize',12);
    text(x(i,2)-.015, 5, Conditions{2}{i}, 'Rotation', 90,'Interpreter','latex','FontSize',12);
end

hold off

function  id = condition2Id(budget, against, vrs_setting)
    if budget == 0.0625
        if strcmp(against,'MARRR')
            id = [1,1];
        else
            if strcmp(vrs_setting,'X1_PER_4X4_PIXELS')
                id = [1,2];
            else
                id = [1,3];
            end
        end
    elseif budget == .125
        if strcmp(against,'MARRR')
            id = [2,1];
        else
            if strcmp(vrs_setting,'X1_PER_4X2_PIXELS')
                id = [2,2];
            else
                id = [2,3];
            end
        end
    end
end