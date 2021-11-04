% This script visualizes validation experiment results as a heatmap

if ~exist('myBinomTest','file')
    addpath('../utils/')
end

observer = 'all';
t = readtable('alsarr_validation_experiment.csv');
if strcmp(observer,'all') ~= 1
    t = t(string(t.id)==observer ,:);
end


disp("% Wins = " + sum(t.answer)/height(t));

% Separate Motion Paths
Conditions = {"6.25% | MARRR" , "6.25% | 4x4 120FPS", "6.25% | 2x2 30FPS", "12.5% | MARRR" , "12.5% | 4x2 120FPS", "12.5% | 2x2 60FPS"};
Scenes = {"Scifi-P","Scifi-R","Scifi-F","Scifi-U", "Sally-P", "Sally-R", "Sally-F", "Sally-U", "JoLem-P", "JoLem-R", "JoLem-F", "JoLem-U"};

total_trials = zeros(length(Scenes),length(Conditions));
passed_trials = zeros(length(Scenes),length(Conditions));

for ii = 1:height(t)
    sceneID = sceneName2Id(t(ii,:).scene{1});
    conditionID = condition2Id(t(ii,:).budget, t(ii,:).against, t(ii,:).uniformSetting);
    total_trials(sceneID,conditionID) = total_trials(sceneID,conditionID) + 1;
    if t(ii,:).answer == 1
        passed_trials(sceneID,conditionID) = passed_trials(sceneID,conditionID) + 1;
    end    
end

prob = passed_trials ./ total_trials;
figure(1)
heatmap(Conditions,Scenes,prob)
title(observer)

% Aggregated Motion Paths
Scenes = {{'Scifi'}, {'Country'}, {'Haunted'}};

total_trials = zeros(length(Scenes),length(Conditions));
passed_trials = zeros(length(Scenes),length(Conditions));

for ii = 1:height(t)
    [~,sceneID] = sceneName2Id(t(ii,:).scene{1});
    conditionID = condition2Id(t(ii,:).budget, t(ii,:).against, t(ii,:).uniformSetting);
    total_trials(sceneID,conditionID) = total_trials(sceneID,conditionID) + 1;
    if t(ii,:).answer == 1
        passed_trials(sceneID,conditionID) = passed_trials(sceneID,conditionID) + 1;
    end    
end

prob = passed_trials ./ total_trials;
figure(2)
heatmap(Conditions,Scenes,prob)
title(observer)
title('Probability of picking our method')
xlabel('Budget | Rendering method')
ylabel('Scene')



function [id,sid] = sceneName2Id(sceneName)
    switch sceneName
        case 'SciFi-FlyThroughParallax'
            id = 1;
            sid = 1;
        case 'SciFi-FlyThroughRandom'
            id = 2;
            sid = 1;
        case 'SciFi-FlyThroughForward'
            id = 3;
            sid = 1;
        case 'SciFi-FirstPerson'
            id = 4;
            sid = 1;
        case 'Sally-FlyThroughParallax'
            id = 5;
            sid = 2;
        case 'Sally-FlyThroughRandom'
            id = 6;
            sid = 2;
        case 'Sally-FlyThroughForward'
            id = 7;
            sid = 2;
        case 'Sally-FirstPerson'
            id = 8;
            sid = 2;
        case 'LemonChanged-FlyThroughParallax'
            id = 9;
            sid = 3;
        case 'LemonChanged-FlyThroughRandom'
            id = 10;
            sid = 3;
        case 'LemonChanged-FlyThroughForward'
            id = 11;
            sid = 3;
        case 'LemonChanged_ThirdPerson'
            id = 12;           
            sid = 3;
    end            
end

function  id = condition2Id(budget, against, vrs_setting)
    if budget == 0.0625
        if strcmp(against,'MARRR')
            id = 1;
        else
            if strcmp(vrs_setting,'X1_PER_4X4_PIXELS')
                id = 2;
            else
                id = 3;
            end
        end
    elseif budget == .125
        if strcmp(against,'MARRR')
            id = 4;
        else
            if strcmp(vrs_setting,'X1_PER_4X2_PIXELS')
                id = 5;
            else
                id = 6;
            end
        end
    end
end