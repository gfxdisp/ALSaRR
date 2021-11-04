% This script scales the pairwise comparision data to JND units
% Clone and add this repository to your matlab path:
% https://github.com/mantiuk/pwcmp

% Experiment details
observers = ["aj577", "rkm38", "am2442", "fz261", "pmh64","ma912","fg405","aoy20", "am2806", "frs29", "R02","R01","R07"];%, "aj577_2", "aj577_3", "aj577_4", "aj577_5"];
configs = ["BlockPC", "BlockMobile", "BlockVR"];
resolutions = [1,1/2,1/4;1,1/2,1/4;1,1/2,1/4];
texture_types = 0:3;
velocities = [0,10,30;3,10,20;10,45,75];
M = zeros([9,9,size(velocities,2),length(configs)]);
M_observers = {{},{},{};{},{},{};{},{},{}};

% Load results into table
exp_results = readtable('vrs_exp_data.csv');

% Build comparision matrix per observer from table
for cc = 1:length(configs)
    for oo = 1:length(observers)     
        tab = exp_results(strcmp(t.Config,configs(cc))~=0 & strcmp(t.Observer,observers(oo))~=0,:);
        if height(tab) == 0
            continue;
        end
        % initialize comparision matrix
        for vv=1:size(velocities,2)    
            if isempty(M_observers{cc,vv})
                M_observers{cc,vv} = zeros([9,9]);
            end
            M_observers{cc,vv}(:,:,oo) = zeros([9,9]);
        end
        for ii = 1:height(tab)
            c1 = stim2c(tab(ii,:).Stimulus1_type, tab(ii,:).Stimulus1_resolution);
            c2 = stim2c(tab(ii,:).Stimulus2_type, tab(ii,:).Stimulus2_resolution);
            assert(tab(ii,:).Stimulus1_xvelolicty == tab(ii,:).Stimulus2_xvelolicty);
            vv = find(velocities(cc,:) == tab(ii,:).Stimulus1_xvelolicty,1);
            if(tab(ii,:).User_Selection == 1)
                M(c1,c2,vv,cc) = M(c1,c2,vv,cc) + 1;
                M_observers{cc,vv}(c1,c2,oo) = M_observers{cc,vv}(c1,c2,oo)+1;
            elseif(tab(ii,:).User_Selection == 2)
                M(c2,c1,vv,cc) = M(c2,c1,vv,cc) + 1;
                M_observers{cc,vv}(c2,c1,oo) = M_observers{cc,vv}(c2,c1,oo)+1;
            end
        end
    end
end

% Initialize JND array
jnds = zeros([9 length(velocities) length(configs)]);
GT = zeros([4,3,length(velocities) length(configs)]);
ErrorLow = GT;
ErrorHigh = GT;

% Scale PWC matrix to JND scale. Generate confidence intervals by
% bootstrapping
for cc=1:length(configs)
    for vv=1:size(velocities,2)
        [jnds(:,vv,cc),~] = pw_scale(M(:,:,vv,cc));
        jnds(:,vv,cc) = jnds(:,vv,cc) - jnds(1,vv,cc);
        j = jnds(:,vv,cc);
        GT(:,:,vv,cc) = [[j(1),j(2),j(3)];[j(1),j(4),j(5)];[j(1),j(6),j(7)];[j(1),j(8),j(9)]];        
        M_boot = reshape(M_observers{cc,vv},[1 81 size(M_observers{cc,vv},3)]); M_boot = squeeze(M_boot);M_boot = M_boot';
        [jnds_boot,stats] = pw_scale_bootstrp(M_boot,100,{'alpha',0.05,'prior','bounded'});
        stats.jod_high = stats.jod_high' - jnds_boot(1);
        stats.jod_low = stats.jod_low' - jnds_boot(1);
        jnds_boot = jnds_boot' - jnds_boot(1);
        j = jnds_boot;
        GT(:,:,vv,cc) = [[j(1),j(2),j(3)];[j(1),j(4),j(5)];[j(1),j(6),j(7)];[j(1),j(8),j(9)]];
        el = abs(stats.jod_low - jnds_boot);
        ErrorLow(:,:,vv,cc) = [[el(1),el(2),el(3)];[el(1),el(4),el(5)];[el(1),el(6),el(7)];[el(1),el(8),el(9)]];
        eh = abs(stats.jod_high - jnds_boot);
        ErrorHigh(:,:,vv,cc) = [[eh(1),eh(2),eh(3)];[eh(1),eh(4),eh(5)];[eh(1),eh(6),eh(7)];[eh(1),eh(8),eh(9)]];
    end
end

% Plot results
pltNdx = 1;
for cc=1:length(configs)
    for vv=1:size(velocities,2)
        subplot(3,3,pltNdx);
        for tt = 1:4
            errorbar(resolutions(cc,:),GT(tt,:,vv,cc),ErrorLow(tt,:,vv,cc),ErrorHigh(tt,:,vv,cc));hold on; 
        end
        hold off;
        xlabel("Shading Resolution")
        ylabel("JND")
        titleStr = sprintf("Config%d | %d deg/s",cc,velocities(cc,vv));
        title(titleStr)
        axis([0.1 1 -12 2])
        legend({'Checkerboard', 'Gradient', 'Grass', 'Noise'},'Location','southeast')
    pltNdx = pltNdx + 1;
    end
end

% save vrs_exp_data_jnd.mat GT ErrorLow ErrorHigh;


function cn = res2c(res)
    if res == 1
        cn = 0;
    elseif res == .5
        cn = 1;
    elseif res == .25
        cn = 2;
    end
end

function texid = tex2c(tex)
    if strcmp(tex,'Checkerboard')
        texid = 0;
    elseif strcmp(tex,'Gradient')
        texid = 1;
    elseif strcmp(tex,'Grass')
        texid = 2;
    elseif strcmp(tex,'Noise')
        texid = 3;
    end
end

%Maps the stimulus to 9x9 PWC matrix
function cn = stim2c(tex_type, res)
    if res == 1
        cn = 1;
    else
        cn = 1 + 2 * tex2c(tex_type) + res2c(res);
    end
end