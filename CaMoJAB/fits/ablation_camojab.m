ablation_masks = [[1 1 1]; [1 1 0]; [0 0 1]; [0 1 0]; [1 0 0]; [1 0 1]; [0 1 1]]; %;

for ii = 1:length(ablation_masks)
    ablation_mask = ablation_masks(ii,:);
    [ws,wt, beta] = train_camojab(ablation_mask);
    fprintf("[%d %d %d] ws = %f | wt = %f | beta = %f\n",ablation_mask(1),ablation_mask(2),ablation_mask(3),ws,wt,beta);
end