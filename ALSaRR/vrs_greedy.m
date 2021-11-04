function [vrs_img,Q] = vrs_greedy(B, b, q, N, numBlocksX, numBlocksY)
%GREEDYVRS Calculates close to optimal VRS map
%   Calculates the optimal VRS map given the bandwidth budget and all
%   quality values of all possible VRS states
% B     : Maximum allowed bandwidth (scalar)
% b     : bandwidth of each shading rate of each tile (Nx7 array)
% q     : quality of each shading rate of each tile (Nx7 array)
% n     : total number of tiles
% (numBlocksX,numBlocksY)  : resolution of VRS map
% vrs_img : colour coded VRS map 
% Q    : sum of quality of each tile at selected state

% Refer to paper for an explanation of the below code

% Initialize an array of ratios of each possible shading rates of each tile
ratio = zeros([N*(size(b,2)-1),3]);
count = 1;
for j = 1:N
    for k = 2:size(b,2)
        ratio(count,:) = [q(j,k) / b(j,k), j, k]; % 2nd and 3rd enteries to track the tiles
        count = count + 1;
    end
end

% Sort ration array in non-increasing order
[~,idx] = sort(ratio(:,1),'descend'); % sort just the first column
ratio = ratio(idx,:);   % sort the whole matrix using the sort indices

% Assign lowest shading rate (4x4) to each tile and update remaining budget
selection = ones([N,1]);
c = B - N; % remaining total bandwidth
r = 1; % moving pointer
Q = 0; % total quality
for j = 1:N
    Q = Q + q(j,1); % assigned 4x4 to every tile
end

% Iterate through ratio array and keep updating the shading rate of highest
% ratio tile till we run out of remaining bandwidth
while r < size(ratio,1)
    if q(ratio(r,2),ratio(r,3)) > q(ratio(r,2),selection(ratio(r,2)))
        c = c - b(ratio(r,2),ratio(r,3)) + b(ratio(r,2),selection(ratio(r,2)));
        if c < 0
            break;
        end
        Q = Q + q(ratio(r,2),ratio(r,3)) - q(ratio(r,2),selection(ratio(r,2)));
        selection(ratio(r,2)) = ratio(r,3);        
    end
    r = r + 1;
end

% Generate colour coded VRS map
vrs_img = zeros([numBlocksY,numBlocksX,3]);
k = 1;
for j = 1:numBlocksX
    for i = 1:numBlocksY
        if selection(k) == 1   %4x4 
            vrs_img(i,j,:) = [0.9294117647058824, 0.1098039215686275, 0.1411764705882353].^ 2.2;
        elseif selection(k) == 2 % 2x4
            vrs_img(i,j,:) = [0.9333333333333333, 0.5294117647058824, 0.1333333333333333].^2.2;
        elseif selection(k) == 3 % 4x2
            vrs_img(i,j,:) = [0.9843137254901961, 0.6901960784313725, 0.2509803921568627].^2.2;
        elseif selection(k) == 4 % 2x2
            vrs_img(i,j,:) = [0.5529411764705882, 0.7764705882352941, 0.2470588235294118].^2.2;
        elseif selection(k) == 5 % 1x2
            vrs_img(i,j,:) = [0.1686274509803922, 0.7137254901960784, 0.4509803921568627].^2.2;
        elseif selection(k) == 6 % 2x1
            vrs_img(i,j,:) = [0, 0.6549019607843137, 0.615686274509803].^2.2;
        elseif selection(k) == 7 % 1x1
            vrs_img(i,j,:) = [0.1098039215686275, 0.4588235294117647, 0.7372549019607843].^2.2;
        else                % Invalid vrs state
            vrs_img(i,j,:) = [0,0,0];
        end
k = k + 1;
    end
end

end
