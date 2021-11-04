function [vrs_img,TQ] = vrs_dp(B, b, q, N, numBlocksX, numBlocksY)
%KNAPSACKVRS Calculates the optimal VRS map
%   Calculates the optimal VRS map given the bandwidth budget and all
%   quality values of all possible VRS states
% B     : Maximum allowed bandwidth (scalar)
% b     : bandwidth of each shading rate of each tile (Nx7 array)
% q     : quality of each shading rate of each tile (Nx7 array)
% n     : total number of tiles
% (numBlocksX,numBlocksY)  : resolution of VRS map
% vrs_img : colour coded VRS map 
% TQ    : sum of quality of each tile at selected state

% See supplementary for explanantion of below code

% add an extra shading rate and an extra tile with 0 bandwidth for ease of
% implementation of memoization table
b = [-1*ones([size(b,1),1]) , b];
q = [-1*ones([size(q,1),1]) , q];
b = [-1*ones([1,size(b,2)]); b];
q = [-1*ones([1,size(q,2)]); q];
% Q(j,c)= Maximum quality possible when only tiles 1 to j are considered 
% and the maximum bandwidth is at most c
Q = zeros([N+1,B+1]); 
% Q_ndx(j,c) = shading rate that gives maximum quality for tile j and
% bandwidth c
Q_ndx = ones([N+1,B+1]);

% Build table Q[][] in bottom up manner
for j = 1:N+1
    for c = 1:B+1
        if j==1
            Q(j,c) = 0; % base case
        elseif c==1
            Q(j,c) = -inf;  % base case
        else
            maxVal = -inf;
            maxNdx = 1;
            for k = 2:size(b,2) % For all shading rates
                if (c-b(j,k) >= 1) && (maxVal < (q(j,k) + Q(j-1,c-b(j,k))))
                    maxVal = q(j,k) + Q(j-1,c-b(j,k));
                    maxNdx = k;
                    Q(j,c) = maxVal;
                    Q_ndx(j,c) = maxNdx;
                end
            end
        end
    end
end

TQ = Q(N+1,B+1);

% Get the assigned shading rates
c = B;
result = zeros([N+1,1]);
for j = N+1:-1:2
    k = Q_ndx(j,c);
    if k<=1
        disp("No Solution")
        break
    end
    result(j) = k;
    c = c - b(j,k);
end

% Generate colour coded VRS map
vrs_img = zeros([numBlocksY,numBlocksX,3]);
k = 2;
for j = 1:numBlocksX
    for i = 1:numBlocksY
        if result(k) == 2   %4x4 
            vrs_img(i,j,:) = [0.9294117647058824, 0.1098039215686275, 0.1411764705882353].^ 2.2;
        elseif result(k) == 3 % 2x4
            vrs_img(i,j,:) = [0.9333333333333333, 0.5294117647058824, 0.1333333333333333].^2.2;
        elseif result(k) == 4 % 4x2
            vrs_img(i,j,:) = [0.9843137254901961, 0.6901960784313725, 0.2509803921568627].^2.2;
        elseif result(k) == 5 % 2x2
            vrs_img(i,j,:) = [0.5529411764705882, 0.7764705882352941, 0.2470588235294118].^2.2;
        elseif result(k) == 6 % 1x2
            vrs_img(i,j,:) = [0.1686274509803922, 0.7137254901960784, 0.4509803921568627].^2.2;
        elseif result(k) == 7 % 2x1
            vrs_img(i,j,:) = [0, 0.6549019607843137, 0.615686274509803].^2.2;
        elseif result(k) == 8 % 1x1
            vrs_img(i,j,:) = [0.1098039215686275, 0.4588235294117647, 0.7372549019607843].^2.2;
        else                % Invalid vrs state
            vrs_img(i,j,:) = [0,0,0];
        end
        k = k + 1;
    end
end

end

