%% Main file for Computer Vision Project 3
clear;

% Load Laplacian Data
load("laplacian_4Layer.mat");
% Variables:
%   img = database of images in double format (640 x 640 x 59)
%   kur_min and kur_max = minimum and maximum kurtosis values
%   mean_min & mean_max = " "
%   skew_max & skew_min = " "
%   V = feature vectors of laplacian blocks
%   var_min & var_max = minimum and maximum variance values

% Load blocks to test
load("blocks.mat");
%   blocks = blocks to test in double format (64 x 64 x 5900)

% Correct stores the number of correct classifications
correct = 0;

% Test each block and classify it based on feature vectors
for i = 1:5900
    % feature vector calc.
    Bm = mean(block(:,:,i), "all");           % calculates the mean
    Bv = var(block(:,:,i), 0, "all");         % calculates the variance
    Bs = skewness(block(:,:,i), 1, "all");    % calculates the skewness
    Bk = kurtosis(block(:,:,i), 1, "all");    % calculates the kurtosis
    
    % Normalized Vector using my custom-made function, val_norm.
    U(i,:) = [val_norm(Bm, mean_min, mean_max), val_norm(Bv, var_min, ...
        var_max), val_norm(Bs, skew_min, skew_max), val_norm(Bk, ...
        kur_min, kur_max)];

    % label will store the estimated label. Row corresponds to block while
    % column corresponds to layer. In other words, row 2 column 1
    % corresponds to block 2, layer 1 estimation. 
    label = []; 
    % Calculate Euclidean Distance using laplace layer 1 of texture and 
    % label it. 
    best_L1 = [];

    % Loop through textures
    for j = 1:59
        % Euclidean Distance calculation using layer 1 row j (i.e. texture
        % j since rows correspond to textures).
        dist = sqrt(sum((V(j,:,1)-U(i,:)) .^2));
        if j == 1
            best_L1 = dist;
            label(i,1) = j;
        end
        if (dist < best_L1)
            best_L1 = dist;
            label(i,1) = j;
        end
    end
%     
%     % Calculate Euclidean Distance using laplace layer 2 of texture and 
%     % label it. 
%     best_L2 = [];
%     for j = 1:59
%         dist = sqrt(sum((V-U) .^2));
%         if j == 1
%             best_L2 = dist;
%         end
%         if (dist < best_L2)
%             best_L2 = dist;
%             label(i,2) = j;
%         end
%     end
% 
%     % Calculate Euclidean Distance using laplace layer 3 of texture and 
%     % label it. 
%     best_L3 = [];
%     for j = 1:59
%         dist = sqrt(sum((V-U) .^2));
%         if j == 1
%             best_L3 = dist;
%         end
%         if (dist < best_L3)
%             best_L3 = dist;
%             label(i,3) = j;
%         end
%         
%     end
% 
%     % Calculate Euclidean Distance using laplace layer 4 of texture and 
%     % label it. 
%     best_L4 = [];
%     for j = 1:59
%         dist = sqrt(sum((V-U) .^2));
%         if j == 1
%             best_L4 = dist;
%         end
%         if (dist < best_L4)
%             best_L4 = dist;
%             label(i,4) = j;
%         end
%     end


    % Check if classification is correct for layer 1; increment if true.
    if (block_label(i) == label(i,1))
        correct = correct + 1;
    end
end

PCC_L1 = (correct ./ 5900) .* 100;
disp(PCC_L1);