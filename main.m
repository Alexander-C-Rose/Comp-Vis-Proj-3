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
%   U is the normalized feature vector of the Laplacian
%   block_label = stores the label of the block for classification

%% Laplacian Classification

% Correct stores the number of correct classifications
correct = 0;

% Test each block and classify it based on feature vectors
for i = 1:5900

    % label will store the estimated label. Row corresponds to block while
    % column corresponds to layer. In other words, row 2 column 1
    % corresponds to block 2, layer 1 estimation. 
    label = []; 
    % Calculate Euclidean Distance using laplace layer 1 of texture and 
    % label it. 
    best_L1 = [];

    % Loop through textures
    for j = 1:59
        % Euclidean Distance calculation
        temp = (V(j,2,2:5) - U(i,2,2:5)); % Use variance of layers 2 to 5
        temp2 = sum(temp .^2) + (V(j,1,1) - U(i,1,1)).^2; % mean of layer 1
        dist = sqrt(temp2);

        % Initialize label estimate
        if j == 1
            best_L1 = dist;
            label(i) = j;
        end

        % If the new distance is closer than the old, assign new label
        if (dist < best_L1)
            best_L1 = dist;
            label(i) = j;
        end
    end


    % Check if classification is correct; increment if true.
    if (block_label(i) == label(i))
        correct = correct + 1;
    end
end

% Percent Correctly Classified Laplacian
PCC_L = (correct ./ 5900) .* 100;
disp(PCC_L);

%% Gabor Classification

%Load Gabor Data
load("Gabor.mat");

% Correct stores the number of correct classifications with Gabor
correct_G = 0;

% Test each block and classify it based on feature vectors
total = 5900;
for i = 1:total

    label = []; 
    best_L1 = [];

    % Loop through textures
    for j = 1:59
        % Euclidean Distance calculation
        temp = (V_Gab(j,1) - U_Gab(i,1));
        temp2 = sum(temp .^2) + (V_Gab(j,1) - U_Gab(i,1)).^2;
        dist = sqrt(temp2);

        % Initialize label estimate
        if j == 1
            best_L1 = dist;
            label(i) = j;
        end

        % If the new distance is closer than the old, assign new label
        if (dist < best_L1)
            best_L1 = dist;
            label(i) = j;
        end
    end


    % Check if classification is correct; increment if true.
    if (block_label(i) == label(i))
        correct_G = correct_G + 1;
    end
end

% Percent Correctly Classified Laplacian
PCC_G = (correct_G ./ total) .* 100;
disp(PCC_G);
