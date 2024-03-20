%% Main file for Computer Vision Project 3
clear;

% Add pyramid toolbox to path
addpath(genpath('.\Laplacian Pyramid Toolbox'));

% Load images for some figure generation
load("img.mat");

% Load blocks to test
load("blocks.mat");
%   blocks = blocks to test in double format (64 x 64 x 5900)
%   U is the normalized feature vector of the Laplacian
%   block_label = stores the label of the block for classification

%% Laplacian Classification
run = 0;
for pyr_num = 4
run = run+1;
fprintf("Number of Layers = %i", pyr_num);
% Load Laplacian Data
D = sprintf("Laplacian_%i.mat", pyr_num);
load(D);
% Variables:
%   img = database of images in double format (640 x 640 x 59)
%   kur_min and kur_max = minimum and maximum kurtosis values
%   mean_min & mean_max = " "
%   skew_max & skew_min = " "
%   V = feature vectors of laplacian blocks
%   var_min & var_max = minimum and maximum variance values

% Correct stores the number of correct classifications
correct = 0;
label = []; 

% Stores mislabeled textures; Column will be the correct label, row will
% be what it was mislabeled to. The value at row x column is how many times
% it was mislabeled to this label. 
mislabel = zeros(59, 59);
mislabel_G = zeros(59, 59);

% Test each block and classify it based on feature vectors
for i = 1:5900

    % Calculate Euclidean Distance using laplace layer 1 of texture and 
    % label it. 
    best_L1 = [];

    % Loop through textures
    for j = 1:59
        % Euclidean Distance calculation
        % Pick which statistic to use with row being the texture, column
        % being the statistic, and dim3 being the layer. Layer 1 is the
        % smallest and layer 5 is the largest. 
        v{1} = (V(j,1,1) - U(i,1,1)).^2;
        v{2} = (V(j,2,:) - U(i,2,:)).^2;
        v{3} = (V(j,3,2:4) - U(i,3,2:4)).^2;
        v{4} = (V(j,4,2:4) - U(i,4,2:4)).^2;
        temp = [v{1}(:); v{2}(:); v{3}(:); v{4}(:)];
        temp2 = sum(temp);
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
    
    else
        corr_texture = block_label(i);
        mislabel(label(i), corr_texture) = mislabel(label(i), corr_texture) +1;
    end
end
% Percent Correctly Classified Laplacian
PCC_L = (correct ./ 5900) .* 100;
disp(PCC_L);
if run == 1
    Laplacian_optimal = round(PCC_L, 2);
    Layer_optimal = pyr_num;
end
if PCC_L > Laplacian_optimal
    Laplacian_optimal = round(PCC_L, 2);
    Layer_optimal = pyr_num;
end
end
D = sprintf("Best Laplacian PCC is %f%% accurate with %i Layers", ...
    Laplacian_optimal, Layer_optimal);
disp(D);
[mislabel Ind] = max(mislabel);


%% Gabor Classification
run = 0;
% Run through every number of scales and orientations
for nscale = 4
for norient = 4
run = run + 1;
fprintf("scale = %i; orientation = %i", nscale, norient);

%Load Gabor Data file for a given scale and orientation
D = sprintf('Gabor_%i_%i.mat', nscale, norient); % String representing the filename
load(D);

% Correct stores the number of correct classifications with Gabor
correct_G = 0;

% Test each block and classify it based on feature vectors
total = 5900;
for i = 1:total

    label = []; 
    best_L1 = [];

    % Loop through textures and check the euclidean distance of each
    % texture with the block. Assign the label of the block to the texture
    % with the closest euclidean distance. 
    for j = 1:59 % Loop through textures 1 to 59

        % Euclidean Distance calculation
        % Pick which statistic to use with each cell being a different 
        % texture. Each column corresponds to a different statistic while
        % the number of rows are determined by the number of scales x
        % orientations (that is, a 4 scale 4 orientation filter has 16 rows
        % with each containing a set of statistics)
        v{1} = (V_Gab{j}(:,1) - U_Gab{i}(:,1)).^2; %mean
        v{2} = (V_Gab{j}(:,2) - U_Gab{i}(:,2)).^2; %var
        %v{3} = (V_Gab{j}(:,3) - U_Gab{i}(:,3)).^2; %skew
        % v{4} = (V_Gab{j}(:,4) - U_Gab{i}(:,4)).^2; %kur
        temp = [v{1}(:); v{2}(:)]; %; v{3}(:); v{4}(:)];
        temp2 = sum(temp);
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


    % Check if classification is correct; increment correct_G if true.
    if (block_label(i) == label(i))
        correct_G = correct_G + 1;
    else
        corr_texture = block_label(i);
        mislabel_G(label(i), corr_texture) = mislabel_G(label(i), corr_texture) +1;
    end
end

% Percent Correctly Classified Gabor
PCC_G = (correct_G ./ total) .* 100;
disp(PCC_G);
if run == 1
    Gabor_optimal = round(PCC_G, 2);
    norient_optimal = norient;
    nscale_optimal = nscale;
end
if PCC_G > Gabor_optimal
    Gabor_optimal = round(PCC_G, 2);
    norient_optimal = norient;
    nscale_optimal = nscale;
end
end
end
[mislabel_G Ind_G] = max(mislabel_G);

D = sprintf("Best Gabor PCC is %f%% accurate with %i scales and %i orientations", ...
    Gabor_optimal, nscale_optimal, norient_optimal);
disp(D);

%% Plot an image of laplacian pyramid of a texture
figure(Color="white");
subplot(2,2,1);
pyr = lpd((img(:,:,32)/256), 'Burt', 4);

% Function I made to convert pyramid to image. Must be 4 layers at any size.
toshow = pyr_to_img4(pyr);
imshow(toshow);
title("Texture 32 pyramid");

% arbitrarily pick a block from within texture 51
pyr_bl = lpd(img(:,:,34)/256, 'Burt', 4);
toshow_bl = pyr_to_img4(pyr_bl);

subplot(2,2,2);
imshow(toshow_bl);
title("Texture 34 pyramid");

subplot(2,2,3);
imshow(mat2gray(img(:,:,32)));
title("Texture 32")

subplot(2,2,4);
imshow(mat2gray(img(:,:,34)));
title("Texture 34");

figure(Color="White");
minWaveLength = 3;
mult = 2;
sigmaOnf = 0.65;
dThetaOnSigma = 1.5;
nscale = 4;
norient = 4;
E0 = gaborconvolve(img(:,:,15),nscale, norient, minWaveLength, mult, ...
        sigmaOnf, dThetaOnSigma);

[X,Y]=size(img(:,:,15));		% Size of the texture image
for i=1:nscale
    for j=1:norient
        ind=(i-1)*norient+j;      		% Calculate the index of each sub-plot
        subplot(4,4,ind);           		% Create a multi-figure plot
        Mi=abs(E0{i,j});            		% Create the magnitude for each Gabor channel
        imshow(mat2gray(Mi));              		% Show the Gabor filter output
        title(sprintf("S=%i; O=%i", i, j));
        Miv{ind}=reshape(Mi,X*Y,1); 	% Reshape the matrix data to vector data
    end
end


figure(Color="White");
E1 = gaborconvolve(img(:,:,14),nscale, norient, minWaveLength, mult, ...
        sigmaOnf, dThetaOnSigma);

[X,Y]=size(img(:,:,14));		% Size of the texture image
for i=1:nscale
    for j=1:norient
        ind=(i-1)*norient+j;      		% Calculate the index of each sub-plot
        subplot(4,4,ind);           		% Create a multi-figure plot
        Mi=abs(E1{i,j});            		% Create the magnitude for each Gabor channel
        imshow(mat2gray(Mi));              		% Show the Gabor filter output
        title(sprintf("S=%i; O=%i", i, j));
        Miv{ind}=reshape(Mi,X*Y,1); 	% Reshape the matrix data to vector data
    end
end



