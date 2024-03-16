%% This file generates the feature vectors for the Laplacian Pyramid
% The data is stored in a .mat file to save time when running other code
clear;
load("blocks.mat");
load("img.mat");
% Variables to store max values
mean_max = [];
var_max = [];
skew_max = [];
kur_max = [];

% Variables to store min values
mean_min = [];
var_min = [];
skew_min = [];
kur_min = [];

% Read in the textures, calculate the laplacian pyramid for each texture,
% calculate and enter statistics of each layer of the pyramid into a 
% feature vector, and finally, find the max/min value of statistics across 
% all feature vectors. 

pyr_num = 5; % Number of Pyramid layers to use
for i = 1:59
    % For this Pyramid toolbox, the image must be normalized by 256
    img(:,:,i) = img(:,:,i)/256; 

    % generate laplacian pyramid of image
    pyr = lpd(img(:,:,i), 'Burt', pyr_num);

    % For each layer...
    for j = 1:pyr_num
        % feature vector calc.
        Dm = mean(pyr{j}(:));        % calculates the mean
        Dv = var(pyr{j}(:));         % calculates the variance
        Ds = skewness(pyr{j}(:));    % calculates the skewness
        Dk = kurtosis(pyr{j}(:));    % calculates the kurtosis
        
        % Store values in feature vector
        V(i,:,j) = [Dm, Dv, Ds, Dk];

        % set initial max and min values
        if (i == 1)
            % initialize starting value for each layer
            % j corresponds to the layer this value is for.
            mean_max(j) = Dm;
            mean_min(j) = Dm;
            var_max(j) = Dv;
            var_min(j) = Dv;
            skew_max(j) = Ds;
            skew_min(j) = Ds;
            kur_max(j) = Dk;
            kur_min(j) = Dk;
            
        end

        % Update max/min values every time a new image is read in.
        if (Dm> mean_max(j))
            mean_max(j) = Dm;
        end
        if (Dv > var_max(j))
            var_max(j) = Dv;
        end
        if (Ds > skew_max(j))
            skew_max(j) = Ds;
        end
        if (Dk > kur_max(j))
            kur_max(j) = Dk;
        end

        % Update min values
        if (Dm < mean_min(j))
            mean_min(j) = Dm;
        end
        if (Dv < var_min(j))
            var_min(j) = Dv;
        end
        if (Ds < skew_min(j))
            skew_min(j) = Ds;
        end
        if (Dk < kur_min(j))
            kur_min(j) = Dk;
        end
    end
end


%% Normalize the feature vectors
for i = 1:59
    for j = 1:pyr_num
        % val_norm is a function I made to normalize values
        m = val_norm(V(i,1,j), mean_min(j), mean_max(j));  %norm mean
        v = val_norm(V(i,2,j), var_min(j), var_max(j));%norm variance
        s = val_norm(V(i,3,j), skew_min(j), skew_max(j));  %norm skew 
        k = val_norm(V(i,4,j), kur_min(j), kur_max(j));%norm kurtosis

        % Replace the vector with the normalized vector
        V(i,:,j) = [m, v, s, k];
    end
end


%% Calculate statistics of block files
for i=1:5900
    % generate laplacian pyramid of the block
    block(:,:,i) = block(:,:,i)/256;
    pyr_block = lpd(block(:,:,i), 'Burt', pyr_num);

    for j = 1:pyr_num % number of pyramid layers used
        % feature vector calc.
        Dm = mean(pyr_block{j}, "all");           % calculates the mean
        Dv = var(pyr_block{j}, 0, "all");         % calculates the variance
        Ds = skewness(pyr_block{j}, 1, "all");    % calculates the skewness
        Dk = kurtosis(pyr_block{j}, 1, "all");    % calculates the kurtosis
        
        % val_norm is a function I made to normalize values
        Dm = val_norm(Dm, mean_min(j), mean_max(j));  %norm mean
        Dv = val_norm(Dv, var_min(j), var_max(j));%norm variance
        Ds = val_norm(Ds, skew_min(j), skew_max(j));  %norm skew 
        Dk = val_norm(Dk, kur_min(j), kur_max(j));%norm kurtosis

        % Store Normalized values in vector
        U(i,:,j) = [Dm, Dv, Ds, Dk];
    end
end

%% Save other data to laplacian data file
% Save remaining variables to this file
save("laplacian_4Layer.mat", "V", "U");

% Clear remaining variables from workspace
clear
disp("Laplacian file Generated");