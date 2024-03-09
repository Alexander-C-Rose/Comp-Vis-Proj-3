%% This file generates the feature vectors for the Laplacian Pyramid
% The data is stored in a .mat file to save time when running other code
clear;

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

pyr_num = 4; % Number of Pyramid layers to use
for i = 1:59
    % read in the file
    D = sprintf('D%i.bmp', i);
    img(:,:,i) = double(imread(D));

    % generate laplacian pyramid of image
    pyr = genPyr( img, 'lap', pyr_num);
    for j = 1:pyr_num
        % feature vector calc.
        Dm = mean(pyr{j}, "all");           % calculates the mean
        Dv = var(pyr{j}, 0, "all");         % calculates the variance
        Ds = skewness(pyr{j}, 1, "all");    % calculates the skewness
        Dk = kurtosis(pyr{j}, 1, "all");    % calculates the kurtosis
        
        % Store values in feature vector
        V(i,:,j) = [Dm, Dv, Ds, Dk];

        % set initial max and min values
        if ((i == 1) && (j == 1))
            mean_max = Dm;
            mean_min = Dm;
            var_max = Dv;
            var_min = Dv;
            skew_max = Ds;
            skew_min = Ds;
            kur_max = Dk;
            kur_min = Dk;
        end

        % Update max values
        if (Dm > mean_max)
            mean_max = Dm;
        end
        if (Dv > var_max)
            var_max = Dv;
        end
        if (Ds > skew_max)
            skew_max = Ds;
        end
        if (Dk > kur_max)
            kur_max = Dk;
        end

        % Update min values
        if (Dm < mean_min)
            mean_min = Dm;
        end
        if (Dv < var_min)
            var_min = Dv;
        end
        if (Ds < skew_min)
            skew_min = Ds;
        end
        if (Dk < kur_min)
            kur_min = Dk;
        end
    end
end

%% Normalize the feature vectors
for i = 1:59
    for j = 1:pyr_num
        % val_norm is a function I made to normalize values
        mean = val_norm(V(i,1,j), mean_min, mean_max);  %norm mean
        variance = val_norm(V(i,2,j), var_min, var_max);%norm variance
        skew = val_norm(V(i,3,j), skew_min, skew_max);  %norm skew 
        kurtosis = val_norm(V(i,4,j), kur_min, kur_max);%norm kurtosis

        % Replace the vector with the normalized vector
        V(i,:,j) = [mean, variance, skew, kurtosis];    
    end
end

% Save image data to separate file (used by multiple programs)
save("img.mat", "img");   

% Clear unneeded variables from memory before saving
clear D Dk Dm Ds Dv pyr i j pyr_num mean variance skew kurtosis img

% Save remaining variables to this file
save("laplacian_4Layer.mat");

% Clear remaining variables from workspace
clear
disp("Laplacian file Generated");