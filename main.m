clear;
%% Main file for Computer Vision Project 3

% Steps

% 1 - derive each image's mean, variance, skewness, kurtosis
% 2 - find the highest and lowest mean, variance, skewness, and kurtosis
% and store them for later use in Classification step
% 3 - normalize all feature vectors and store in a texture library
% 4 - divide each texture image in 100 blocks of size 64x64 (590 total?)
% 5 - compute a feature vector for all blocks
% 6 - 


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

% Feature vector array
V = zeros(59, 4, 8);

% Read in the textures, calculate the laplacian pyramid for each texture,
% calculate and enter statistics of each layer of the pyramid into a 
% feature vector, and finally, find the max/min value of statistics across 
% all feature vectors. 
for i = 1:59
    % read in the file
    D = sprintf('D%i.bmp', i);
    img = double(imread(D));

    % generate laplacian pyramid of image
    pyr = genPyr( img, 'lap', 8);
    for j = 1:8
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

%% Normalize Feature vectors

for i = 1:59
    for j = 1:8
        V(i,:,j) = [V(i,1,j), V(i,2,j), V(i,3,j), V(i,4)];
    end
end
