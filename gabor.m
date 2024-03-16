%% This file generates the feature vectors for the Laplacian Pyramid
% The data is stored in a .mat file to save time when running other code
clear;

% Load in file data
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

% Read in the textures, calculate the Gabor filter for each texture,
% calculate and enter statistics of each filter result into a 
% feature vector, and finally, find the max/min value of statistics across 
% all filtered images.

%   Variable       Suggested   Description
%   name           value
%  ----------------------------------------------------------
%    im                        Image to be convolved.
%    nscale          = 4;      Number of wavelet scales.
%    norient         = 6;      Number of filter orientations.
%    minWaveLength   = 3;      Wavelength of smallest scale filter.
%    mult            = 2;      Scaling factor between successive filters.
%    sigmaOnf        = 0.65;   Ratio of the standard deviation of the
%                              Gaussian describing the log Gabor filter's transfer function 
%	                       in the frequency domain to the filter center frequency.
%    dThetaOnSigma   = 1.5;    Ratio of angular interval between filter orientations
%			       and the standard deviation of the angular Gaussian
%			       function used to construct filters in the
%                              freq. plane.
%
%        Hence:
%        abs(EO{s,o}) returns the magnitude of the convolution over the
%                     image at scale s and orientation o.
%        angle(EO{s,o}) returns the phase angles.

nscale = 5;
norient = 6;
minWaveLength = 3;
mult = 2;
sigmaOnf = 0.65;
dThetaOnSigma = 1.5;

for i = 1:59
    % Filter the texture with Gabor filtering
    E0 = gaborconvolve(img(:,:,i),nscale, norient, minWaveLength, mult, ...
        sigmaOnf, dThetaOnSigma);

    % Create magnitude feature vector
    Mag = [];
    for s = 1:nscale
        for o = 1:norient
            Mag_temp = abs(E0{s,o});
            Mag = [Mag; Mag_temp];
        end
    end

    % statistics calc.
    Dm = mean(Mag(:));        % calculates the mean
    Dv = var(Mag(:));         % calculates the variance
    Ds = skewness(Mag(:));    % calculates the skewness
    Dk = kurtosis(Mag(:));    % calculates the kurtosis
    
    % Store values in feature vector
    V(i,:) = [Dm, Dv, Ds, Dk];

    % set initial max and min values
    if (i == 1)
        % initialize starting value for each layer
        % j corresponds to the layer this value is for.
        mean_max = Dm;
        mean_min = Dm;
        var_max = Dv;
        var_min = Dv;
        skew_max = Ds;
        skew_min = Ds;
        kur_max = Dk;
        kur_min = Dk;
        
    end

    % Update max/min values every time a new image is read in.
    if (Dm> mean_max)
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


%% Normalize the feature vectors
for i = 1:59
    % val_norm is a function I made to normalize values
    m = val_norm(V(i,1), mean_min, mean_max);  %norm mean
    v = val_norm(V(i,2), var_min, var_max);%norm variance
    s = val_norm(V(i,3), skew_min, skew_max);  %norm skew 
    k = val_norm(V(i,4), kur_min, kur_max);%norm kurtosis

    % Replace the vector with the normalized vector
    V(i,:) = [m, v, s, k];
end


%% Load texture blocks and save Gabor information to file
% I'm loading the blocks generated from the laplacian file but without the
% feature vector from said file. 

% Load blocks
load("blocks.mat", "block");
%   blocks = blocks to test in double format (64 x 64 x 5900)


for i=1:5900
    fprintf('Processing block %d \n', i);
    % generate laplacian pyramid of the block
    % Filter the texture with Gabor filtering
    E1 = gaborconvolve(block(:,:,i),nscale, norient, minWaveLength, mult, ...
        sigmaOnf, dThetaOnSigma);

    % Create magnitude feature vector
    Mag1 = [];
    for s = 1:nscale
        for o = 1:norient
            Mag_temp = abs(E0{s,o});
            Mag1 = [Mag1; Mag_temp];
        end
    end

    % statistics calc.
    Dm = mean(Mag1(:));           % calculates the mean
    Dv = var(Mag1(:));         % calculates the variance
    Ds = skewness(Mag1(:));    % calculates the skewness
    Dk = kurtosis(Mag1(:));    % calculates the kurtosis
    
    % val_norm is a function I made to normalize values
    Dm = val_norm(Dm, mean_min, mean_max);      %norm mean
    Dv = val_norm(Dv, var_min, var_max);        %norm variance
    Ds = val_norm(Ds, skew_min, skew_max);      %norm skew 
    Dk = val_norm(Dk, kur_min, kur_max);        %norm kurtosis

    % Store Normalized values in vector
    U(i,:) = [Dm, Dv, Ds, Dk];
end
% save("gabor.mat", "block", "block_label", "U");
% disp("blocks file generated");
% 
% %% Save other data to laplacian data file
% % Save image data to separate file (in case it's used by other programs)
% save("img.mat", "img");   
% 
% % Clear unneeded variables from memory before saving
% clear D Dk Dm Ds Dv E0 pyr_block i j pyr_num mean variance skew kurtosis img
% 
% % Save remaining variables to this file
% save("laplacian_4Layer.mat");
% 
% % Clear remaining variables from workspace
% clear
% disp("Laplacian file Generated");