%% This file generates the feature vectors for the Laplacian Pyramid
% The data is stored in a .mat file to save time when running other code
tic
clear;

% Load in file data
load("blocks.mat");
load("img.mat");

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

minWaveLength = 3;
mult = 2;
sigmaOnf = 0.65;
dThetaOnSigma = 1.5;
for nscale = 1:4
for norient = 1:6

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

% Feature vector cell
V_Gab = {};

% Increment textures and calculate feature vectors
for i = 1:59
    % Filter the texture with Gabor filtering
    E0 = gaborconvolve(img(:,:,i),nscale, norient, minWaveLength, mult, ...
        sigmaOnf, dThetaOnSigma);

    V_temp = [];
    V_temp2 = [];

    % Create magnitude feature vector
    Mag = [];
    for s = 1:nscale
        for o = 1:norient
            Mag = abs(E0{s,o});
            % statistics calc.
            Dm = mean(Mag(:));        % calculates the mean
            Dv = var(Mag(:));         % calculates the variance
            Ds = skewness(Mag(:));    % calculates the skewness
            Dk = kurtosis(Mag(:));    % calculates the kurtosis
            
            % Store values in feature vector
            V_temp = [Dm, Dv, Ds, Dk];
            V_temp2 = [V_temp2; V_temp]; 

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
    end

    % Assign the vector for this texture into a cell containing all
    % textures feature vectors
    V_Gab{i} = V_temp2;
end


%% Normalize the feature vectors
for i = 1:59
    [r, c] = size(V_Gab{i});
    for j = 1:r
        % val_norm is a function I made to normalize values
        m = val_norm(V_Gab{i}(j,1), mean_min, mean_max);  %norm mean
        v = val_norm(V_Gab{i}(j,2), var_min, var_max);%norm variance
        s = val_norm(V_Gab{i}(j,3), skew_min, skew_max);  %norm skew 
        k = val_norm(V_Gab{i}(j,4), kur_min, kur_max);%norm kurtosis
    
        % Replace the vector with the normalized vector
        V_Gab{i}(j,:) = [m, v, s, k];
    end
end


%% Load texture blocks and save Gabor information to file

% Load blocks
load("blocks.mat", "block");
%   blocks = blocks to test in double format (64 x 64 x 5900)

U_Gab = {};
% Run through blocks
for i=1:5900
    fprintf('Processing block %d \n', i);
    % generate laplacian pyramid of the block
    % Filter the texture with Gabor filtering
    E1 = gaborconvolve(block(:,:,i),nscale, norient, minWaveLength, mult, ...
        sigmaOnf, dThetaOnSigma);
    
    % Create magnitude feature vector
    Mag1 = [];
    U_temp = [];
    U_temp2 = [];
    for s = 1:nscale
        for o = 1:norient
            Mag1 = abs(E1{s,o});
            
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
        
            % Store values in feature vector
            U_temp = [Dm, Dv, Ds, Dk];
            U_temp2 = [U_temp2; U_temp];
        end
    end

    % Assign the vector for this texture into a cell containing all
    % textures feature vectors
    U_Gab{i} = U_temp2;
end
%% Save data to laplacian data file
% Save variables to a file
D = sprintf('Gabor_%i_%i.mat', nscale, norient); % String representing the filename
save(D, "V_Gab", "U_Gab");

% Clear remaining variables from workspace
disp("Gabor file Generated");
end
end

clear

toc