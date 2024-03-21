%% Generates a laplacian pyramid of Notre Dame
clear

% Add pyramid toolbox to path
addpath(genpath('.\Laplacian Pyramid Toolbox'));

% image under creative commons license and sourced from:
% https://www.flickr.com/photos/michalo/6094164096
img = double(im2gray(imread("Notre_Dame.jpg")));

% Plot an image of laplacian pyramid of a texture
figure(Color="white");

subplot(1,2,1);
imshow("Notre_Dame.jpg");

subplot(1,2,2);
pyr = lpd((img/256), 'Burt', 4);

% Function I made to convert pyramid to image. Must be 4 layers at any size.
toshow = pyr_to_img4(pyr);
imshow(toshow);
title("Notre Dame Pyramid");