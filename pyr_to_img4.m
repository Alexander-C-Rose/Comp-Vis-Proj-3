function [toshow] = pyr_to_img4(pyr)
%UNTITLED2 Summary of this function goes here
[r, c] = size(pyr{4});
toshow = zeros(r,c + c/2);
toshow(1:r,1:c) = mat2gray(pyr{4});
toshow(1:r/2,(c+1):(c+c/2)) = mat2gray(pyr{3});
toshow((r/2+1):(r/2 + r/4),(c+1):(c+c/4)) = mat2gray(pyr{2});
toshow((r/2 + r/4 + 1):(r/2 + r/4 + r/8),(c+1):(c+c/8)) = mat2gray(pyr{1});

end