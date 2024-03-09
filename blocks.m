%% Divide texture into blocks and save to file
clear
load("img.mat")

% count 
count = 1;
for k = 1:59
    temp = double(img(:,:,k));
    for i = 0:64:639
        for j = 0:64:639
            % Store image to a block
            block(:,:,count) = temp(i+1:i+64, j+1:j+64);

            % Store the block label for testing. 
            block_label(count) = k;
            count = count + 1;  % increase count for next image to store

        end
    end
end


save("blocks.mat", "block", "block_label");
disp("blocks file generated");
clear