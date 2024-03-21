# Comp-Vis-Proj-3

project_report.pdf and .pptx contain the same information. 

blocks.m reads and converts images into usable formats for the files, storing the data in the blocks.mat file and img.mat file.

main.m tries to classify the textures and determine the optimal number of layers for the Laplacian and scales/orientations for the Gabor, importing data from the many .mat files.

Laplacian.m makes the laplacian pyramid data for each texture and the blocks and stores to several .mat files. 

Gabor.m performs the gabor filtering on the textures and blocks and stores the data for a wide range of scales and orientation configurations in .mat files.

