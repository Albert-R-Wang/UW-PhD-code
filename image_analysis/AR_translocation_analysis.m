% ------------------------------------------------------------------------------
% Title:         Nucleus-to-Cytoplasm Signal Quantification Analysis
% Author:        Albert Wang
% Last updated:  2025-07-17
% ------------------------------------------------------------------------------

%% Summary
% Androgen receptor (AR) is typically located in the cytoplasm in its inactive state 
% and translocates to the nucleus upon binding to androgens, 
% where it functions as a transcription factor.
% Quantifying AR nuclear translocation is essential for assessing androgen receptor activation, 
% as its movement from the cytoplasm to the nucleus is a key step in regulating gene expression 
% involved in development, differentiation, and disease progression, particularly in prostate cancer.
% 
% This script reads cytoplasm, nucleus, and AR channel TIFF images selected by the user,
% converts each RGB image into a summed grayscale intensity map,
% creates binary masks for the cytoplasm and nucleus regions via user-defined thresholding,
% computes the AR signal sum and mean intensity in each region,
% calculates both raw and intensity-normalized nucleus-to-cytoplasm (NC) ratios,
% outlines and visualizes the segmented regions on original images,
% counts nuclei while excluding small objects based on area,
% and exports the AR segmentation figure with overlays to a PNG file.
%
% Example usage: https://doi.org/10.1016/j.gendis.2023.07.001
% Wang, A. R., Baschnagel, A. M., Ni, Z., Brennan, S. R., Newton, H. K., Buehler, D., . . . Iyer, G. (2024). 
% Network analyses: Inhibition of androgen receptor signaling reduces inflammation 
% in the lung through AR-MAF-IL6 signaling axes. Genes & Diseases, 11(3), 101072.


%% Script Outline
% 1. Data Import and Intensity Conversion
% 2. Thresholding and Mask Generation
% 3. Segmentation Visualization
% 4. AR Signal Quantification and NC Ratio Calculation
% 5. Integrated Density and Normalized NC Ratio
% 6. Nucleus Counting (Excluding Small Objects)
% 7. Result Display and Export

%% =============================================================================
% Required toolbox: Image Processing Toolbox

%% =============================================================================
% 1. Data Import and Intensity Conversion

close all;
clear all;
clc;

% Select input images (cytoplasm, AR, nucleus)
[files, pathname, ~] = uigetfile('*.tif', 'Select all image files', 'MultiSelect', 'on');
cyto = imread(files{1});
AR    = imread(files{2});
nuc   = imread(files{3});


% Convert each image to grayscale by summing RGB channels
A = {cyto, nuc, AR};
I = cell(1,3); % Grayscale images

for n = 1:3
    rmat = A{n}(:,:,1);
    gmat = A{n}(:,:,2);
    bmat = A{n}(:,:,3);
    
    figure(n);
    subplot(2,2,1), imshow(rmat); title('Red Channel');
    subplot(2,2,2), imshow(gmat); title('Green Channel');
    subplot(2,2,3), imshow(bmat); title('Blue Channel');
    subplot(2,2,4), imshow(A{n});  title('Original Image');
    
    I{n} = rmat + gmat + bmat;
end

% Display all raw and grayscale images side by side
figure(4);
subplot(2,3,1), imshow(cyto); title('Cytoplasm original image');
subplot(2,3,2), imshow(nuc);  title('Nucleus original image');
subplot(2,3,3), imshow(AR);   title('AR original image');
subplot(2,3,4), imshow(I{1}); title('Cytoplasm');
subplot(2,3,5), imshow(I{2}); title('Nucleus');
subplot(2,3,6), imshow(I{3}); title('AR');

% Display file names selected
disp(files{1});
disp(files{2});
disp(files{3});

%% =============================================================================
% 2. Thresholding and Mask Generation

% Create binary masks for cytoplasm and nucleus
mask   = cell(1,2);
binary = cell(1,2);

cyto_threshold = input('Enter threshold for cytoplasm (0-1): '); % default: 0.2
nuc_threshold  = input('Enter threshold for nucleus (0-1): ');
threshold = [cyto_threshold, nuc_threshold];

for n = 1:2
    bw = imbinarize(I{n}, threshold(n));     % Apply threshold
    bw = bwareaopen(bw, 100);                 % Remove small objects
    BWfill = imfill(bw, 'holes');            % Fill holes
    
    figure; 
    subplot(1,2,1), imshow(bw); title('Binary Mask')
    subplot(1,2,2), imshow(BWfill); title('Filled Binary Mask');
    
    binary{n} = BWfill;
    mask{n} = find(BWfill);                  % Store linear indices
end

%% =============================================================================
% 3. Segmentation Visualization

% Outline segmentation boundaries
outline_cyto = bwperim(binary{1});
outline_nuc  = bwperim(binary{2});

% Overlay outlines on cytoplasm and nucleus images
Segout_cyto = I{1}; Segout_cyto(outline_cyto) = 255;
Segout_nuc  = I{2}; Segout_nuc(outline_nuc)  = 255;

figure;
subplot(1,2,1), imshow(Segout_cyto); title('Cytoplasm Segmentation');
subplot(1,2,2), imshow(Segout_nuc);  title('Nucleus Segmentation');

% Overlay outlines on AR image
Segout_ARcyto = I{3}; Segout_ARcyto(outline_cyto) = 255;
Segout_ARnuc  = I{3}; Segout_ARnuc(outline_nuc)  = 255;

figure;
subplot(1,2,1), imshow(Segout_ARcyto); title('AR - Cytoplasm Segmentation');
subplot(1,2,2), imshow(Segout_ARnuc);  title('AR - Nucleus Segmentation');
%% =============================================================================
% 4. AR Signal Quantification and NC Ratio Calculation
% 
% This section extracts AR signal intensity from the cytoplasm and nucleus
% using the previously generated binary masks. It calculates the total AR signal 
% (summed pixel intensities) within each compartment, then computes the 
% nucleus-to-cytoplasm (NC) ratio based on these sums. This total NC ratio 
% reflects the overall distribution (ratio) of AR protein between nucleus and cytoplasm.

% Extract AR signal in masked regions
ARsignal_cyto = I{3}(mask{1});
ARsignal_nuc  = I{3}(mask{2});
ARsum_cyto = sum(ARsignal_cyto);
ARsum_nuc  = sum(ARsignal_nuc);
NCratio = ARsum_nuc / ARsum_cyto;



%% =============================================================================
% 5. Integrated Density and Normalized NC Ratio
% 
% This section calculates the average AR intensity (mean pixel value) within
% the cytoplasm and nucleus by dividing total AR signal by the number of pixels 
% in each region. It then computes a normalized NC ratio using these average 
% values. This area-normalized metric provides a more accurate comparison 
% of AR concentration between compartments, particularly when the nucleus and 
% cytoplasm differ in size or cell morphology varies across images.


% Calculate area (number of pixels in each mask)
area = [length(mask{1}), length(mask{2})];
int_cyto = ARsum_cyto / area(1);
int_nuc  = ARsum_nuc  / area(2);
NCratio_int = int_nuc / int_cyto;

%% =============================================================================
% 6. Nucleus Counting (Excluding Small Objects)

% Label connected components in nucleus mask
cc = bwconncomp(binary{2}, 26);
S = regionprops(cc);
number = cc.NumObjects;

% Exclude small nuclei (area < 300 px)
numOfSmall = sum([S.Area] < 300);
numOfObject = number - numOfSmall;

%% =============================================================================
% 7. Result Display and Export

% Print NC ratio results
fprintf('\n=== Results ===\n');
fprintf('Nucleus-to-Cytoplasm AR Signal Ratio (Total): %.4f\n', NCratio);
fprintf('Nucleus-to-Cytoplasm AR Signal Ratio (Intensity-normalized): %.4f\n', NCratio_int);
fprintf('Number of nuclei (area ≥ 300 px): %d\n', numOfObject);

% Prepare filename for figure export
[~, name, ~] = fileparts(files{1});
name_clean = regexprep(name, '_ch\d{2}$', ''); % Remove _ch00-like suffix
exportName = fullfile(pathname, [name_clean '_AR_segmentation.png']);

% Export final AR segmentation figure
saveas(gcf, exportName);
fprintf('Figure exported to: %s\n', exportName);


