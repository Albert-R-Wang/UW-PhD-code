% ------------------------------------------------------------------------------
% Title:         Nuclear Foci Quantification Pipeline
% Author:        Albert Wang
% Last updated:  2025-07-17
% ------------------------------------------------------------------------------

%% Summary
% γH2AX is the phosphorylated form of the histone variant H2AX, rapidly formed 
% at sites of DNA double-strand breaks as an early cellular response to genotoxic stress. 
% Quantifying γH2AX foci provides a sensitive and reliable indicator of DNA 
% double-strand breaks, allowing researchers to assess both the extent of 
% DNA damage and the cellular response to genotoxic stress, such as 
% radiation therapy.
% 
% This script reads nucleus and foci TIFF images from user selection,
% converts each to a summed-intensity grayscale map,
% applies either dynamic or universal thresholding,
% filters out small nuclei and removes their associated foci,
% computes per-nucleus intensity statistics and foci counts,
% writes results to an Excel file,
% and saves an annotated segmentation JPEG with colored nucleus indices.
%
% Example usage: https://doi.org/10.1158/1535-7163.MCT-17-0897
% Elsaid, M., Shahi De, A., Wang, A., Baiu, D., Li, C., Werner, L., . . . Otto, M. (2018). 
% Enhanced Radiosensitivity in Solid Tumors using a Tumor-Selective Alkyl Phosphocholine Analog. 
% Mol Cancer Ther, 17, molcanther.0897.2017.


%% Script Outline
% 1. Data Import and Intensity Conversion
% 2. Preview Images
% 3. Threshold Selection
% 4. Nucleus Masking
% 5. Foci Masking and Filtering
% 6. Per-Cell Metric Computation
% 7. Export Results to Excel
% 8. Segmentation Visualization and JPEG Export

%% =============================================================================
% Required toolbox: Image Processing Toolbox

%% =============================================================================
% 1. Data Import and Intensity Conversion

close all;        % Close all open figure windows
clear all;        % Clear workspace variables
clc;              % Clear command window output

% Prompt user to select two TIFF files: first nucleus, then foci
[files, pathname] = uigetfile('*.tif', ...
    'Select nucleus then foci images', 'MultiSelect', 'on');
% Read images into arrays
nuc  = imread(fullfile(pathname, files{1}));
foci = imread(fullfile(pathname, files{2}));

% Combine RGB channels to get summed-intensity maps for segmentation
A = {nuc, foci};
I = cell(1, 2);
for n = 1:2
    r = A{n}(:,:,1);  % Red channel
    g = A{n}(:,:,2);  % Green channel
    b = A{n}(:,:,3);  % Blue channel
    % Sum channels and convert to double for processing
    I{n} = double(r) + double(g) + double(b);
end

%% =============================================================================
% 2. Preview Images
% Display original and intensity-converted images for user verification
figure;
subplot(2,2,1), imshow(nuc),    title('Nucleus original');
subplot(2,2,2), imshow(foci),   title('Foci original');
subplot(2,2,3), imshow(I{1},[]),title('Nucleus intensity');
subplot(2,2,4), imshow(I{2},[]),title('Foci intensity');

%% =============================================================================
% 3. Threshold Selection
% Prompt user for thresholds
nuc_thresh  = input('Enter nucleus threshold (0–255): ');    % Fraction for imbinarize
foci_thresh = input('Enter foci threshold (0–255): ');    % Intensity cutoff

% The pipeline supports two thresholding modes: dynamic mode applies nucleus-specific thresholds 
% to account for local intensity variation, while universal mode uses a global threshold for all nuclei 
% to ensure consistency across the dataset.
mode_sel    = input('Threshold mode: dynamic (1) or universal (2)? ');

%% =============================================================================
% 4. Nucleus Masking
% Create binary mask of nuclei, remove small noise, fill holes
bw_nuc = imbinarize(I{1}, nuc_thresh);   % Threshold nucleus intensity
bw_nuc = bwareaopen(bw_nuc, 50);        % Remove objects smaller than 50 pixels
bw_nuc = imfill(bw_nuc, 'holes');       % Fill any holes inside nuclei

figure;
subplot(1,2,1), imshow(I{1}, []), title('Original necleus image')
subplot(1,2,2), imshow(bw_nuc), title('Nucleus mask')


% Label connected components in the nucleus mask
cc_nuc      = bwconncomp(bw_nuc, 26);
% Extract area and centroid of each nucleus
S_nuc       = regionprops(cc_nuc, 'Area', 'Centroid');
% Create a label matrix for per-pixel nucleus index
L_nuc       = labelmatrix(cc_nuc);

% Define filtering parameters
area_thresh = 1000;   % Minimum area (pixels) to keep a nucleus
area_cover  = 90;     % Half-size of square region around small nuclei to clear foci
% Identify small (to discard) and kept nuclei by index
small_idx   = find([S_nuc.Area] < area_thresh);
keep_idx    = find([S_nuc.Area] >= area_thresh);

%% =============================================================================
% 5. Foci Masking and Filtering
% Build a final foci mask by per-nucleus threshold
bw_foci_raw = (I{2} > foci_thresh);
bw_foci     = false(size(I{2}));

for i = keep_idx
    mask_i  = (L_nuc == i);
    pixVals = I{2}(mask_i);

    if mode_sel == 1
        thr_i = mean(pixVals)*2;
        thr_i = max(min(thr_i, 200), 50);
    else
        thr_i = foci_thresh;
    end

    foci_mask_i = mask_i & (I{2} > thr_i);
    foci_mask_i = bwareaopen(foci_mask_i, 2);
    bw_foci     = bw_foci | foci_mask_i;
end

% Remove any foci inside small/discarded nuclei regions
[xmax, ymax] = size(bw_foci);
for i = small_idx
    cx = round(S_nuc(i).Centroid(1));
    cy = round(S_nuc(i).Centroid(2));
    x1 = max(1, cx-area_cover); x2 = min(xmax, cx+area_cover);
    y1 = max(1, cy-area_cover); y2 = min(ymax, cy+area_cover);
    bw_foci(y1:y2, x1:x2) = false;
end

% Display raw vs final mask
figure;
subplot(1,3,1), imshow(I{2}, []), title('Original foci image')
subplot(1,3,2), imshow(bw_foci_raw), title('Initial foci mask');
subplot(1,3,3), imshow(bw_foci),     title('Final foci mask');

%% =============================================================================
% 6. Per-Cell Metric Computation
perCell = [];
count   = 0;

for i = keep_idx
    count   = count + 1;
    pixVals = I{2}(L_nuc == i);
    mu      = mean(pixVals);
    med     = median(pixVals);
    mo      = mode(pixVals);
    sd      = std(double(pixVals));

    cc_fi   = bwconncomp((L_nuc == i) & bw_foci, 8);
    stats_f = regionprops(cc_fi,'Area');
    fcount  = sum([stats_f.Area] >= 2);

    perCell(count,:) = [i, S_nuc(i).Area, mu, med, mo, sd, fcount];
end

%% =============================================================================
% 7. Export Results to Excel
% Convert to table with descriptive column names
T = array2table(perCell, ...
    'VariableNames',{'Index','Area','MeanInt','MedianInt','ModeInt','StdInt','FociCount'});

% Remove channel suffix (e.g., _ch00, _ch01, etc.) from the base file name
[~, baseName, ~] = fileparts(files{1});
baseNameClean = regexprep(baseName, '_ch\d+$', '');

% Construct output filename based on cleaned image name
outExcel = fullfile(pathname, [baseNameClean '_perCellResults.xlsx']);

% Write table to Excel file
writetable(T, outExcel);
disp(['Per-cell results written to: ' outExcel]);

%% =============================================================================
% 8. Segmentation Visualization and JPEG Export
% Create outlines of nuclei and filtered foci via perimeter extraction
outline_n = bwperim(bw_nuc);
outline_f = bwperim(bw_foci);
% Overlay outlines on intensity images
SegNuc = I{1}; SegNuc(outline_n) = max(SegNuc(:));
SegFoc = I{2}; SegFoc(outline_f) = max(SegFoc(:));

% Plot 2x2 panel with colored index map
figure;
subplot(6,6,[19,20,21,25,26,27,31,32,33]), imagesc(L_nuc),       title('Nucleus index');
colormap(colorcube);                   % Use vibrant colors per nucleus
axis square off;                       % Square aspect, hide axes
% Overlay nucleus index numbers at centroids
for a = keep_idx
    text(S_nuc(a).Centroid(1), S_nuc(a).Centroid(2), sprintf('%d',a), ...
         'Color','k','FontSize',6,'HorizontalAlignment','center', 'FontWeight','bold');
end

subplot(6,6,[1,2,3,7,8,9,13,14,15]), imshow(SegNuc,[]),   title('Nucleus segment');

subplot(6,6,[4,5,6,10,11,12,16,17,18]), imshow(SegFoc,[]),   title('Foci segment');

subplot(6,6,[22,23,24,28,29,30,34,35,36]), imshow(foci),         title('Original foci image');


% Ensure consistent font sizing across panels
set(findall(gcf,'-property','FontSize'),'FontSize',8);

% Save the composite figure as a high-resolution JPEG
outJpeg = fullfile(pathname, [baseNameClean '_segmentation.jpg']);
print(gcf, outJpeg, '-djpeg','-r600');
disp(['Segmentation JPEG saved to: ' outJpeg]);
