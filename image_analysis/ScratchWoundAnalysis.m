% ------------------------------------------------------------------------------
% Title:         Scratch Wound Assay Quantification
% Author:        Albert Wang
% Last updated:  2025-07-19
% ------------------------------------------------------------------------------


%% Summary
% The scratch wound assay is a widely used in vitro method for studying cell migration. 
% In this assay, a scratch or gap is created in a confluent monolayer of cells 
% using a pipette tip or other tools. Researchers then observe how cells move into the gap over time. 
% Under the microscope, the scratch typically appears as a clear, cell-free region 
% between two dense regions of cells, which gradually closes as the cells migrate inward.
%
% This script quantifies the area and average length of a manually defined scratch 
% wound from TIFF images. Prior to running the script, users must draw two white 
% boundary lines (see ScratchWoundAnalysis.pdf for example) that enclose the wound region. 
% These lines can be in vertical or horizontal orientation and are used to define the edges of the wound. 
% 
% The script processes each image by thresholding and skeletonizing the white 
% boundary lines, detecting their endpoints, identifying which endpoints belong 
% to the same boundary, and connecting them to form a closed region of interest (ROI). 
% It then fills the ROI and calculates the wound area in pixels. Additionally, 
% it computes the distance between the midpoints of the two connecting boundaries 
% as a representative measure of wound length, which can be used for normalization 
% or to enable consistent comparisons across images. Optional cropping can be applied 
% to remove unwanted rows (for vertical scratches) or columns (for horizontal scratches) 
% at the edges of each image.
%
% For each image, the script outputs an annotated version with the ROI highlighted 
% and the wound length line drawn in red. A summary Excel file is also generated, 
% compiling the area and length measurements across all processed images.


%% Script Outline
% 1. Image Selection and User Input
% 2. Image Preprocessing and Binary Line Extraction
% 3. Endpoint and Line Labeling
% 4. Connecting Boundary Endpoints and ROI Filling
% 5. Wound Area Calculation and Overlay Visualization
% 6. Wound length (Midpoint Distance) Measurement
% 7. Export of Annotated Image and Summary Table
% 8. Utility Function for Drawing Lines

%% =============================================================================
% Required toolbox: Image Processing Toolbox

%% =============================================================================
% 1. Image Selection and User Input

clear all
clc

% Select multiple image files
[files, pathname, filterindex] = uigetfile('*.tif', 'Select all image files', 'MultiSelect', 'on');

% Check if multiple files were selected
if iscell(files)
    num_files = length(files);
else
    num_files = 1;
    files = {files}; % Convert to cell for consistency
end

% Prompt user to enter the scratch orientation
orient = input("Enter the orientation of the scratch ('vertical' or 'horizontal'): ", 's');

% Validate orientation input
if ~strcmpi(orient, "vertical") && ~strcmpi(orient, "horizontal")
    warning("Invalid entry for scratch orientation. Defaulting to 'vertical'.");
    orient = "vertical";
end

% Prompt for optional cropping along the manually drawn boundaries
fprintf("The following prompts allow you to apply additional cropping if needed.\n");
fprintf("Cropping is measured in pixels from the edges of the image. Enter 0 for none.\n");

crop_top = input("Enter the number of pixels to crop from the top (for vertical scratch) or left (for horizontal scratch) side: ");
crop_bottom = input("Enter the number of pixels to crop from the bottom (for vertical scratch) or right (for horizontal scratch) side: ");



%% =============================================================================
% 2. Image Preprocessing and Binary Line Extraction

% Initialize an array to store the areas
areas = zeros(num_files, 1);
distances = zeros(num_files, 1);

for i = 1:num_files
    % Read the image
    img = imread(fullfile(pathname, files{i}));
    % Convert the image to grayscale
    gray_img = rgb2gray(img);
    
    % Threshold to extract white boundary lines (value may need adjustment)
    binary_img = gray_img > 230;

    % Apply vertical cropping if specified
    if orient == "vertical"
        binary_img(1:crop_top, :) = 0; % Exclude rows from 1 to crop_top
        binary_img(end - crop_bottom + 1:end, :) = 0; % Exclude rows from crop_bottom to end
    elseif orient == "horizontal"
        binary_img(:, 1:crop_top) = 0; % Exclude columns from 1 to crop_top
        binary_img(:, end - crop_bottom + 1:end) = 0; % Exclude columns from crop_bottom to end
    end


    % Thin white lines to skeletons
    thinned_lines = bwmorph(binary_img, 'thin', Inf);
    %imshow(thinned_lines)

    % Check that the skeleton contains exactly two connected components
    ccLines = bwconncomp(thinned_lines, 8);    % 8‑connected to avoid diagonal breaks
    if ccLines.NumObjects ~= 2
        error("Scratch‑wound QC: expected 2 boundary lines, but detected %d. " + ...
            "Review the threshold or redraw the lines.", ccLines.NumObjects);
    end
    

    % Detect endpoints of the two manual boundaries
    endpoints = bwmorph(thinned_lines, 'endpoints');

    % Check that exactly four endpoints exist (two per line)
    nEndpoints = nnz(endpoints);

    if nEndpoints ~= 4
        error("Scratch‑wound QC: expected 4 line endpoints, but detected %d. " + ...
            "Ensure each boundary line is continuous and has two visible ends.", nEndpoints);
    end

    % Extract 4 endpoints
    [end_y, end_x] = find(endpoints);

    end1 = [end_y(1), end_x(1)];
    end2 = [end_y(2), end_x(2)];
    end3 = [end_y(3), end_x(3)];
    end4 = [end_y(4), end_x(4)];

% =============================================================================
% 3. Endpoint and Line Labeling

    % Find and label connected components to determine which endpoints belong together
    % (Label line 1 and line 2 in the binary matrix)
    cc = bwconncomp(thinned_lines);
    labels = labelmatrix(cc);

    % Extract the labels for the give endpoints
    label_end1 = labels(end1(1), end1(2));
    label_end2 = labels(end2(1), end2(2));
    label_end3 = labels(end3(1), end3(2));
    label_end4 = labels(end4(1), end4(2));

    % Identify which two endpoints belong to the same boundary line
    if label_end1 == label_end2
        % If end1 and end2 belong to the same line (line1)
        line1_end1 = end1; line1_end2 = end2;
    elseif label_end1 == label_end3
        % If end1 and end3 belong to the same line (line1)
        line1_end1 = end1; line1_end2 = end3;
    else
        % If end1 and end4 belong to the same line (line1)
        line1_end1 = end1; line1_end2 = end4;
    end

    % For the second line (line2)
    remaining_ends = [end2; end3; end4];
    remaining_labels = [label_end2, label_end3, label_end4];

    % Find which two of the remaining endpoints belong to line2
    line2_idx = find(remaining_labels ~= label_end1);

    % Assign the endpoints to line2
    line2_end1 = remaining_ends(line2_idx(1), :);
    line2_end2 = remaining_ends(line2_idx(2), :);

% =============================================================================
% 4. Connecting Boundary Endpoints and ROI Filling

    % Connect endpoints to form the two connecting boundaries

    % Compute distances between line endpoints
    dists = [
        pdist2(line1_end1, line2_end1), pdist2(line1_end1, line2_end2);
        pdist2(line1_end2, line2_end1), pdist2(line1_end2, line2_end2)
    ];

    [min_dist1, min_idx1] = min(dists(1,:)); % min distance for line1_end1 to line2
    [min_dist2, min_idx2] = min(dists(2,:)); % min distance for line1_end2 to line2

    % Determine which pair is the closest pair
    if min_idx1 == 1
        closest_pair1 = [line1_end1; line2_end1];
    elseif min_idx1 == 2
        closest_pair1 = [line1_end1; line2_end2];
    end

    if min_idx2 == 1
        closest_pair2 = [line1_end2; line2_end1];
    elseif min_idx2 == 2
        closest_pair2 = [line1_end2; line2_end2];
    end


    % Draw connecting lines between boundary endpoints
    % Initialize the line image and mask
    line_img = thinned_lines;
    line_mask = false(size(line_img));

    % Draw the first line using closest_pair1
    [rr1, cc1] = drawline(closest_pair1(1, 2), closest_pair1(1, 1), closest_pair1(2, 2), closest_pair1(2, 1));
    line_mask(sub2ind(size(line_img), rr1, cc1)) = true;

    % Draw the second line using closest_pair2
    [rr2, cc2] = drawline(closest_pair2(1, 2), closest_pair2(1, 1), closest_pair2(2, 2), closest_pair2(2, 1));
    line_mask(sub2ind(size(line_img), rr2, cc2)) = true;

    % Combine the original thinned lines with the new lines
    line_img = line_img | line_mask;

    % Fill the enclosed region (ROI) formed by the two manual and two connecting boundaries
    filled_img = imfill(line_img, 'holes');

% =============================================================================
% 5. Wound Area Calculation and Overlay Visualization

    % Calculate the area of the ROI
    area_in_pixels = sum(filled_img(:));
    disp(['Area of ROI in ', files{i}, ': ', num2str(area_in_pixels), ' pixels']);
    
    % Store the area in the array
    areas(i) = area_in_pixels;

    % Highlight the ROI on the original image
    red_channel = img(:,:,1);
    green_channel = img(:,:,2);
    blue_channel = img(:,:,3);

    red_channel(filled_img) = 0;   % highlight the ROI with a specified color. Change if needed
    green_channel(filled_img) = 0;
    blue_channel(filled_img) = 128;

    highlighted_img = cat(3, red_channel, green_channel, blue_channel);


%% =============================================================================
% 6. Wound length (Midpoint Distance) Measurement

    % Calculate midpoints of line1 and line2 (two possibilities)
    midpoint_line1_1 = round([(line1_end1(2) + line2_end1(2)) / 2, (line1_end1(1) + line2_end1(1)) / 2]);
    midpoint_line2_1 = round([(line1_end2(2) + line2_end2(2)) / 2, (line1_end2(1) + line2_end2(1)) / 2]);
    
    midpoint_line1_2 = round([(line1_end1(2) + line2_end2(2)) / 2, (line1_end1(1) + line2_end2(1)) / 2]);
    midpoint_line2_2 = round([(line1_end2(2) + line2_end1(2)) / 2, (line1_end2(1) + line2_end1(1)) / 2]);


    % Calculate the distance between the midpoints
    distance_between_midpoints_1 = pdist2(midpoint_line1_1, midpoint_line2_1);
    distance_between_midpoints_2 = pdist2(midpoint_line1_2, midpoint_line2_2);

    if distance_between_midpoints_2 > distance_between_midpoints_1
        distance_between_midpoints = distance_between_midpoints_2;
        midpoint_line1 = midpoint_line1_2;
        midpoint_line2 = midpoint_line2_2;
    else
        distance_between_midpoints = distance_between_midpoints_1;
        midpoint_line1 = midpoint_line1_1;
        midpoint_line2 = midpoint_line2_1;
    end
    
    % Display the distance
    disp(['Distance between the midpoints of top and bottom lines: ', num2str(distance_between_midpoints), ' pixels']);
    
    % Store the distance in an array
    distances(i) = distance_between_midpoints;
    
    % Draw the line connecting the midpoints in red
    [rr_mid, cc_mid] = drawline(midpoint_line1(1), midpoint_line1(2), midpoint_line2(1), midpoint_line2(2));
    highlighted_img(sub2ind(size(highlighted_img), rr_mid, cc_mid, 1*ones(size(rr_mid)))) = 255; % Red channel
    highlighted_img(sub2ind(size(highlighted_img), rr_mid, cc_mid, 2*ones(size(rr_mid)))) = 0;   % Green channel
    highlighted_img(sub2ind(size(highlighted_img), rr_mid, cc_mid, 3*ones(size(rr_mid)))) = 0;   % Blue channel

% =============================================================================
% 7. Export of Annotated Image and Summary Table
    
    figure; 
    subplot(1,2,1); imshow(img); title('Original image')
    subplot(1,2,2); imshow(highlighted_img); title('ROI')

    % Get the original file name without extesnsion
    [~, name, ~] = fileparts(files{i});
    % Export image
    imwrite(highlighted_img, fullfile(pathname, [name '_highlightROI.tif']), 'tiff', 'Resolution', [300 300]);
end

% Derive common prefix for output summary
common_to_use = files{1}(1:find(any(diff(char(files(:))),1),1,'first')-1);
if isempty(common_to_use)
   common_to_use = name;
end

% Export table of area and distance
output_excel_filename = fullfile(pathname, [common_to_use '_result.xlsx']);
T = table(files(:), areas, distances, ...
    'VariableNames', {'Filename', 'ROI_Area_Pixels', 'ROI_Distance_Pixels'});
writetable(T, output_excel_filename);
disp(['Summary saved to ', output_excel_filename]);

%% =============================================================================
% 8. Utility Function for Drawing Lines

function [rr, cc] = drawline(x1, y1, x2, y2)
    % Calculate the pixels that form a straight line between (x1, y1) and (x2, y2)
    dx = abs(x2 - x1); dy = abs(y2 - y1);
    sx = sign(x2 - x1); sy = sign(y2 - y1);
    err = dx - dy;
    rr = y1; cc = x1;

    while x1 ~= x2 || y1 ~= y2
        e2 = 2 * err;
        if e2 > -dy, err = err - dy; x1 = x1 + sx; end
        if e2 < dx, err = err + dx; y1 = y1 + sy; end
        rr = [rr; y1]; cc = [cc; x1];
    end
end









