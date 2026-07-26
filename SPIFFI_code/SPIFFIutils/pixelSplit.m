%   Split a 3D image into four stacks using checkerboard pixel sampling
%  and reconstruct missing pixels by linear interpolation.
%
%   Input:
%       img3d  - A 3D image array of size (rows × cols × frames), where each
%                slice represents one frame of the image sequence.
%   Output:
%       img1   - Interpolated image stack using pixels from even positions
%                (X + Y is even) in the original coordinate system.
%       img2   - Interpolated image stack using pixels from odd positions
%                (X + Y is odd) in the original coordinate system.
%       img3   - Interpolated image stack using even-position pixels with
%                swapped coordinates (Y, X) to simulate rotated sampling.
%       img4   - Interpolated image stack using odd-position pixels with
%                swapped coordinates (Y, X).

% ********************************************************************************
% Copyright (c) [2025] [Wei Guo, et al. @ EPFL]. All Rights Reserved.
% wei.guo@epfl.ch
%
% This repository accompanies the manuscript: 
% "Spatial Polarization-Induced Fluorescence Fluctuation Imaging (SPIFFI) Enables Single-shot Super-Resolution and Multidimensional Imaging"
% bioRxiv 2025. DOI: 10.64898/2025.12.12.693764
% 
% RESTRICTIONS:
% **UNAUTHORIZED** copying, modification, distribution, or use for any other 
% **COMMERCIAL PURPOSE** is strictly **PROHIBITED**.
% ********************************************************************************

function [img1,img2,img3,img4]=pixelSplit(img3d)
[rows, cols, frames] = size(img3d);
z = frames;
[X, Y] = meshgrid(1:cols, 1:rows);
maskA = mod(X + Y, 2) == 0; % even
maskB = ~maskA;             % odd

img1 = zeros(rows, cols, frames);
img2 = zeros(rows, cols, frames);
img3 = zeros(rows, cols, frames);
img4 = zeros(rows, cols, frames);
%frames
for k = 1:frames
    slice = double(img3d(:,:,k));
    
    % img1 & img2:original coordination
    img1(:,:,k) = interpolate_logic(slice, maskA, X, Y);
    img2(:,:,k) = interpolate_logic(slice, maskB, X, Y);
    
    % img3 & img4: swith X and Y to simulate the rotation,
    % avoid the image shape problem
    img3(:,:,k) = interpolate_logic(slice, maskA, Y, X); 
    img4(:,:,k) = interpolate_logic(slice, maskB, Y, X);
end

imgSq = cat(3, img1, img2, img3, img4);

% Invalid value
if any(~isfinite(imgSq(:)))
    minFiniteValue = min(imgSq(isfinite(imgSq)));
    imgSq(~isfinite(imgSq)) = minFiniteValue;
end

% results
img1 = imgSq(:,:,1:z);
img2 = imgSq(:,:,1+z:2*z);
img3 = imgSq(:,:,1+z*2:3*z);
img4 = imgSq(:,:,1+z*3:end);

function img_filled = interpolate_logic(img_slice, mask, gridX, gridY)
    % known value
    known_x = gridX(mask);
    known_y = gridY(mask);
    known_values = img_slice(mask);
    
    % unknown value
    unknown_x = gridX(~mask);
    unknown_y = gridY(~mask);
    
    % griddata
    interpolated_values = griddata(known_x, known_y, known_values, unknown_x, unknown_y, 'linear');
    
    % fill
    img_filled = img_slice;
    img_filled(~mask) = interpolated_values;
end

end