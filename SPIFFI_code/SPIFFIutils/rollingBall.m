% rollingBall  Background subtraction using the rolling ball algorithm.
%   Input:
%       I            - Input 2D/3D image (rows × cols × slices).
%       radius       - Radius of the rolling ball.
%       isParaboloid - Logical flag; true uses a paraboloid approximation
%                      (default), false uses a spherical ball.
%   Output:
%       output       - Background-corrected image.
%       background   - Estimated background image.

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

function [output, background] = rollingBall(I, radius, isParaboloid)

    % Parabolic surface is used by default
    if nargin < 3, isParaboloid = true; end
    
    [rows, cols, numSlices] = size(I);
    output = zeros(size(I), 'like', I);
    background = zeros(size(I), 'like', I);

    % Shrink Factor
    shrinkFactor = 1;
    if radius > 100
        shrinkFactor = 4; % Further scaling when the radius is extremely large
    elseif radius > 30
        shrinkFactor = 2;
    end
    
    % Construct a structuring element (only needs to be constructed once).
    r_small = radius / shrinkFactor;
    [X, Y] = meshgrid(-fix(r_small):fix(r_small), -fix(r_small):fix(r_small));
    distSq = X.^2 + Y.^2;
    mask = distSq <= r_small^2;
    
    if isParaboloid
        heights = -distSq / (2 * r_small);
    else
        heights = sqrt(max(0, r_small^2 - distSq)) - r_small;
    end
    heights(~mask) = -Inf;
    se = offsetstrel('arbitrary', heights);

    % stack
    for i = 1:numSlices
        I_slice = double(I(:,:,i)); 
        %Scaling
        if shrinkFactor > 1
            I_small = imresize(I_slice, 1/shrinkFactor, 'bilinear');
            bg_small = imopen(I_small, se);
            bg_slice = imresize(bg_small, [rows, cols], 'bilinear');
        else
            bg_slice = imopen(I_slice, se);
        end
        
        % Subtract background and save
        res_slice = I_slice - bg_slice;
        res_slice(res_slice < 0) = 0;
        
        background(:,:,i) = cast(bg_slice, 'like', I);
        output(:,:,i) = cast(res_slice, 'like', I);
    end
end