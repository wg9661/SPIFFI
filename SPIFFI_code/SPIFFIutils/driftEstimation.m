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


% driftEstimation: Input a 3D stack, specify the number of segments,
% and directly output the frame-by-frame subpixel drift trajectory.
%
% Input:
% stack - Image stack of [H x W x Frames]
% numBins - Total number of segments to divide into (e.g., 4)
% upsampling - Subpixel precision multiplier (default 10, i.e., 0.1 pixels)
% Output:
% driftX - The x-direction drift data (in pixels) for each frame, with the same length as the number of frames.
% driftY - The y-direction drift data (in pixels) for each frame, with the same length as the number of frames.
% driftInd - The index of frames

function [driftX, driftY,driftInd] = driftEstimation(stack, numBins, upsampling)

    if nargin < 3 || isempty(upsampling), upsampling = 10; end

    [h, w, totalFrames] = size(stack);
    
    binSize = floor(totalFrames / numBins);
    
    if binSize < 1
        error(' numBins is larger than frames!');
    end

    binDrift = zeros(numBins, 2);
    
    % Calculate the reference bin - using the overlay of the first bin as a reference.
    refBinImg = sum(stack(:,:, 1:binSize), 3);
    F_ref = fft2(double(refBinImg));

    % Calculate the displacement of each Bin relative to the first Bin in a loop.
    for b = 1:numBins
        % Get the frame range of the current segment
        idx = (b-1)*binSize + (1:binSize);
        currentBinImg = sum(stack(:,:, idx), 3);
        F_tar = fft2(double(currentBinImg));

        % Pixel-level cross-correlation (rough localization) 
        R = (F_ref .* conj(F_tar));
        CC = ifft2(R);
        [~, max_idx] = max(abs(CC(:)));
        [y_peak, x_peak] = ind2sub([h, w], max_idx);
        
        % Periodic displacement correction
        dx_p = x_peak - 1; if dx_p > w/2, dx_p = dx_p - w; end
        dy_p = y_peak - 1; if dy_p > h/2, dy_p = dy_p - h; end

        % Matrix DFT Subpixel Refinement (Fine Positioning)
        [dx_sub, dy_sub] = dftRefinementCore(R, dy_p, dx_p, upsampling, h, w);
        
        binDrift(b, :) = [dx_p + dx_sub, dy_p + dy_sub];
    end

    % Calculate the frame position at the center of each Bin.
    binMidPoints = ((1:numBins)' - 0.5) * binSize + 0.5;
    queryPoints = (1:totalFrames)';
    
    driftX = interp1(binMidPoints, binDrift(:,2), queryPoints, 'pchip', 'extrap');
    driftY = interp1(binMidPoints, binDrift(:,1), queryPoints, 'pchip', 'extrap');
    driftInd=queryPoints;
    fprintf('Drift estimated!\n');
end


function [dx_sub, dy_sub] = dftRefinementCore(R, row_shift, col_shift, us, h, w)
    grid_size = 1.5; % Search within a 1.5 pixel range around the pixel peak.
    sample_points = ceil(grid_size * us);
    d_grid = (0:sample_points)' / us - grid_size/2;
    
    u = [0:ceil(w/2)-1, -floor(w/2):-1]' / w;
    v = [0:ceil(h/2)-1, -floor(h/2):-1]' / h;
    
    mat_col = exp(2j*pi * u * (col_shift + d_grid'));
    mat_row = exp(2j*pi * v * (row_shift + d_grid'));
    
    CC_sub = abs(mat_row' * R * mat_col);
    [~, loc] = max(CC_sub(:));
    [r_idx, c_idx] = ind2sub(size(CC_sub), loc);
    
    dx_sub = d_grid(c_idx);
    dy_sub = d_grid(r_idx);
end