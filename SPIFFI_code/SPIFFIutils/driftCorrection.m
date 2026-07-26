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


%  Performs subpixel correction on images using fiducial drift data
% % Input:
% Stack - The stack of raw images, with dimensions [height, width, nFrames]
% driftX - The x-direction drift data (in pixels) for each frame, with the same length as the number of frames.
% driftY - The y-direction drift data (in pixels) for each frame, with the same length as the number of frames.
% driftInd - The index of frames
% Output:
% correctedStack - The stack of images after drift correction.
% % Note:
% This function assumes driftX and driftY are absolute drift data, using the first frame as a reference.
%For the i-th frame, the corrected displacement is calculated as follows:
%shiftX = -(driftX(i) - driftX(1))
%shiftY = -(driftY(i) - driftY(1))
% Frequency domain shift correction is achieved by calling the internal function fourierShift, thereby achieving sub-pixel level alignment.


function correctedStack = driftCorrection(stack, driftX, driftY,driftInd)
    [h, w, nFrames] = size(stack);
    correctedStack = zeros(h, w, nFrames);
    
    driftX=driftX-driftX(1);
    driftY=driftY-driftY(1);
    driftInd = truncateAtFirstNaN(driftInd);
    driftX = truncateAtFirstNaN(driftX);
    driftY = truncateAtFirstNaN(driftY);

    driftX = double(driftX);
    driftY = double(driftY);
    
    driftX=interp1(driftInd,driftX,1:nFrames);
    driftY=interp1(driftInd,driftY,1:nFrames);

    refDriftX = driftX(1);
    refDriftY = driftY(1);
    
    for i = 1:nFrames
        shiftX = -(driftX(i) - refDriftX);
        shiftY = -(driftY(i) - refDriftY);
        % Subpixel correction is achieved using frequency domain shifting (shift format: [vertical, horizontal])
        correctedStack(:,:,i) = fourierShift(stack(:,:,i), [shiftY, shiftX]);
    end

    fprintf('Drift corrected!\n');
end

function shiftedImg = fourierShift(img, shift)
% performs subpixel translation on an image using the Fourier translation theorem
% % Input:
% img - Input 2D image
% shift - [vertical, horizontal] Translation amount (unit: pixels, can be a non-integer)
% % Output:
% shiftedImg - The translated image
% % Explanation:
% This function first calculates the FFT of the image, then multiplies it by a phase factor to achieve translation,
% Finally, it uses the inverse FFT to recover the shifted image.

    [h, w] = size(img);
    
    % mesh a frequency grid (X corresponds to the horizontal direction, Y corresponds to the vertical direction).
    [X, Y] = meshgrid(-floor(w/2):(w-floor(w/2)-1), -floor(h/2):(h-floor(h/2)-1));
    
    fftImg = fftshift(fft2(img));
    
    % phase factor：exp(-i*2*pi*(X*dx/w + Y*dy/h))
    % shift format [dy, dx]
    dx = shift(2);
    dy = shift(1);
    phaseFactor = exp(-1i * 2 * pi * ((X * dx)/w + (Y * dy)/h));

    shiftedFFT = fftImg .* phaseFactor;
    shiftedImg = real(ifft2(ifftshift(shiftedFFT)));
end

