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

function psf=PSFcalculation(pixelSize, wavelength, NA,focusDepth,n)

if n<=256
    gridSize=ceil(n/4);
else
    gridSize=64;
end

sinThetaSquared = (1 - (1 - NA^2)) / 2;
opticalPathDiff = 8 * pi * focusDepth * sinThetaSquared / wavelength;

gridRange = -gridSize * pixelSize : pixelSize : gridSize * pixelSize;
[gridX, gridY] = meshgrid(gridRange, gridRange);
[~, radialDistance] = cart2pol(gridX, gridY);

apertureMask = radialDistance <= 1;
radialDistance(~apertureMask) = 0; 

validRadii = radialDistance(apertureMask);
integrationPoints = linspace(0, 1, 100); 
[radiusGrid, integrationGrid] = meshgrid(validRadii, integrationPoints); 
kernelFunction = 2 * exp(1i * opticalPathDiff * (integrationGrid.^2) / 2) ...
                 .* besselj(0, 2 * pi * radiusGrid * NA / wavelength .* integrationGrid);

integratedValues = trapz(integrationPoints, kernelFunction, 1);
fieldAmplitude = zeros(size(radialDistance));
fieldAmplitude(apertureMask) = real(integratedValues); 


psf = abs(fieldAmplitude).^2;
psf = psf / sum(psf(:)); 
end
