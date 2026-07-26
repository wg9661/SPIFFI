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


function [Image, Info]=readAll(filePath)
    Info = imfinfo(filePath);
    Zstacks = size(Info,1);
    Width = Info.Width;
    Height = Info.Height;
    Image = zeros(Height, Width, Zstacks);
    for i=1:Zstacks
        Image(:,:,i)=imread(filePath,i);    
    end
end