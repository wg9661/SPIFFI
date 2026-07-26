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

%%
addpath(genpath('./SPIFFIutils'));clear
%% load
fpath='./SPIFFIdata/singleshot';
fname = 'demo_microtubules.tif';
[Img,Info] = readAll([fpath,'/',fname]);
ImgWidth = Info.Width; ImgHeight = Info.Height;

%% images
Img_0=Img(:,:,1);
Img_90=Img(:,:,2);
Img_45=Img(:,:,3);
Img_135=Img(:,:,4);
Img_WF=Img_0+Img_90+Img_45+Img_135; 
%% parameters
cfg.pxSize = 55 *1e-9;  %pixel size
cfg.lambda = 600 *1e-9; % emission wavelength
cfg.NA = 1.49; %NA
cfg.M = 2;  %upsampling, Recommendation: >= 2
cfg.ifPreDeconv = true; %RL deconv
cfg.deconvPre = 9;  %pre iterations
cfg.deconvPost = 3; %post iterations

cfg.ifbkg = false; % remove background by rollingball

cfg.ifT=false; % temporal continuity
cfg.tframes=1; % frames for the reconstruction

cfg.order = 3; % 2 or 3, usually choose 3
cfg.ifAC = true; %true: autocorrelation (slow); false: moment (faster)
cfg.ifAC2 = false; % temporal fluctuation from multiframes
cfg.ACscale = 0.7; % AC1=0-0.8 

cfg.PSFshrink = 1.3;% range 1.2-1.5, usually keep 1.3 for a well calibrated microscope
%%%%%%%%%%%%%%%%%

pstack=[];
idxFrames=[1];
for i = idxFrames
    disp(['Figure index: ', num2str(i)]);
    isWithinTimeRange = (i - cfg.tframes) > 0 && (i + cfg.tframes) <= size(Img_45,3);
    if cfg.ifT && isWithinTimeRange
        disp(['Temporal continuity: ±',num2str(cfg.tframes),' frames'])
        frameRange = i - cfg.tframes : i + cfg.tframes;
    else
        frameRange = i;
    end

    % subimage squence
    img0 = Img_0(:, :, frameRange);
    img45 = Img_45(:, :, frameRange);
    img90 = Img_90(:, :, frameRange);
    img135 = Img_135(:, :, frameRange);
   
    % Pixel Split
    [img01, img02, img03, img04] = pixelSplit(img0);
    [img901, img902, img903, img904] = pixelSplit(img90);
    [img451, img452, img453, img454] = pixelSplit(img45);
    [img1351, img1352, img1353, img1354] = pixelSplit(img135);

    imgSq = cat(3, img0, img45, img90, img135, ...
                   img01, img901, img451, img1351, ...
                   img02, img902, img452, img1352, ...
                   img03, img903, img453, img1353, ...
                   img04, img904, img454, img1354);
    rng(61, 'philox'); %Random Seed
    imgSq = imgSq(:, :, randperm(size(imgSq, 3)));

    [pstack(:,:,i)] = SPIFFI_functions(imgSq(:,:,:),cfg);
end


%% linearization 
% The linearization coefficient could be set to `cfg.order-1` or `cfg.order` 
% It is typically determined by whether the background is subtracted; 
% However, we recommend visually comparing the two in actual use.
SPIFFI=zeros(size(pstack));
bkgRatio=0.01;
for i =idxFrames
    if cfg.ifbkg
        ptemp = pstack(:,:,i).^(1/(cfg.order));
    else
        ptemp = pstack(:,:,i).^(1/(cfg.order-1));
    end
    ptemp(ptemp < cfg.order * bkgRatio * max(ptemp(:))) = 0;
    SPIFFI(:,:,i)=ptemp;
end

%% display
figure;
subplot(1,2,1)
imagesc(Img_WF);axis equal;axis tight;colormap(fire)
title('WF')
subplot(1,2,2);
imagesc(SPIFFI);axis equal;axis tight;colormap(fire)
title('SPIFFI')

