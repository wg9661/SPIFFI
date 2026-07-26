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
fpath='./SPIFFIdata/orientation';
fname = 'demo_GUV.tif';
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
cfg.NA = 1.49;
cfg.M = 1;  %upsampling, Recommendation: >= 2
cfg.ifPreDeconv = true; %RL deconv
cfg.deconvPre = 9;  %pre iterations
cfg.deconvPost = 3; %post iterations

cfg.ifbkg = true; %remove background by rollingball
cfg.bkgradius = 10; %radius of rollingball

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
SPIFFI=zeros(size(pstack));
bkgRatio=0.03;
for i =idxFrames
    ptemp = pstack(:,:,i).^(1/(cfg.order));
    ptemp(ptemp < cfg.order * bkgRatio * max(ptemp(:))) = 0;
    SPIFFI(:,:,i)=ptemp;
end
%% intensity calibriation for angular calculation
%cali 0, 90
Ipi=[130, 30; 26, 150]'; %measured intensity
Ipo=[1.0, 0.0; 0.0, 1.0]'; %truth intensity

Ipi=Ipi./sum(Ipi);
cal.M090=Ipo*pinv(Ipi);

% cali 45 135
Ipi2=[220,46; 60, 280]'; %measured intensity
Ipo2=[1, -0; -0, 1]'; %truth intensity
Ipi2=Ipi2./sum(Ipi2);
cal.M45135=Ipo2*pinv(Ipi2);

[angle_colored,angle_raw,DOLP_colored,DOLP_raw]=plotAngle(Img_0(:,:,i), Img_45(:,:,i), Img_90(:,:,i), Img_135(:,:,i), SPIFFI, cal, cfg, 0, 0);

%% display
figure;
subplot(2,2,1)
imagesc(Img_WF);axis equal;axis tight;colormap(fire)
title('WF')
subplot(2,2,2);
imagesc(SPIFFI);axis equal;axis tight;colormap(fire)
title('SPIFFI')

subplot(2,2,3);
imshow(angle_colored);axis equal;axis tight;
title('Angle')
subplot(2,2,4);
imshow(DOLP_colored);axis equal;axis tight;
title('DoLP')
