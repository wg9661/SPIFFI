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
addpath(genpath('./SPIFFIutils')); clear
%% load
fpath='./SPIFFIdata/twostage';
[Img_0,~] = readAll([fpath,'/','demo_80nmNanoruler_0.tif']);
[Img_45,~] = readAll([fpath,'/','demo_80nmNanoruler_45.tif']);
[Img_90,~] = readAll([fpath,'/','demo_80nmNanoruler_90.tif']);
[Img_135,Info] = readAll([fpath,'/','demo_80nmNanoruler_135.tif']);
ImgWidth = Info.Width; ImgHeight = Info.Height;
Img_WF=Img_0+Img_90+Img_45+Img_135;

%% parameters
cfg.pxSize = 55 *1e-9;  %pixel size
cfg.lambda = 600 *1e-9; % emission wavelength
cfg.NA = 1.49; %NA
cfg.M = 3;  %upsampling, Recommendation: >= 2
cfg.ifPreDeconv = true; %RL deconv
cfg.deconvPre = 10;  %pre iterations
cfg.deconvPost = 3; %post iterations

cfg.ifbkg = false; % remove background by rollingball

cfg.ifT=false; % temporal continuity
cfg.tframes=1; % frames for the reconstruction

cfg.order = 3; % 2 or 3, usually choose 3
cfg.ifAC = true; %true: autocorrelation (slow); false: moment (faster)
cfg.ifAC2 = false; % temporal fluctuation from multiframes
cfg.ACscale = 0.7; % AC1=0-0.8 

cfg.PSFshrink = 1.3;% range 1.2-1.4, usually keep 1.3 for a well calibrated microscope
%%%%%%%%%%%%%%%%%

pstack=[];
idxFrames=[1:1000];  % frames to be calculated
parfor i = idxFrames
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

%% drift correction
[driftX, driftY,driftInd] = driftEstimation(pstack, 4, 3);
SPIFFI = driftCorrection(pstack, driftX, driftY,driftInd);
figure
plot(driftInd,driftX)
hold on
plot(driftInd,driftY)
legend('DriftX','DriftY')
xlabel('frames')
ylabel('pixels')
title('drift cruves')
%% config for sofi
cfg.ifAC2=true; %Execution temporal flunctuation analysis
cfg.M2 = 1; %magification
cfg.order2 = 3; %order
cfg.deconvPost2 = 1;
cfg.ACscale = 1; % AC2==1==SOFI
disp(['SPIFFI + SOFI with dirft correction']);
[pAC] = SPIFFI_functions(SPIFFI,cfg);

%% SPIFFI+SOFI without drift correction
cfg.ifAC2=true;
cfg.M2 = 1;
cfg.order2 = 3;
cfg.deconvPost2 = 1;
cfg.ACscale = 1; % AC2==1==SOFI
disp(['SPIFFI + SOFI without dirft correction']);
[pACd] = SPIFFI_functions(pstack,cfg);

%% linearization
bkgRatio = 0.01;
SRimg = pAC.^(1/cfg.order2);
SRimg(SRimg < cfg.order2 * bkgRatio * max(SRimg(:))) = 0;

SRimgd = pACd.^(1/cfg.order2);
SRimgd(SRimgd < cfg.order2 * bkgRatio * max(SRimgd(:))) = 0;

%% visulization
subplot(1,3,1);
imagesc(sum(SPIFFI(:,:,idxFrames),3));axis equal;axis tight;colormap(fire)
title('SPIFFIs')

roi1 = [95, 47, 35, 35]; 
roi2 = [15, 130, 35, 35]; 
subplot(1,3,2);imagesc(SRimg);
colormap(fire);axis equal;axis tight
title('SPIFFI + SOFI (Drift Corrected)')
hold on; 
rectangle('Position', roi1, 'EdgeColor', 'w', 'LineWidth', 2, 'LineStyle', '-.');
rectangle('Position', roi2, 'EdgeColor', 'r', 'LineWidth', 2, 'LineStyle', '-');

subplot(1,3,3);imagesc(SRimgd);
colormap(fire);axis equal;axis tight
title('Without Drift Correction')
hold on; 
rectangle('Position', roi1, 'EdgeColor', 'w', 'LineWidth', 2, 'LineStyle', '-.');
rectangle('Position', roi2, 'EdgeColor', 'r', 'LineWidth', 2, 'LineStyle', '-');