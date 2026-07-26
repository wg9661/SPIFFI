% calculate angle and dolp from stokes vector
%   Input:
%       I - raw image
%       cal - intensity calibration
%       SR - superresolution image
%       cfg - config from SPIFFI
%       thresI - intensity thres of SR to show the ROI
%       thresDOP - DOP thres to show the ROI
%   Output:
%       angle_colored - rendered orientation image
%       theta - calcualted orientation 
%       DOP_colored - rendered Dolp image
%       DOP - calcualted DoLP 

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

function [angle_colored,theta,DOP_colored,DOP]=plotAngle(I0, I45, I90, I135, SR, cal, cfg, thresI, thresDOP)
    % remove bkg
    imgstack=cat(3,I0, I45, I90, I135);
    if cfg.ifbkg
        if isfield(cfg, 'bkgradius')
            bkgradius = cfg.bkgradius;
        else
            bkgradius=10*ceil((cfg.lambda)/(2*cfg.NA)/cfg.pxSize);
        end
        bkgradius=min(bkgradius, min(ceil(size(imgstack,[1,2])/2)));
        [imgstack,~] = rollingBall(imgstack, bkgradius,  true);
    end
    I0 = imgstack(:,:,1); I45 = imgstack(:,:,2);I90 = imgstack(:,:,3);I135 = imgstack(:,:,4);

    maxi=max([I0,I90,I45,I135],[],"all");
    I0=I0/maxi;
    I90=imhistmatch(I90/maxi,I0);
    I45=imhistmatch(I45/maxi,I0);
    I135=imhistmatch(I135/maxi,I0);

    %smooth the noise
    psf = PSFcalculation(cfg.pxSize*cfg.PSFshrink, cfg.lambda, cfg.NA, 0, min(size(I0,1),size(I0,2)));
    maskpsf=ceil(size(psf,2)/2)-8:ceil(size(psf,2)/2)+8;
    psf=psf(maskpsf,maskpsf); psf=psf./sum(psf(:));
    I0=deconvlucy(I0, psf,2);
    I90=deconvlucy(I90, psf,2);
    I45=deconvlucy(I45, psf,2);
    I135=deconvlucy(I135, psf,2);
    I0=imfilter(I0, psf.^(1/sqrt(cfg.order)));
    I90=imfilter(I90,psf.^(1/sqrt(cfg.order)));
    I45=imfilter(I45, psf.^(1/sqrt(cfg.order)));
    I135=imfilter(I135, psf.^(1/sqrt(cfg.order)));

    %magnification
    rows = round(size(I0,1) * cfg.M);
    cols = round(size(I0,2) * cfg.M);
    I0 = imresize(I0, [rows, cols], 'bilinear'); 
    I90 = imresize(I90, [rows, cols], 'bilinear'); 
    I45 = imresize(I45, [rows, cols], 'bilinear'); 
    I135 = imresize(I135, [rows, cols], 'bilinear'); 
    
    %cali
    I090=cal.M090*([I0(:),I90(:)]./sum([I0(:),I90(:)],2))';
    I45135=cal.M45135*([I45(:),I135(:)]./sum([I45(:),I135(:)],2))';
    I0=reshape(I090(1,:),size(I0,1),size(I0,2));
    I90=reshape(I090(2,:),size(I90,1),size(I90,2));
    I45=reshape(I45135(1,:),size(I45,1),size(I45,2));
    I135=reshape(I45135(2,:),size(I135,1),size(I135,2));

    %stokes
    S1=(I0-I90)./(I0+I90);
    S2=(I45-I135)./(I45+I135);

    % angle
    theta = rad2deg(0.5 * atan2(S2, S1));
    theta(theta<0)=theta(theta<0)+180;
    % DOLP
    DOP=sqrt(S1.^2 + S2.^2);
    DOP=max(DOP,0);
    DOP=min(DOP,1);
    theta_resized=theta;
    DOP_resized=DOP;
    
    % Brightness threshold mask: Angular color is only displayed when the brightness of the super-resolution image exceeds a threshold.
    mask_DOP = SR > thresI*max(SR(:)); 
    mask_DOP = mask_DOP.*SR;
    target_v = 1.5;
    v =target_v*mask_DOP./max(mask_DOP(:)); % 
    % DOP）
    DOP_colormap = isolum; 
    DOP_image=uint8(DOP_resized*255);
    DOP_colored = ind2rgb(DOP_image, DOP_colormap); %  DOP2RGB 
    DOP_colored=DOP_colored.* v;
    mask_angle=(DOP_resized>thresDOP)&(mask_DOP>0);
    mask_angle=mask_angle.*SR;

    % angle2rgb
    sfactor=1.15;
    vfactor=1.15;
    h = mod(theta_resized/180, 1); 
    s = sfactor * ones(size(h)); 
    v = vfactor*mask_angle./max(mask_angle(:)); % 
    % hsv
    hsv_image = cat(3, h, s, v); 
    angle_colored = hsv2rgb(hsv_image);
    
