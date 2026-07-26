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

function [results]=SPIFFI_functions(imgstack, cfg)

if cfg.ifAC2
    cfg.ifPreDeconv=false;
    cfg.pxSize = cfg.pxSize./cfg.M; 
    cfg.ifAC=false;
    
    Mag=cfg.M2;
    ACorder=cfg.order2;
    deconvPost=cfg.deconvPost2;
    PSFshrink=1;
    deconvACorder=cfg.order*cfg.order2;

else
    Mag=cfg.M;
    ACorder=cfg.order;
    deconvPost=cfg.deconvPost;
    PSFshrink=cfg.PSFshrink;
    deconvACorder=cfg.order;

end
if cfg.ifT && size(imgstack,3)>20
    cfg.ifAC=false;
end
%% remove bkg
if cfg.ifbkg
    if isfield(cfg, 'bkgradius')
        bkgradius = cfg.bkgradius;
    else
        bkgradius=8*ceil((cfg.lambda)/(2*cfg.NA)/cfg.pxSize);
    end
    bkgradius=min(bkgradius, min(ceil(size(imgstack,[1,2])/2)));
    [imgstack,~] = rollingBall(imgstack, bkgradius,  true);
    smartDisp('Background subtraction')
end
%% hist match
maxi=max(imgstack,[],"all");
stack=zeros(size(imgstack));
imgstack(:,:,1)=imgstack(:,:,1)/maxi;
for ff=2:size(imgstack,3)
    imgstack(:,:,ff)=imhistmatch(imgstack(:,:,ff)/maxi,imgstack(:,:,1),'method','uniform');
end
imgstack(imgstack<0)=0;
for i = 1 : size(imgstack,3)
    stack(:,:,i) = imgstack(:,:,i) - min(min(imgstack(:,:,i)));
end

%% psf   
psf = PSFcalculation(cfg.pxSize*PSFshrink, cfg.lambda, cfg.NA, 0, min(size(stack,1),size(stack,2)));
psf = psf.^(sqrt(deconvACorder./ACorder));
psf = psf./sum(psf(:));

psfAC = PSFcalculation(cfg.pxSize*PSFshrink/Mag, cfg.lambda, cfg.NA, 0, min(size(stack,1),size(stack,2)));
psfAC = psfAC.^(sqrt(deconvACorder));
psfAC = psfAC./sum(psfAC(:));
%% pre deconvolution
stackD=zeros(size(stack));
if cfg.ifPreDeconv
    smartDisp(['Pre deconvolution'])
    for i = 1:size(stack,3)
        stackD(:,:,i)=deconvlucy(stack(:,:,i),psf,cfg.deconvPre);
    end
else 
    stackD=stack;
end
clear stack
%% interpolation

if Mag == 1
    stackI = stackD;
else
    stackI = fourierInterpolation(stackD, [Mag, Mag, 1], 'lateral');
    stackI(stackI < 0) = 0;
    smartDisp('Fourier interpolation');
end
clear stackD
%% pAC
N=size(stackI,3);
stackI=stackI-cfg.ACscale*mean(stackI,3);
if cfg.ifAC
    switch ACorder
        case 3
            AC3=[]; 
            for k1 = 0:N-1
                for k2 = 0:N-1
                    if k1+k2<N
                        ACtemp=sum(stackI(:,:,1:N-max(k1,k2)).*stackI(:,:,k1+1:N-max(k1,k2)+k1).*stackI(:,:,k2+1:N-max(k1,k2)+k2),3);
                        AC3=cat(3,AC3,ACtemp);
                    end
                end
            end
            smartDisp(['Result Calculation: AC3'])
            pAC=abs(mean(AC3,3));
        case 2
            AC2 = [];
            for k=0:N-1
                AC2=cat(3,AC2,sum(stackI(:,:,1:N-k) .* stackI(:,:,k+1:N),3));
            end
            smartDisp(['Result Calculation: AC2'])
            pAC=abs(mean(AC2,3));
    end
else
    switch ACorder
        case 2
            pAC=abs(mean(stackI(:,:,1:end-1).*stackI(:,:,2:end),3));
            smartDisp(['Result Calculation: Moment-2'])

        case 3
            pAC = abs(mean(stackI(:,:,1:end-2).*stackI(:,:,2:end-1).*stackI(:,:,3:end),3));
            smartDisp(['Result Calculation: Moment-3'])
    end
end
%% AC deconv
smartDisp(['Post deconvolution'])
pACtemp = deconvlucy(pAC,psfAC, deconvPost);

%%
results=pACtemp;
smartDisp('-------------------------')


end
