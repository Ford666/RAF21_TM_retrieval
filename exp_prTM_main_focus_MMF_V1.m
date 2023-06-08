% Phase retrieval for complex TM recovery of a MMF
% Compare prVBEM, GGS2-1, RAF, RAF2-1 in focusing performance
% This includes single-focusing boxplot and multi-focusing image.

% Author: Shengfu Cheng
% Updated Date: 08/03/203

addpath('./utils');
addpath('./DMDCode');
addpath('./algs/');
close all

DMD_close;
clearvars -except SLMmsk CGH_calib

%% DMD allocation or halt 
% load DMD
if ~libisloaded('DMD')
    loadlibrary('C:\Program Files\ALP-4.3\ALP-4.3 API\x64\alp4395.dll','C:\Program Files\ALP-4.3\ALP-4.3 API\alp.h','alias','DMD');  
end

if ~exist('deviceid', 'var') %DMD allocation for the first exp
    deviceid = alpdevicealloc;
else                                          %DMD halt for the next exp               
    alp_returnvalue = alpprojhalt(deviceid);
    if exist('FocSeqid', 'var')
        alp_returnvalue = AlpSeqFree(deviceid, FocSeqid); %free memory of existing sequence
        clearvars FocSeqid;
    end
end


%% Parameters
bitplanes = 1;
picturetime = 4000; % us
repeattimes = 999999; 

N = 32; %Input size
gamma = 5; %sampling ratio
P = gamma * N^2; %measurement number
Lxy = 28;
alpha = 1/4; beta = 1/4;

expType = 2; %1: single-focusing; 2: multi-focusing

fiberAlign = 0;

%% camera ready
imaqreset;
camVid = videoinput('mwspinnakerimaq', 1, 'Mono16');
camSrc = getselectedsource(camVid);
if fiberAlign
    ROI_pos = [0 0 720 540];
    
else
    xoffset = 328; yoffset = 296; M = 140; ROI_pos = [xoffset yoffset M M];
end

camSrc.ExposureTime = 10002;
camVid.ROIPosition = ROI_pos;


%% fiber alignment
if ~exist('AllOnSeqid', 'var')
    CGH_allon = 128 * CudaLee_1920_1080(zeros(N, N, 'single'),Lxy, alpha, beta, 1);
    AllOnSeqid = alpseqalloc(deviceid,bitplanes,1);
    alp_returnvalue = alpseqput(deviceid,AllOnSeqid,CGH_allon);
end

alp_returnvalue = alpseqtiming(deviceid,AllOnSeqid,picturetime,picturetime-100);
alp_returnvalue = alpseqcontrol(deviceid,AllOnSeqid,repeattimes); %repeat sequence displaying
alp_returnvalue = alprojstart(deviceid,AllOnSeqid);
Icam_allon = im2double(getsnapshot(camVid));


if exist('deviceid', 'var') %DMD halt
    alp_returnvalue = alpprojhalt(deviceid);
end

if fiberAlign
    I1 = medfilt2(Icam_allon, [5, 5]); %median filtering
    I2 = imbinarize(I1);
    R = regionprops(bwlabel(I2),'Centroid','Area');
    [~, maxRegion] = max(cat(1,R.Area));
              
    centers = round(R(maxRegion).Centroid);
    radius = sqrt(1/pi*R(maxRegion).Area);
    
    img_size = radius*2.6;
    M = round(img_size); 
    xoffset = round(centers(1) - img_size/2);
    yoffset = round(centers(2) - img_size/2);  
end

if mod(M, 4) ~= 0; M = 4*round(M/4); end
if mod(xoffset, 4) ~= 0; xoffset = 4*round(xoffset/4); end
if mod(yoffset, 4) ~= 0; yoffset = 4*round(yoffset/4); end

[X, Y] = meshgrid(1:M);

if fiberAlign
    figure('color', [1 1 1], 'position', [150 150 700 480]),
    imagesc(uint8(255*Icam_allon)); title('Fiber alignment')
    colorbar('off'); colormap(gray); hold on
    plot([xoffset, xoffset, xoffset+M, xoffset+M, xoffset], ...
        [yoffset, yoffset+M, yoffset+M, yoffset, yoffset], 'r-', 'linewidth', 1.5); hold off

else
    figure('color', [1 1 1], 'position', [150 150 550 480]),
    imagesc(uint8(255*Icam_allon)); title('Fiber alignment')
    colorbar('off'); colormap(gray); 
end
pause(1)

fprintf("Fiber alignment done, with an area of %dx%d \n", M, M);

%% Indices of focusing
if expType==1 %single-focusing
    Nscan = 20;
    iX = linspace(ceil((1/2-sqrt(2)/4)*M), ceil((1/2+sqrt(2)/4)*M), Nscan+2); 
    iX = ceil(iX(2:end-1)); iY = iX;
    [X, Y] = meshgrid(iX, iY); focInd = sub2ind([M, M], Y, X); 
    focNum = numel(focInd);
    fprintf("Experiment %d: single focusing boxplot\n", expType);

elseif expType==2 %multi-focusing
    cordStar = {[81, 48], [76, 60], [73, 72], [60, 72], [48, 72], [58, 80], [67, 87], [64, 96], [60 110], [71, 103], [81, 96], ...
        [90, 103], [102, 110], [97, 97], [94, 87], [103, 80], [113, 72], [100, 72], [89, 72], [85, 60]};
    focInd = [];
    id2 = 1;
    for id = 1:length(cordStar)
        cord = round(cordStar{id}./160*M); 
        focInd(id2) = (cord(1)-1)*M+cord(2); id2=id2+1;%center 
    end

    focNum = 1; 
    fprintf("Experiment %d: multi-focusing image\n", expType);

end


%% Camera settings
imaqreset;
camVid = videoinput('mwspinnakerimaq', 1, 'Mono16');
camSrc = getselectedsource(camVid);
ROI_pos = [xoffset yoffset M M];
camVid.ROIPosition = ROI_pos;

trigNum = P + 4*focNum;
triggerconfig(camVid, 'hardware');
camSrc.TriggerMode = 'On';
camVid.FramesPerTrigger = 1;
set(camVid, 'TriggerRepeat', trigNum-1);
set(camVid, 'Timeout', 200); 

camSrc.TriggerSource = 'Line2';
camSrc.AdcBitDepth = 'Bit12';
camSrc.ExposureAuto = 'Off';
camSrc.ExposureTime = 53;
camSrc.GainAuto = 'Off';
camSrc.Gain = 0;
camSrc.GammaEnable = 'false';


start(camVid);
fprintf('Camera for TM is ready\n');

%% Load CGH onto DMD and measure TM
if ~(exist('CGH_calib', 'var') && size(CGH_calib, 3) == P)
    SLMmsk = exp(1i*2*pi*rand(N, N, P, 'single'));
    CGH_calib = 128 * CudaLee_1920_1080(angle(SLMmsk), Lxy, alpha, beta, 0);
    SLMmsk = reshape(SLMmsk, [], P);
end

% upload TM calibration patterns
[~,~,picnum] = size(CGH_calib);
if ~(exist('TMSeqid', 'var') && exist('CGH_calib', 'var'))  %upload H basis only at 1st exp.
    TMSeqid = alpseqalloc(deviceid,bitplanes,picnum);
    alp_returnvalue = alpseqput(deviceid,TMSeqid,CGH_calib);
end

fprintf ('Calibration CGH into DMD done\n');

alp_returnvalue = alpseqtiming(deviceid,TMSeqid,picturetime,picturetime-100);
alp_returnvalue = alprojstart(deviceid,TMSeqid);
datav = getdata(camVid, picnum);

fprintf ('Camera for TM calibration done\n');

Icam_speck = single(datav(:, :, 1:picnum)); clear datav;
Icam_speck = Icam_speck ./ max(Icam_speck,  [], 'all');
SpekAmp = reshape(sqrt(Icam_speck), [], picnum);


%% TM retrieval via different algorithms
disp('Phase retrieval for TM')
iter_VBEM = 100;
nWorker = 16;
optsRAF21 = struct;
optsRAF21.nWorker = nWorker;
optsRAF21.mu = 3;
optsRAF21.useGPU = 1;
optsRAF21.iters = 350;

optsGGS21 = struct;
optsGGS21.nWorker = nWorker;
optsGGS21.useGPU = 1;
optsGGS21.iters = 350;


subSpekAmp = SpekAmp(focInd, :);

%prVBEM
[TM_vbem_sub, ~, ~] = prVBEM_tm(subSpekAmp, SLMmsk, iter_VBEM, 0);
      
%GGS2-1
[TM_ggs21_sub, ~, ~] = GGS21_tm(subSpekAmp, SLMmsk, optsGGS21);

% RAF
[TM_raf_sub, ~, ~] = RAF_tm(subSpekAmp.^2, SLMmsk, optsRAF21);

% RAF2-1
[TM_raf21_sub, ~, ~] = RAF21_tm(subSpekAmp.^2, SLMmsk, optsRAF21);


%% Comparison of focusing performance
if expType==1
    focPha_vbem = reshape(TM_vbem_sub', N, N, []); 
    focPha_ggs21 = reshape(TM_ggs21_sub', N, N, []); 
    focPha_raf = reshape(TM_raf_sub', N, N, []); 
    focPha_raf21 = reshape(TM_raf21_sub', N, N, []); 

elseif expType==2
    focPha_vbem = reshape(sum(TM_vbem_sub', 2), N, N);
    focPha_ggs21 = reshape(sum(TM_ggs21_sub', 2), N, N);
    focPha_raf = reshape(sum(TM_raf_sub', 2), N, N);
    focPha_raf21 = reshape(sum(TM_raf21_sub', 2), N, N);
end

CGH_focus = 128 * CudaLee_1920_1080(cat(3, angle(focPha_vbem), angle(focPha_ggs21), ...
    angle(focPha_raf), angle(focPha_raf21)), Lxy, alpha, beta, 0); 

% Focusing
camSrc.ExposureTime = 11;
 
FocSeqid = alpseqalloc(deviceid,bitplanes,4*focNum);
alp_returnvalue = alpseqput(deviceid,FocSeqid,CGH_focus);
fprintf ('Focus CGHs into DMD done\n');

alp_returnvalue = alpseqtiming(deviceid,FocSeqid,picturetime,picturetime-100);
alp_returnvalue =  alpseqcontrol(deviceid,FocSeqid,20);
alp_returnvalue = alprojstart(deviceid,FocSeqid);

focImg = getdata(camVid,4*focNum); 
Ifoc_vbem = single(focImg(:, :,1:focNum)); 
Ifoc_ggs21 = single(focImg(:, :,focNum+1:2*focNum));  
Ifoc_raf = single(focImg(:, :,2*focNum+1:3*focNum));
Ifoc_raf21 = single(focImg(:, :,3*focNum+1:4*focNum)); 

delete(camVid)
clear camVid

fprintf ('Focusing done\n');


%% Show results
pbr_ref = pi/4*(N^2-1)+1;

if expType==1 %single-focusing

    % PBR list
    Ifoc_vbem_max = max(Ifoc_vbem, [], 3); 
    Ifoc_ggs21_max = max(Ifoc_ggs21, [], 3);
    Ifoc_raf_max = max(Ifoc_raf, [], 3); 
    Ifoc_raf21_max = max(Ifoc_raf21, [], 3);
    pbrs_vbem = squeeze(max(Ifoc_vbem,[], [1, 2])./mean(Ifoc_vbem,[1, 2])); 
    pbrs_ggs21 = squeeze(max(Ifoc_ggs21,[], [1, 2])./mean(Ifoc_ggs21,[1, 2])); 
    pbrs_raf = squeeze(max(Ifoc_raf,[], [1, 2])./mean(Ifoc_raf,[1, 2])); 
    pbrs_raf21 = squeeze(max(Ifoc_raf21,[], [1, 2])./mean(Ifoc_raf21,[1, 2]));
        
    clearvars Ifoc_vbem Ifoc_ggs21 Ifoc_raf Ifoc_raf21


    fig1 = figure('color', 'w', 'position', [150 150 600 500]); %box plot of focus intensity
    myBoxPlot([pbrs_vbem./pbr_ref, pbrs_ggs21./pbr_ref, pbrs_raf./pbr_ref, pbrs_raf21./pbr_ref], {'prVBEM', 'GGS2-1', 'RAF', 'RAF21'});
    set(gca,'FontSize',14); ylabel('Algorithms', 'fontsize', 18); ylabel('Normalized focus PBR', 'fontsize', 18);
    %title('Focus PBR distribution by different algorithms', 'fontsize', 20);
    
    fig2 = figure('color', 'w', 'position', [150 150 1200 250]); %focus scanning results
    clim1 = min([min(Ifoc_vbem_max, [], [1, 2]), min(Ifoc_ggs21_max, [], [1, 2]), min(Ifoc_raf_max, [], [1, 2]), min(Ifoc_raf21_max, [], [1, 2])]);
    clim2 = max([max(Ifoc_vbem_max, [], [1, 2]), max(Ifoc_ggs21_max, [], [1, 2]), max(Ifoc_raf_max, [], [1, 2]), max(Ifoc_raf21_max, [], [1, 2])]);
    subplot(141), imagesc(Ifoc_vbem_max); title(sprintf('prVBEM PBR=%.1f', mean(pbrs_vbem)), 'fontsize', 14), caxis([clim1, clim2]);
    subplot(142), imagesc(Ifoc_ggs21_max); title(sprintf('GGS2-1 PBR=%.1f', mean(pbrs_ggs21)), 'fontsize', 14), caxis([clim1, clim2]);
    subplot(143), imagesc(Ifoc_raf_max); title(sprintf('RAF PBR=%.1f', mean(pbrs_raf)), 'fontsize', 14), caxis([clim1, clim2]);
    subplot(144), imagesc(Ifoc_raf21_max); title(sprintf('RAF2-1 PBR=%.1f', mean(pbrs_raf21)), 'fontsize', 14), caxis([clim1, clim2]);
    colormap(gray); 
    h=colorbar('eastoutside','fontsize',15);
    set(h,'Position', [0.92 0.10 0.02 0.75]);


    % save data
    d = datetime; d.Format = 'dd-MM-yyyy, HH-mm'; 
    foldName = sprintf('./Data/Single-focusing/%s', d);
    mkdir(foldName);
    save(sprintf('%s/exp_data.mat', foldName), 'Ifoc_vbem_max', 'pbrs_vbem', ...
                                                    'Ifoc_ggs21_max', 'pbrs_ggs21', ...
                                                    'Ifoc_raf_max', 'pbrs_raf', ...
                                                    'Ifoc_raf21_max', 'pbrs_raf21', 'N', 'gamma');
    saveas(fig1, sprintf('./Data/Single-focusing/%s/boxplot.fig', d))
    saveas(fig2, sprintf('./Data/Single-focusing/%s/focSum.fig', d))

elseif expType==2 %multi-focusing
    % PBR sum & focus uniformity
    pbr_vbem = sum(Ifoc_vbem(focInd) ./ mean2(Ifoc_vbem)); unif_vbem = calcUnif(Ifoc_vbem, focInd);
    pbr_ggs21 = sum(Ifoc_ggs21(focInd) ./ mean2(Ifoc_ggs21)); unif_ggs21 = calcUnif(Ifoc_ggs21, focInd); 
    pbr_raf = sum(Ifoc_raf(focInd) ./ mean2(Ifoc_raf)); unif_raf = calcUnif(Ifoc_raf, focInd); 
    pbr_raf21 = sum(Ifoc_raf21(focInd) ./ mean2(Ifoc_raf21)); unif_raf21 = calcUnif(Ifoc_raf21, focInd);
    

    fig3 = figure('color', 'w', 'position', [150 150 1500 320]); %focus results
	maxFoc = max([max(Ifoc_vbem, [], [1, 2]), max(Ifoc_ggs21, [], [1, 2]), max(Ifoc_raf, [], [1, 2]), max(Ifoc_raf21, [], [1, 2])]);
    subplot(141), imshow(Ifoc_vbem./ maxFoc, []); title('prVBEM'); 
    subplot(142), imshow(Ifoc_ggs21./ maxFoc, []); title('GGS2-1');
    subplot(143), imshow(Ifoc_raf./ maxFoc, []); title('RAF');
    subplot(144), imshow(Ifoc_raf21./ maxFoc, []); title('RAF21');
    colormap(hot);
    h=colorbar('eastoutside','fontsize',15);
    set(h,'Position', [0.92 0.10 0.02 0.85]);
    
    fprintf(sprintf('Multi-focusing pattern, sampling rate = %d, N = %d ', gamma, N^2)); fprintf('\n');
    fprintf(sprintf('prVBEM| PBR efficiency=%.2f, uniformity=%.2f', pbr_vbem/pbr_ref, 100*unif_vbem)); fprintf('\n');
    fprintf(sprintf('GGS21| PBR efficiency=%.2f, uniformity=%.2f', pbr_ggs21/pbr_ref, 100*unif_ggs21)); fprintf('\n');
    fprintf(sprintf('RAF| PBR efficiency=%.2f, uniformity=%.2f', pbr_raf/pbr_ref, 100*unif_raf)); fprintf('\n');
    fprintf(sprintf('RAF21| PBR efficiency=%.2f, uniformity=%.2f', pbr_raf21/pbr_ref, 100*unif_raf21)); fprintf('\n');

    % save data
    d = datetime; d.Format = 'dd-MM-yyyy, HH-mm'; 
    foldName = sprintf('./Data/Multi-focusing/%s', d);
    mkdir(foldName);
    save(sprintf('%s/exp_data.mat', foldName), 'Ifoc_vbem', 'pbr_vbem', 'unif_vbem', ...
                                                'Ifoc_ggs21', 'pbr_ggs21', 'unif_ggs21', ...
                                                'Ifoc_raf', 'pbr_raf', 'unif_raf', ...
                                                'Ifoc_raf21', 'pbr_raf21', 'unif_raf21', 'pbr_ref', 'N', 'gamma');
    saveas(fig3, sprintf('./Data/Multi-focusing/%s/MultiFoc.fig', d))

end