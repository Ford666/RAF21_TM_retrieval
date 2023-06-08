% Phase retrieval for complex TM recovery of a MMF
% Focusing PBR map using the TM retrieved by RAF2-1, at a certain sampling ratio 

% Author: Shengfu Cheng
% Updated Date: 08/03/203


addpath('./utils');
addpath('./DMDCode');
addpath('./algs/');
close all


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
picturetime = 3500; % us
repeattimes = 999999; 

N = 32; %Input size
gamma = 8; %sampling ratio
P = gamma * N^2; %measurement number
Lxy = 28;
alpha = 1/4; beta = 1/4;

ifScan = 1;
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

    if ifScan
        focInd = find((X-M/2).^2+(Y-M/2).^2 < radius^2);
    else
        focInd = M/2 * M + M/2; 
    end

else
    figure('color', [1 1 1], 'position', [150 150 550 480]),
    imagesc(uint8(255*Icam_allon)); title('Fiber alignment')
    colorbar('off'); colormap(gray); 

    if ifScan
        focInd = find((X-M/2).^2+(Y-M/2).^2 < (3/8*M)^2);
    else
        focInd = M/2 * M + M/2; 
    end

end
pause(1)
fprintf("Fiber alignment done\n");
focNum = numel(focInd);

fprintf("%d points to be scanned inside %dx%d square \n", focNum, M, M);

%% Camera settings
imaqreset;
camVid = videoinput('mwspinnakerimaq', 1, 'Mono16');
camSrc = getselectedsource(camVid);
ROI_pos = [xoffset yoffset M M];
camVid.ROIPosition = ROI_pos;

trigNum = P + focNum;
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

fprintf ('Cmaera for TM calibration done\n');

Icam_speck = single(datav(:, :, 1:picnum)); clear datav;
Icam_speck = Icam_speck ./ max(Icam_speck,  [], 'all');
SpekAmp = reshape(sqrt(Icam_speck), [], picnum);


%% TM retrieval via RAF21
disp('RAF21 for TM retrieval')

nWorker = 16;
optsRAF21 = struct;
optsRAF21.nWorker = nWorker;
optsRAF21.mu = 3;
optsRAF21.useGPU = 1;
optsRAF21.iters = 300;


subSpekAmp = SpekAmp(focInd, :);


% RAF2-1
[TM_raf21, ~, ~] = RAF21_tm(subSpekAmp.^2, SLMmsk, optsRAF21);


%% Focus scanning 
camSrc.ExposureTime = 11;

Infocus = reshape(TM_raf21', N, N, focNum);
CGH_focus = 128 * CudaLee_1920_1080(angle(Infocus), Lxy, alpha, beta, 0);  

picN = size(CGH_focus, 3);
FocSeqid = alpseqalloc(deviceid,bitplanes,picN);
alp_returnvalue = alpseqput(deviceid,FocSeqid,CGH_focus);

fprintf ('Focus CGHs into DMD done\n');

alp_returnvalue = alpseqtiming(deviceid,FocSeqid,picturetime,picturetime-100);
alp_returnvalue =  alpseqcontrol(deviceid,FocSeqid,10); %repeat sequence displaying
alp_returnvalue = alprojstart(deviceid,FocSeqid);
datav = getdata(camVid, picN);
Icam_focus = single(datav(:,:,1:picN));

clear datav;
delete(camVid)


clear camVid

fprintf ('Focus scanning done\n');

%% Show results
if ~ifScan
    figure('color', [1 1 1], 'position', [150 150 800 350]),
    subplot(121), imagesc(Icam_focus(:, :, 1));
    subplot(122), semilogy(Icam_focus(M/2, M/2-30:M/2+30)./max(Icam_focus(M/2, M/2-30:M/2+30))); title('Focus profile')  
    pbr = Icam_focus(focInd)/mean2(Icam_focus); 
    sgtitle(sprintf('exp-RAF21-singleFoc, N = %d, Max = %.2f, PBR=%.2f', N^2, max(max(Icam_focus(:, :))), pbr)); 
else
    Ifoc_raf21_max = max(Icam_focus, [], 3); 
    pbrs_raf21 = squeeze(max(Icam_focus,[], [1, 2])./mean(Icam_focus,[1, 2])); 
    unif_raf21 = calcUnif(Ifoc_raf21_max, focInd);

    PBRMap_raf21 = zeros(M, M); PBRMap_raf21(focInd) = pbrs_raf21;
    fig = figure('color','w', 'position', [120 150 620 500]);
    imagesc(PBRMap_raf21); colorbar
    title(sprintf('exp-RAF21-FocScan-PBRMap, PBR=%.2f, uniformity=%.2f%%', mean(pbrs_raf21), 100*unif_raf21), 'fontsize', 14)

    % save data
    d = datetime; d.Format = 'dd-MM-yyyy, HH-mm'; 
    foldName = sprintf('./Data/RAF21-MMF-FocScan/gamma=%d/%s', gamma, d);
    mkdir(foldName);
    save(sprintf('%s/exp_data.mat', foldName), 'Ifoc_raf21_max', 'PBRMap_raf21', 'pbrs_raf21', 'unif_raf21', 'focInd');
    saveas(fig, sprintf('%s/exp-RAF21-FocScan-PBRMap.fig', foldName)) 
end


fprintf ('TM retrieval experiment is done\n');
