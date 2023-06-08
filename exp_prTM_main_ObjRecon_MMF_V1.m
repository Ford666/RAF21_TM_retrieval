% Phase retrieval for complex TM recovery of a MMF
% Image transmission through MMF with the retrieved TM, at different
% sampling rates.

% Author: Shengfu Cheng
% Updated Date: 11/03/203


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
gamma_list = [3, 4, 5, 6, 7, 8]; %sampling ratio
gammaLen = length(gamma_list);
P = max(gamma_list) * N^2;
Lxy = 28;
alpha = 1/4; beta = 1/4;

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

    focInd = find((X-M/2).^2+(Y-M/2).^2 < radius^2);

else
    figure('color', [1 1 1], 'position', [150 150 550 480]),
    imagesc(uint8(255*Icam_allon)); title('Fiber alignment')
    colorbar('off'); colormap(gray); 

    focInd = find((X-M/2).^2+(Y-M/2).^2 < (3/8*M)^2);
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

trigNum = P + 1;
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

Icam_speck = single(datav(:, :, 1:picnum)); clearvars datav;
SpekAmp = reshape(sqrt(Icam_speck), [], picnum); clearvars Icam_speck;
SpekAmp = SpekAmp ./ max(SpekAmp,  [], 'all');

%% Image transmission through MMF 
if ~exist('guang', 'var'); load('./target/guang32x32.mat'); end
Object = guang;

SLMmsk_obj = exp(1i*pi*Object); 
CGH_obj = 128 * CudaLee_1920_1080(angle(SLMmsk_obj), Lxy, alpha, beta, 0);

FocSeqid = alpseqalloc(deviceid,bitplanes,1);
alp_returnvalue = alpseqput(deviceid,FocSeqid,CGH_obj);

fprintf ('Obj CGH into DMD done\n');

alp_returnvalue = alpseqtiming(deviceid,FocSeqid,picturetime,picturetime-100);
alp_returnvalue =  alpseqcontrol(deviceid,FocSeqid,repeattimes); %repeat sequence displaying
alp_returnvalue = alprojstart(deviceid,FocSeqid);
datav = getdata(camVid, 1);
Icam_obj = single(datav(:,:,1));
SpekAmp_obj = reshape(sqrt(Icam_obj(focInd)), [], 1);
SpekAmp_obj = SpekAmp_obj ./ max(SpekAmp_obj, [], 'all');

clearvars datav
delete(camVid)
clear camVid

fprintf ('Object measurement done\n');


%% TM retrieval via RAF21
disp('RAF21 for TM retrieval')

nWorker = 32;
optsRAF21 = struct;
optsRAF21.nWorker = nWorker;
optsRAF21.mu = 3;
optsRAF21.useGPU = 1;
optsRAF21.iters = 300;


subSpekAmp = SpekAmp(focInd, :); clearvars SpekAmp;

TM_raf21 = cell(1, gammaLen);

gammaInd = 1;
for gamma = gamma_list
    
    % RAF2-1
    P_temp = N^2 * gamma;
    TM_raf21{gammaInd} = RAF21_tm(subSpekAmp(:, 1:P_temp).^2, SLMmsk(:, 1:P_temp), optsRAF21);
    fprintf('gamma=%d done\n', gamma);
    pause(5)
    
    gammaInd = gammaInd+1;
end
% clearvars subSpekAmp;


%% RAF21 for object reconstruction
fprintf('Image transmission with one speckle pattern after TM scaling\n');
optsRAF_objRec = struct;
optsRAF_objRec.signalType = 'real';
optsRAF_objRec.iters = 100;
optsRAF_objRec.obj = Object; 
optsRAF_objRec.ratio = 2/3;


ObjRe_raf21 = cell(1, gammaLen);
ts_raf21 = cell(1, gammaLen);
errCurve_raf21 = cell(1, gammaLen); 
corrCurve_raf21 = cell(1, gammaLen); 
corr_raf21 = zeros(1, gammaLen);
gammaInd = 1;
for gamma = gamma_list
    
    % RAF2-1
    [ObjRe_raf21_temp, ts_raf21{gammaInd}, errCurve_raf21{gammaInd}, ...
        corrCurve_raf21{gammaInd}] = RAF_objRe(SpekAmp_obj.^2, TM_raf21{gammaInd}, optsRAF_objRec);
    
    ObjRe_raf21_temp = real(exp(-1i * angle(trace(Object(:)' * ObjRe_raf21_temp(:)))) * ObjRe_raf21_temp);
    ObjRe_raf21{gammaInd} = reshape((ObjRe_raf21_temp-min(ObjRe_raf21_temp))./(max(ObjRe_raf21_temp)-min(ObjRe_raf21_temp)), N, N);
    
    corr_raf21(gammaInd) = corr2(ObjRe_raf21{gammaInd}, Object);

    gammaInd = gammaInd+1;
end
fprintf ('Object reconstruction done\n');


%% Show results

% image
fig1 = figure('color', 'w', 'position', [50 50 1000 600]); 
gammaInd = 1;
for gamma = gamma_list
    subplot(2, 3, gammaInd), imagesc(ObjRe_raf21{gammaInd}), axis off, colormap(gray), 
    title(sprintf('\\gamma=%d, PCC=%.2f', gamma, corr_raf21(gammaInd))); 
    gammaInd = gammaInd+1;
end
h=colorbar('eastoutside','fontsize',15);
set(h,'Position', [0.92 0.10 0.02 0.85]);

% correlation curve
fig2 = figure('color', 'w', 'position', [150 150 600 550]); 
samInd = floor(linspace(1, optsRAF_objRec.iters, 15));
color_list = ["k", "#77AC30", "#0072BD", "#D95319", "#7E2F8E", "r"];
marker_list = ['x', 'd', 's', '^', 'v', 'o'];
for gammaInd = 1:gammaLen
    plot(1:optsRAF_objRec.iters, corrCurve_raf21{gammaInd}, 'Color', color_list(gammaInd), 'linewidth', 2.2, ...
        'Marker', marker_list(gammaInd), 'MarkerIndices', samInd, 'MarkerSize', 8); hold on
end
hold off
legend('\gamma=3', '\gamma=4', '\gamma=5', '\gamma=6', '\gamma=7', '\gamma=8', ...
    'location', 'northwest', 'fontsize', 14);
ylim([0 1]); set(gca,'FontSize',18, 'LineWidth', 1); 
xlabel('Iterations', 'fontsize', 21), ylabel('PCC', 'fontsize', 21);
RemoveSubplotWhiteArea(1, 1, 1, 1)


% error curve
maxErr = -1;
for gamma = gamma_list
    maxErr = max(maxErr, max(errCurve_raf21{gammaInd}, [], 'all'));
end

fig3 = figure('color', 'w', 'position', [150 150 600 550]); 
for gammaInd = 1:gammaLen
    plot(1:optsRAF_objRec.iters, errCurve_raf21{gammaInd}./maxErr, 'linewidth', 2); hold on
end
hold off
legend('\gamma=3', '\gamma=4', '\gamma=5', '\gamma=6', '\gamma=7', '\gamma=8', ...
    'location', 'northeast', 'fontsize', 12);
ylim([0 1]); grid on; set(gca,'FontSize',18, 'LineWidth', 1); 
xlabel('Iterations', 'fontsize', 18), ylabel('Normalzied error', 'fontsize', 18);
title(sprintf('RAF21| Error curve when N = %d', N^2), 'fontsize', 18);


% save data
d = datetime; d.Format = 'dd-MM-yyyy, HH-mm'; 
foldName = sprintf('./Data/ObjRec/RAF21/%s', d);
mkdir(foldName);
save(sprintf('%s/exp_data.mat', foldName), 'ObjRe_raf21', 'ts_raf21', 'errCurve_raf21', ...
    'corrCurve_raf21', 'corr_raf21', 'gamma_list', 'Object', 'N');
saveas(fig1, sprintf('%s/image.fig', foldName));
saveas(fig2, sprintf('%s/corrCurve.fig', foldName));
saveas(fig3, sprintf('%s/errCurve.fig', foldName));


fprintf ('TM retrieval experiment is done\n');
