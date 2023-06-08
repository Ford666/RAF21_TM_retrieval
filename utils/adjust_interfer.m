% Ensure the DMD is displaying the same CGH repeatedly and the shutter is
% open first, adjust the setting for interference with reference beam
% (e.g.,the angle, the polarizer) and capture the interference pattern
% repeatedly. Show the Fourier domain to monitor if the adjustment is good
% or not.

addpath('Z:\Shengfu\DMDCode');
close all


% load shutter
if ~libisloaded('PiUsb')
    setenv('MW_MINGW64_LOC','C:\TDM-GCC-64');
    loadlibrary('../usbshutter/PiUsb.dll','../usbshutter/PiUsb.h');
end

% shutter open
if ~exist('sptr', 'var')
    serialnumber = 0518;  %should be pre-determined! shutter1: 0520, shutter2: 0518
    [sptr, r] = calllib('PiUsb', 'piConnectShutter', 0, serialnumber);
end
calllib('PiUsb', 'piSetShutterState', 1, sptr);
pause(2)

%% Camera setting
imaqreset;
camVid = videoinput('mwspinnakerimaq', 1, 'Mono16');
camSrc = getselectedsource(camVid);

camSrc.ExposureTime = 198;
repatNum = 999999;

% M1 = 720; M2 = 540;
% xoffset = 0;
% yoffset = 0;
M1 = 160; M2 = 160;
xoffset = 176;
yoffset = 176;

if exist('xoffset', 'var') && exist('yoffset', 'var') 
    ROI_pos = [xoffset yoffset M1 M2];
else
    error('Please input camera ROI position!')
end
camVid = videoinput('mwspinnakerimaq', 1, 'Mono16');
camSrc.TriggerMode = 'Off';
camSrc.AdcBitDepth = 'Bit12';
camSrc.ExposureAuto = 'Off';
camSrc.ExposureTime = 198;
camSrc.Gain = 0;
camSrc.GammaEnable = 'false';

camVid.ROIPosition = ROI_pos;

%% Show spectrum of interfergram 

figure('color', 'w', 'position', [200, 400, 1000, 450]),
for xx = 1:repatNum
    inter = im2double(getsnapshot(camVid));
    [~,M,num] = size(inter);   
    Finters = fftshift(fftshift(fft2(inter), 1), 2); 
    Fh = abs(Finters(:, :, 1)).^2;
    subplot(121), imshow(inter, []); title('Real-time speckle pattern')
    subplot(122), imshow(Fh, [0,max(max(Fh))./1e5]), colormap(pink); title('Frequency spectrum')
    pause(0.5)
end
