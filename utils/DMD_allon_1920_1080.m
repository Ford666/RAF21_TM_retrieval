% DMD all on

close all

%% close DMD    %close DMD at last
% clear TMseqid
% alp_returnvalue = dmdclose(deviceid);
% clear deviceid;

%% prepare
addpath('Z:\Shengfu\utils');
addpath('Z:\Shengfu\DMDCode');

% load DMD
if ~libisloaded('DMD')
    loadlibrary('C:\Program Files\ALP-4.3\ALP-4.3 API\x64\alp4395.dll','C:\Program Files\ALP-4.3\ALP-4.3 API\alp.h','alias','DMD');
end

%% DMD allocation or halt 
if ~exist('deviceid', 'var') %DMD allocation for the first exp
    deviceid = alpdevicealloc;
else                         %DMD halt for the next exp               
    alp_returnvalue = alpprojhalt(deviceid);
end


%% Parameters
bitplanes = 1;
picturetime = 33334; % us.
% picturetime = 4000; % us.
% picturetime = 2500; % us.
repeattimes = 999999;  % repeat should be stop by alpprojhalt 

N = 32; %Input size
Lxy = 28;
alpha = 1/4; beta = 1/4;


%% Load CGH onto DMD for TM calibration
ff = 128 * CudaLee_1920_1080(zeros(N, N, 'single'), Lxy, alpha, beta, 0);
[~,~,picnum] = size(ff);
AllOnSeqid = alpseqalloc(deviceid,bitplanes,picnum);
alp_returnvalue = alpseqput(deviceid,AllOnSeqid,ff);


alp_returnvalue = alpseqtiming(deviceid,AllOnSeqid,picturetime,picturetime-100);
alp_returnvalue =  alpseqcontrol(deviceid,AllOnSeqid,repeattimes);
alp_returnvalue = alprojstart(deviceid,AllOnSeqid);

fprintf ('DMD is all on\n');



