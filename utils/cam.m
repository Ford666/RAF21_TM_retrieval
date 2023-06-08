function out = cam(N)
%CAM Code for creating a video input object.
%   
%   This is the machine generated representation of a video input object.
%   This MATLAB code file, CAM.M, was generated from the OBJ2MFILE function.
%   A MAT-file is created if the object's UserData property is not 
%   empty or if any of the callback properties are set to a cell array  
%   or to a function handle. The MAT-file will have the same name as the 
%   code file but with a .MAT extension. To recreate this video input object,
%   type the name of the code file, cam, at the MATLAB command prompt.
%   
%   The code file, CAM.M and its associated MAT-file, CAM.MAT (if
%   it exists) must be on your MATLAB path.
%   
%   Example: 
%       vidobj = cam;
%   
%   See also VIDEOINPUT, IMAQDEVICE/PROPINFO, IMAQHELP, PATH.
%   

% Check if we can check out a license for the Image Acquisition Toolbox.
canCheckoutLicense = license('checkout', 'Image_Acquisition_Toolbox');

% Check if the Image Acquisition Toolbox is installed.
isToolboxInstalled = exist('videoinput', 'file');

if ~(canCheckoutLicense && isToolboxInstalled)
    % Toolbox could not be checked out or toolbox is not installed.
    error(message('imaq:obj2mfile:invalidToolbox'));
end

% Load the MAT-file containing UserData and CallBack property values.
try
    MATvar = load('Z:\tianting\XGA\testcamera\cam');
    MATLoaded = true;
catch
    warning(message('imaq:obj2mfile:MATload'));
   MATLoaded = false;
end


% Device Properties.
adaptorName = 'mwspinnakerimaq';
deviceID = 1;
vidFormat = 'Mono16';
tag = '';

% Search for existing video input objects.
existingObjs1 = imaqfind('DeviceID', deviceID, 'VideoFormat', vidFormat, 'Tag', tag);

if isempty(existingObjs1)
    % If there are no existing video input objects, construct the object.
    vidObj1 = videoinput(adaptorName, deviceID, vidFormat);
else
    % There are existing video input objects in memory that have the same
    % DeviceID, VideoFormat, and Tag property values as the object we are
    % recreating. If any of those objects contains the same AdaptorName
    % value as the object being recreated, then we will reuse the object.
    % If more than one existing video input object contains that
    % AdaptorName value, then the first object found will be reused. If
    % there are no existing objects with the AdaptorName value, then the
    % video input object will be created.

    % Query through each existing object and check that their adaptor name
    % matches the adaptor name of the object being recreated.
    for i = 1:length(existingObjs1)
        % Get the object's device information.
        objhwinfo = imaqhwinfo(existingObjs1{i});
        % Compare the object's AdaptorName value with the AdaptorName value
        % being recreated.
        if strcmp(objhwinfo.AdaptorName, adaptorName)
            % The existing object has the same AdaptorName value as the
            % object being recreated. So reuse the object.
            vidObj1 = existingObjs1{i};
            % There is no need to check the rest of existing objects.
            % Break out of FOR loop.
            break;
        elseif(i == length(existingObjs1))
            % We have queried through all existing objects and no
            % AdaptorName values matches the AdaptorName value of the
            % object being recreated. So the object must be created.
            vidObj1 = videoinput(adaptorName, deviceID, vidFormat);
        end %if
    end %for
end %if

% Configure properties whose values are saved in Z:\tianting\XGA\testcamera\cam.mat.
if (MATLoaded)
    % MAT-file loaded successfully. Configure the properties whose values
    % are saved in the MAT-file.
    set(vidObj1, 'ErrorFcn', MATvar.errorfcn1);
else
   % MAT-file could not be loaded. Configure properties whose values were
   % saved in the MAT-file to their default value.
    set(vidObj1, 'ErrorFcn', @imaqcallback);
end

% Configure vidObj1 properties.
set(vidObj1, 'FramesPerTrigger', 1);
% set(vidObj1, 'ROIPosition', [263 199 200 200]);
set(vidObj1, 'TriggerRepeat', N-1);
set(vidObj1, 'Timeout', 200); %N/200


% Configure vidObj1 triggering.
triggerconfig(vidObj1, 'hardware', 'DeviceSpecific', 'DeviceSpecific');

% Configure vidObj1's video source properties.
srcObj1 = get(vidObj1, 'Source');
% set(srcObj1(1), 'AasRoiEnable', '(Currently not accessible)');
% set(srcObj1(1), 'AasRoiHeight', '(Currently not accessible)');
% set(srcObj1(1), 'AasRoiOffsetX', '(Currently not accessible)');
% set(srcObj1(1), 'AasRoiOffsetY', '(Currently not accessible)');
% set(srcObj1(1), 'AasRoiWidth', '(Currently not accessible)');
% set(srcObj1(1), 'AcquisitionFrameRate', 766.622);
% set(srcObj1(1), 'AcquisitionResultingFrameRate', 760.389);
% set(srcObj1(1), 'AutoExposureControlPriority', '(Currently not accessible)');
% set(srcObj1(1), 'AutoExposureTargetGreyValue', '(Currently not accessible)');
% set(srcObj1(1), 'BalanceRatio', '(Currently not accessible)');
% set(srcObj1(1), 'BalanceRatioSelector', '(Currently not accessible)');
% set(srcObj1(1), 'BalanceWhiteAuto', '(Currently not accessible)');
% set(srcObj1(1), 'BalanceWhiteAutoDamping', '(Currently not accessible)');
% set(srcObj1(1), 'BalanceWhiteAutoLowerLimit', '(Currently not accessible)');
% set(srcObj1(1), 'BalanceWhiteAutoProfile', '(Currently not accessible)');
% set(srcObj1(1), 'BalanceWhiteAutoUpperLimit', '(Currently not accessible)');
% set(srcObj1(1), 'BlackLevelRaw', '(Currently not accessible)');
% set(srcObj1(1), 'ColorTransformationEnable', '(Currently not accessible)');
% set(srcObj1(1), 'ColorTransformationSelector', '(Currently not accessible)');
% set(srcObj1(1), 'ColorTransformationValue', '(Currently not accessible)');
% set(srcObj1(1), 'ColorTransformationValueSelector', '(Currently not accessible)');
% set(srcObj1(1), 'CounterEventActivation', '(Currently not accessible)');
% set(srcObj1(1), 'CounterResetActivation', '(Currently not accessible)');
% set(srcObj1(1), 'CounterResetSource', '(Currently not accessible)');
% set(srcObj1(1), 'CounterValueAtReset', '(Currently not accessible)');
% set(srcObj1(1), 'DeviceLinkCurrentThroughput', 6.0898e+07);
% set(srcObj1(1), 'DeviceMaxThroughput', 6.13972e+07);
% set(srcObj1(1), 'DeviceTemperature', 48.25);
% set(srcObj1(1), 'DeviceUptime', 73603);
% set(srcObj1(1), 'EventErrorCode', '(Currently not accessible)');
% set(srcObj1(1), 'EventErrorFrameID', '(Currently not accessible)');
% set(srcObj1(1), 'EventErrorTimestamp', '(Currently not accessible)');
% set(srcObj1(1), 'EventExposureEndFrameID', '(Currently not accessible)');
% set(srcObj1(1), 'EventExposureEndTimestamp', '(Currently not accessible)');
% set(srcObj1(1), 'EventSerialData', '(Currently not accessible)');
% set(srcObj1(1), 'EventSerialDataLength', '(Currently not accessible)');
% set(srcObj1(1), 'EventSerialPortReceiveTimestamp', '(Currently not accessible)');
% set(srcObj1(1), 'EventSerialReceiveOverflow', '(Currently not accessible)');
% set(srcObj1(1), 'EventTestTimestamp', '(Currently not accessible)');
% set(srcObj1(1), 'ExposureActiveMode', '(Currently not accessible)');
% set(srcObj1(1), 'FfcEnable', '(Currently not accessible)');
% set(srcObj1(1), 'FfcMode', '(Currently not accessible)');
% set(srcObj1(1), 'FfcUserGain', '(Currently not accessible)');
% set(srcObj1(1), 'FfcUserOffset', '(Currently not accessible)');
% set(srcObj1(1), 'FfcUserTableXCoordinate', '(Currently not accessible)');
% set(srcObj1(1), 'Gamma', '(Currently not accessible)');
% set(srcObj1(1), 'LineStatusAll', 12);
% set(srcObj1(1), 'LinkUptime', 73580);
% set(srcObj1(1), 'PixelDynamicRangeMax', 65535);
% set(srcObj1(1), 'PowerSupplyCurrent', 0.449463);
% set(srcObj1(1), 'RgbTransformLightSource', '(Currently not accessible)');
% set(srcObj1(1), 'Saturation', '(Currently not accessible)');
% set(srcObj1(1), 'SaturationEnable', '(Currently not accessible)');
% set(srcObj1(1), 'SequencerSetActive', '(Currently not accessible)');
% set(srcObj1(1), 'SequencerSetNext', '(Currently not accessible)');
% set(srcObj1(1), 'SequencerTriggerActivation', '(Currently not accessible)');
% set(srcObj1(1), 'Sharpening', '(Currently not accessible)');
% set(srcObj1(1), 'SharpeningAuto', '(Currently not accessible)');
% set(srcObj1(1), 'SharpeningEnable', '(Currently not accessible)');
% set(srcObj1(1), 'SharpeningThreshold', '(Currently not accessible)');
% set(srcObj1(1), 'TransferBlockCount', '(Currently not accessible)');
% set(srcObj1(1), 'TransferOperationMode', '(Currently not accessible)');
% set(srcObj1(1), 'TransferQueueMode', '(Currently not accessible)');
% set(srcObj1(1), 'V3_3Enable', '(Currently not accessible)');


out = vidObj1 ;
