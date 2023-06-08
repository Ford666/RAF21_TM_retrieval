
function [BITPLANESq,PICNUMq,PICTURE_TIMEq,ILLUMINATE_TIMEq,SYNCH_DELAYq,SYNCH_PULSEWIDTHq,ON_TIMEq,OFF_TIMEq] = alpseqinq(deviceid,seqid)
DEFAULT = int32(0);
alp_ok = int32(0);
    
BITPLANES = int32(2200); % Bit depth of the pictures in the sequence
PICNUM = int32(2201);    % Number of pictures in the sequence
PICTURE_TIME = int32(2203); % Time between the start of consecutive pictures in the sequence in microseconds,
                            % the corresponding in frames per second is
                            % picture rate [fps] = 1 000 000 / ALP_PICTURE_TIME []
ILLUMINATE_TIME = int32(2204); % Duration of the display of one picture in microseconds
SYNCH_DELAY = int32(2205);   % Delay of the start of picture display with respect
                               % to the trigger output (master mode) in microseconds
SYNCH_PULSEWIDTH = int32(2206); % Duration of the active trigger output pulse in microseconds
% VD_DELAY = int32(2207); % Delay of the start of picture display with respect to the
%                         % active VD trigger input edge in microseconds
% MAX_TRIGGER_DELAY = int32(2209); % Maximal duration of trigger output to projection delay in microseconds
% MAX_VD_DELAY = int32(2210); % Maximal duration of trigger input to projection delay in microseconds

MIN_PICTURE_TIME = int32(2211);    % Minimum time between the start of consecutive pictures in microseconds
MIN_ILLUMINATE_TIME = int32(2212); % Minimum duration of the display of one picture in microseconds
                                   % depends on ALP_BITNUM and ALP_BIN_MODE
MAX_PICTURE_TIME = int32(2213);    % Maximum value of ALP_PICTURE_TIME
                                   % ALP_PICTURE_TIME = ALP_ON_TIME + ALP_OFF_TIME
                                   % ALP_ON_TIME may be smaller than ALP_ILLUMINATE_TIME
ON_TIME = int32(2214);  % Total active projection time
OFF_TIME = int32(2215); % Total inactive projection time
    
% BITPLANES
    BITPLANESq = int32(0);
    INQUIRETR = libpointer('int32Ptr',BITPLANESq);
    [alp_returnvalue,BITPLANESq] = calllib('DMD', 'AlpSeqInquire', deviceid,seqid,BITPLANES,INQUIRETR);

    if alp_returnvalue ~= alp_ok;

            uiwait(msgbox(sprintf ('something wrong WHEN Seqinquire.')));
    end
    
% PICNUM
    PICNUMq = int32(0);
    INQUIRETR = libpointer('int32Ptr',PICNUMq);
    [alp_returnvalue,PICNUMq] = calllib('DMD', 'AlpSeqInquire', deviceid,seqid,PICNUM,INQUIRETR);

    if alp_returnvalue ~= alp_ok;

            uiwait(msgbox(sprintf ('something wrong WHEN Seqinquire.')));
    end
    
% PICTURE_TIME
    PICTURE_TIMEq = int32(0);
    INQUIRETR = libpointer('int32Ptr',PICTURE_TIMEq);
    [alp_returnvalue,PICTURE_TIMEq] = calllib('DMD', 'AlpSeqInquire', deviceid,seqid,PICTURE_TIME,INQUIRETR);

    if alp_returnvalue ~= alp_ok;

            uiwait(msgbox(sprintf ('something wrong WHEN Seqinquire.')));
    end
% ILLUMINATE_TIME
    ILLUMINATE_TIMEq = int32(0);
    INQUIRETR = libpointer('int32Ptr',ILLUMINATE_TIMEq);
    [alp_returnvalue,ILLUMINATE_TIMEq] = calllib('DMD', 'AlpSeqInquire', deviceid,seqid,ILLUMINATE_TIME,INQUIRETR);

    if alp_returnvalue ~= alp_ok;

            uiwait(msgbox(sprintf ('something wrong WHEN Seqinquire.')));
    end
% SYNCH_DELAY
    SYNCH_DELAYq = int32(0);
    INQUIRETR = libpointer('int32Ptr',SYNCH_DELAYq);
    [alp_returnvalue,SYNCH_DELAYq] = calllib('DMD', 'AlpSeqInquire', deviceid,seqid,SYNCH_DELAY,INQUIRETR);

    if alp_returnvalue ~= alp_ok;

            uiwait(msgbox(sprintf ('something wrong WHEN Seqinquire.')));
    end


% SYNCH_PULSEWIDTH
    SYNCH_PULSEWIDTHq = int32(0);
    INQUIRETR = libpointer('int32Ptr',SYNCH_PULSEWIDTHq);
    [alp_returnvalue,SYNCH_PULSEWIDTHq] = calllib('DMD', 'AlpSeqInquire', deviceid,seqid,SYNCH_PULSEWIDTH,INQUIRETR);

    if alp_returnvalue ~= alp_ok;

            uiwait(msgbox(sprintf ('something wrong WHEN Seqinquire.')));
    end

% ON_TIME
    ON_TIMEq = int32(0);
    INQUIRETR = libpointer('int32Ptr',ON_TIMEq);
    [alp_returnvalue,ON_TIMEq] = calllib('DMD', 'AlpSeqInquire', deviceid,seqid,ON_TIME,INQUIRETR);

    if alp_returnvalue ~= alp_ok;

            uiwait(msgbox(sprintf ('something wrong WHEN Seqinquire.')));
    end
% OFF_TIME
    OFF_TIMEq = int32(0);
    INQUIRETR = libpointer('int32Ptr',OFF_TIMEq);
    [alp_returnvalue,OFF_TIMEq] = calllib('DMD', 'AlpSeqInquire', deviceid,seqid,OFF_TIME,INQUIRETR);

    if alp_returnvalue ~= alp_ok;

            uiwait(msgbox(sprintf ('something wrong WHEN Seqinquire.')));    
    end
    
    
end
