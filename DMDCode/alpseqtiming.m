
function  alp_returnvalue =  alpseqtiming(deviceid,seqid,picturetime,dptime)

    DEFAULT = int32(0);
    alp_ok = int32(0);

    PictureTime = picturetime;
    IlluminateTime = dptime;
    SynchDelay = DEFAULT;
    SynchPulseWidth = SynchDelay + IlluminateTime;
    TriggerInDelay = DEFAULT;
    alp_returnvalue = calllib('DMD','AlpSeqTiming',deviceid,seqid,IlluminateTime,PictureTime,SynchDelay,SynchPulseWidth,TriggerInDelay);
    if alp_returnvalue ~= alp_ok;
        uiwait(msgbox(sprintf ('something wrong with seqtiming control.')));
    end
    
end
