function  deviceid = alpdevicealloc 
    DEFAULT = int32(0);
    alp_ok = int32(0);
    deviceid = int32(0);

    devicenumber = DEFAULT;
    initflag     = DEFAULT;
    deviceidptr = libpointer('uint32Ptr', deviceid);
    [alp_returnvalue, deviceid] = calllib('DMD', 'AlpDevAlloc', devicenumber, initflag, deviceidptr);
    if alp_returnvalue ~= alp_ok

        uiwait(msgbox(sprintf ('ALP device #%i can not be allocated, return value was %i.', deviceid, alp_returnvalue)));
    end
end


