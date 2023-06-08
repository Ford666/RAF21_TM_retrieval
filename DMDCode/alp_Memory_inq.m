function  memory_avil =alp_Memory_inq(deviceid)
    DEFAULT = int32(0);
    alp_ok = int32(0);
    memory_avil = int32(0);
    ALP_AVAIL_MEMORY = int32(2003);
    
    memory_avilptr = libpointer('int32Ptr', memory_avil);
    
    [alp_returnvalue memory_avil]= calllib('DMD','AlpDevInquire',deviceid,ALP_AVAIL_MEMORY,memory_avilptr);
    if alp_returnvalue ~= alp_ok
        uiwait(msgbox(sprintf ('something wrong -- alpseqfree.')));
    end
end