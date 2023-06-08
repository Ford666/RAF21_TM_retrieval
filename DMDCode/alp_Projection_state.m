function  proj_state = alp_Projection_state(deviceid)
    DEFAULT = int32(0);
    alp_ok = int32(0);
    proj_state = int32(0);
    ALP_PROJ_STATE = int32(2400);
    ALP_PROJ_ACTIVE = int32(1200);
    ALP_PROJ_IDLE = int32(1201);
    proj_stateptr = libpointer('int32Ptr', proj_state);
    
    [alp_returnvalue, proj_state]= calllib('DMD','AlpProjInquire', deviceid, ALP_PROJ_STATE,proj_stateptr);
    if alp_returnvalue ~= alp_ok
        uiwait(msgbox(sprintf ('something wrong -- alpseqfree.')));
    end
    if proj_state == ALP_PROJ_ACTIVE
        proj_state = 1;
    elseif proj_state == ALP_PROJ_IDLE
        proj_state = 0;
    end
end