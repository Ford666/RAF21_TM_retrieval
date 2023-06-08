
function alp_returnvalue = dmdclose(deviceid)

DEFAULT = int32(0);
alp_ok = int32(0); 

    alp_returnvalue = calllib('DMD', 'AlpDevHalt', deviceid);
    if alp_returnvalue ~= alp_ok;
        uiwait(msgbox('Something wrong with halt DMD'));
    end
    
    alp_returnvalue = calllib('DMD', 'AlpDevFree', deviceid);
    if alp_returnvalue ~= alp_ok;
        uiwait(msgbox('Something wrong with free DMD'));
    end
end
