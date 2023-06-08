
function alp_returnvalue = alpprojhalt(deviceid)
 
    alp_ok = int32(0);


    alp_returnvalue = calllib('DMD', 'AlpProjHalt', deviceid);

    if alp_returnvalue ~= alp_ok;

            uiwait(msgbox(sprintf ('something wrong with stop projection.')));
    end
end
