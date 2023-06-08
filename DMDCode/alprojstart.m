
function alp_returnvalue = alprojstart(deviceid,seqid);
DEFAULT = int32(0);
alp_ok = int32(0);

alp_returnvalue = calllib('DMD','AlpProjStart',deviceid,seqid);
if alp_returnvalue ~= alp_ok;

   uiwait(msgbox(sprintf ('something wrong -- start proj.')));
end
end

