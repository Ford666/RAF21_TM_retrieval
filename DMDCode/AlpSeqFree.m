
function alp_returnvalue = AlpSeqFree(deviceid,seqid)
DEFAULT = int32(0);
alp_ok = int32(0);
alp_returnvalue = calllib('DMD','AlpSeqFree',deviceid,seqid);
if alp_returnvalue ~= alp_ok

            uiwait(msgbox(sprintf ('something wrong -- alpseqfree.')));
end
end