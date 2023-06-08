
function alp_returnvalue = alpseqput(deviceid,seqid,dmdmask)
DEFAULT = int32(0);
alp_ok = int32(0);


PicOffset = DEFAULT;
[~,~,PicLoad] = size(dmdmask); 
alp_returnvalue = calllib('DMD','AlpSeqPut',deviceid,seqid,PicOffset,PicLoad,dmdmask);
if alp_returnvalue ~= alp_ok;

            uiwait(msgbox(sprintf ('something wrong -- alpseqput.')));
end
end
