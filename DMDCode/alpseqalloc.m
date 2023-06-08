
function seqid = alpseqalloc(deviceid,bitplanes,picnum)
    DEFAULT = int32(0);
    alp_ok = int32(0);
    seqid = uint32(0);
    seqididptr = libpointer('uint32Ptr',seqid);
    [alp_returnvalue,seqid] = calllib('DMD', 'AlpSeqAlloc', deviceid,bitplanes,picnum,seqididptr);

    if alp_returnvalue ~= alp_ok
            uiwait(msgbox(sprintf ('something wrong with alloc seqid.')));
    end
end
