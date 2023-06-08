
function  alp_returnvalue =  alpseqcontrol(deviceid,seqid,repeattimes, first_frame, last_frame)

    DEFAULT = int32(0);
    alp_ok = int32(0);
    
 if nargin < 4   
    ALP_SEQ_REPEAT = int32(2100);
    rt = int32(repeattimes);

    alp_returnvalue = calllib('DMD','AlpSeqControl',deviceid,seqid,ALP_SEQ_REPEAT,rt);
    if alp_returnvalue ~= alp_ok
        uiwait(msgbox(sprintf ('something wrong with seq_rep control.')));
    end
 end
 
 if nargin  == 5
     ALP_SEQ_REPEAT = int32(2100);
     rt = int32(repeattimes);
     alp_returnvalue = calllib('DMD','AlpSeqControl',deviceid,seqid,ALP_SEQ_REPEAT,rt);
     if alp_returnvalue ~= alp_ok
         uiwait(msgbox(sprintf ('something wrong with seq_rep control.')));
     end
     ALP_FIRSTFRAME = int32(2101);
     ALP_LASTFRAME = int32(2102);
     fr = int32(first_frame);
     lr = int32(last_frame);
     alp_returnvalue = calllib('DMD','AlpSeqControl',deviceid,seqid,ALP_FIRSTFRAME,fr);
     if alp_returnvalue ~= alp_ok
         uiwait(msgbox(sprintf ('something wrong with seq_frist_frame_control.')));
     end
          alp_returnvalue = calllib('DMD','AlpSeqControl',deviceid,seqid,ALP_LASTFRAME,lr);
     if alp_returnvalue ~= alp_ok
         uiwait(msgbox(sprintf ('something wrong with seq_last_frame_control.')));
     end
 end
 
     
     

  
    
end
