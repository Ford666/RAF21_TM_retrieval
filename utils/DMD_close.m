% DMD close

close all

%% close DMD     %close DMD at last
if exist('deviceid', 'var')
    alp_returnvalue = dmdclose(deviceid);
    clear deviceid;
    if exist('TMSeqid', 'var')
        clear TMSeqid
    end
    if exist('FocSeqid', 'var')
        clear FocSeqid
    end
    if exist('AllOnSeqid', 'var')
        clear AllOnSeqid
    end
    if exist('AllSeqid', 'var')
        clear AllSeqid
    end
end