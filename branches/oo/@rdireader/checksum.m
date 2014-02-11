function valens=checksum(obj,headpos)
% Check validity of ensemble

fmes=fseek(obj.fid,headpos+2,-1);                                              %Go to byte after two header bytes
if fmes==-1                                                                %Check if move was possible
    valens=false;
    return
end
EnsBytes=fread(obj.fid,1,'uint16=>double');                                    %read number of bytes in ensemble
fmes=fseek(obj.fid,headpos+EnsBytes+2,-1);                                     %Go to the end of the ensemble to check if file doesn't end before
if fmes==-1
    valens=false;
    return
end
fseek(obj.fid,headpos,-1);                                                     %move to beginning of ensemble 
sumbytes=dec2bin(sum(fread(obj.fid,EnsBytes,'uint8')));                        %calculate checksum
checksum=fread(obj.fid,1,'*uint16');                                           %read checksum
if checksum==bin2dec(sumbytes(max((length(sumbytes)-15),1):end))                    %Check checksum
    valens=true;
else
    valens=false;
end
