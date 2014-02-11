function [NDataTypes,DataOffset,DataHeader]=read_head(obj,headpos)
fseek(obj.fid,headpos+5,-1);                                                   %Move to the 'number of data-types'-field
NDataTypes=fread(obj.fid,1,'uint8=>double');                                   %Read the number of data-types
DataOffset=fread(obj.fid,NDataTypes,'uint16');                                 %Read offsets to data-types
DataHeader=ones(NDataTypes,2,'uint16')*65535;                              %Initialize header vector
for cntDataType=1:NDataTypes
    fseek(obj.fid,headpos+DataOffset(cntDataType),-1);                         %Move to data header in ensemble
    DataHeader(cntDataType,:)=fread(obj.fid,2,'*uint16')';                     %read header of data type (And next byte for navigation data)
end