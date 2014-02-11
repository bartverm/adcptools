function headpos=search_head(obj,fpos)
% search in the file for a valid header
headpos=-1;
fmes=fseek(obj.fid,fpos,-1);
if fmes==-1
    return
end
while (~feof(obj.fid))&& headpos==-1                                           %Loop until a end of file or until a header is found
    head=fread(obj.fid,1,'*uint8');                                            %Read 1st byte
    if isempty(head) || head~=127                                          %Check for valid header value
        continue                                                           %If bad continue searching
    end
    head=fread(obj.fid,1,'*uint8');                                            %Read 2nd byte
    if isempty(head)|| head~=127                                           %Check for valid header value
        continue                                                           %If bad continue searching
    end
    headpos=ftell(obj.fid)-2;
end