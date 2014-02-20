function out=parse_blocks(obj,pos,type)

if isempty(pos), out = []; return, end

if strcmp(type,'uint8')
    out=obj.buf(pos);
else
    nb=rdi.PD0.sizeof(type);
    out=typecast(reshape(obj.buf(bsxfun(@plus,pos,0:nb-1))',[],1),type);
end