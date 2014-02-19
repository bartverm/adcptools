function out=parse_blocks(obj,pos,type)

if isempty(pos), out = []; return, end

nb=obj.sizeof(type);
out=typecast(reshape(obj.buf(bsxfun(@plus,pos,0:nb-1))',[],1),type);