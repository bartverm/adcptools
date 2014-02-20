function init(obj)
    reset_all(obj);


    %% Find valid ensembles
    if numel(obj.buf)<4, return, end
    
    % Search for valid ensemble headers
    obj.ens_pos=find([all([obj.buf(1:end-3)==127 obj.buf(2:end-2)==127],2); false])'; % find valid ensemble headers
    obj.ens_bytes=obj.parse_blocks(obj.ens_pos+2,'uint16');
    if isempty(obj.ens_pos), reset_all(obj),return, end
    
    % Remove if header position + ensbytes exceeds size of buffer
    fbad=obj.ens_pos+double(obj.ens_bytes)-1>numel(obj.buf)-2; 
    obj.ens_pos(fbad)=[]; 
    obj.ens_bytes(fbad)=[]; 
    if isempty(obj.ens_pos), reset_all(obj),return, end

    
    % Remove if checksum in file does not match with computed one
    csum=[0 cumsum(double(obj.buf))']; 
    csum=uint16(mod(csum(obj.ens_pos+double(obj.ens_bytes))-csum(obj.ens_pos),65536));
    csumr=typecast(reshape(obj.buf(bsxfun(@plus,obj.ens_pos+double(obj.ens_bytes),[0;1])),[],1),'uint16')';
    fbad=csumr~=csum;
    obj.ens_pos(fbad)=[]; 
    obj.ens_bytes(fbad)=[]; 
    if isempty(obj.ens_pos), reset_all(obj),return, end
    
    % Remove if ens_pos + ens_bytes +2 is higher than start of next
    % ensemble (Ensembles shouldn't overlap, but don't need to be
    % contiguous), and yes, this happens even after checksum checks.
    while numel(obj.ens_pos)>1 % loop while there's still more than 2 ensembles
        fbad=find((obj.ens_pos(1:end-1)+double(obj.ens_bytes(1:end-1))+2)>obj.ens_pos(2:end))+1; % Search for overlapping ensembles
        if isempty(fbad), break, end % if no bad ensemble is found, exit loop
        fbad=fbad([true; diff(fbad)~=1]); % remove only first of series of consecutive overlapping ensembles
        obj.ens_pos(fbad)=[]; % remove
        obj.ens_bytes(fbad)=[]; % remove
    end
    if isempty(obj.ens_pos), reset_all(obj),return, end

    
    obj.n_ensembles=numel(obj.ens_pos);
    
    %% Find position of data in buffer
    obj.ens_ndat=double(obj.buf(obj.ens_pos+5)');
    n_data=sum(obj.ens_ndat);
    obj.n_data=n_data;
    obj.data_ensid=zeros(n_data,1);
    obj.data_ensid(cumsum(obj.ens_ndat)+1)=1;
    obj.data_ensid=cumsum(obj.data_ensid)+1;
    obj.data_ensid(end)=[];
    pos_data_offset=(1:n_data);
    tmp_n_data_in_ens=[0 cumsum(obj.ens_ndat)];
    pos_data_offset=pos_data_offset-tmp_n_data_in_ens(obj.data_ensid);
    obj.data_offset=double(obj.parse_blocks(obj.ens_pos(obj.data_ensid)+6+(pos_data_offset-1)*2,'uint16'))+obj.ens_pos(obj.data_ensid);
    obj.data_headers=obj.parse_blocks(obj.data_offset,'uint16');
    
    %% Change general NMEA header to specific NMEA headers
    f_nmea=obj.data_headers==rdi.headers.NMEA_General;
    obj.data_headers(f_nmea)=obj.parse_blocks(obj.data_offset(f_nmea)+2,'uint16');
    
    
    %% Array sizes
    obj.max_nbins=double(max(obj.nbins));
    obj.max_nbeams=double(max(obj.usedbeams));
    
    %% Make array indexing
    ndat=double(obj.nbins).*double(obj.usedbeams);
    % Ensemble indices
    ensidx=zeros(1,sum(ndat)); 
    ensidx(cumsum(ndat(1:end-1))+1)=1;
    ensidx=cumsum(ensidx)+1;

    % Beam indices
    cnvels=cumsum([0  ndat]);
    beamidx=mod(cumsum(ones(1,sum(ndat)))-cnvels(ensidx)-1,double(obj.usedbeams(ensidx)))+1; % This is slow

    % Cell indices
    cncells=cumsum([0 double(obj.nbins)]);
    cellidx=cumsum(beamidx==1)-cncells(ensidx);

    obj.array_subs=[cellidx;ensidx;beamidx];
    obj.array_idx=sub2ind([double(max(obj.nbins)),obj.n_ensembles,double(max(obj.usedbeams))],cellidx,ensidx,beamidx);
    
end

%% reset function
function reset_all(obj)
    [obj.ens_pos,...
    obj.ens_bytes,...
    obj.ens_ndat,...
    obj.data_ensid,...
    obj.data_offset,...
    obj.data_headers]=deal([]);
end