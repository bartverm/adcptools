function dataout=readADCP(files,varargin)
%READADCP Reads a RDI adcp binary data file
%         [ADCP]=readADCP(FILENAME) reads adcp data and outputs all data
%
%         FILENAME can be a character array with filenames on each row or a
%         cell of strings containing filenames in each element. The
%         function will read multiple files even when configuration differs
%         amongst them. In this case the function will issue a warning.
%         Settings will be stored once for each file
%
%         [ADCP]=readADCP(FILENAME,VARS) reads adcp data given in VARS
%         VARS should be a character array containing one of the following
%         characters:
%         e   ensemble info
%         v   velocity data
%         c   correlation data
%         h   echo intensity
%         p   percentage good
%         b   bottom track data
%         x   external data stored in the adcp file
%         f   external data stored in external winriver ascii-files (The
%             function will look for the files based on the raw-data
%             filenames. Being implemented!
%
%         Example: ADCP=readADCP('data_000.000; data_001.000','vEyx') will
%         read velocity, ensemble information and external data from the
%         two files
%
%         Author:      Bart Vermeulen, David Vermaas, Maximiliano Sassi
%         Last edit:   26-03-2009

%    Copyright 2007-2010 Bart Vermeulen, David Vermaas, Maximiliano Sassi
%
%    This file is part of ADCPTools.
%
%    ADCPTools is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    ADCPTools is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with ADCPTools.  If not, see <http://www.gnu.org/licenses/>.


% global variables are: data_defs, data_out, ens, data, raw_dat


%% Check input
if verLessThan('matlab','7.13')
    error(nargchk(1,2,nargin)) %#ok<NCHKN>
    error(nargoutchk(0,1,nargout)) %#ok<NCHKE>
else
    narginchk(1,2)
    nargoutchk(0,1)
end
inp=inputParser;                                                           % Create an object of the InputParser class
inp.addRequired('files',@(x) (iscellstr(x) | ischar(x)));                  % Add the required variable 'files' and check for right format
% inp.addOptional('flags','vhcpbex',@ischar);                                % Add the optional argument 'flags' and check for right format
inp.parse(files,varargin{:});                                              % Parse input

if ischar(files)
    files=cellstr(files);                                                % Change character array to cell
end
% flags=inp.Results.flags;
clear inp

%Check filenames
nfiles=numel(files);
raw_dat=cell(nfiles,1);
filesize=zeros(nfiles,1);
% disp('Reading files');
for cntfile=1:nfiles                                                       % Loop for files
    [fid,openmessage]=fopen(files{cntfile},'r');                           % Open file in read-only binary mode
    if fid==-1                                                             % If opening file fails
        warning('readADCP:WrongFile',[openmessage,': ',files{cntfile},...
            ', Skipping this file.'])                                      % Show warning
        continue                                                           % Continue to next file
    end
    [raw_dat{cntfile}, filesize(cntfile)]=fread(fid,'*uint8');             % read all data in file
    fclose(fid);                                                           % close file
end
raw_dat=cat(1,raw_dat{:});                                                 % concatenate raw data

% search for ensembles
ens.pos=find_valid_ensembles();
assert(~isempty(ens.pos),'readADCP:NoEnsemble','Could not find any valid ensemble')      % Generate error

% locate the data fields
ens.ndat=double(raw_dat(ens.pos+5)');
data=locate_data_fields();

% parse the fixed leader
array=struct(); % declare here to increase scope of this variable
dataout=struct();
read_data_fields()

n_nonroot_fields=numel(fieldnames(dataout))-1;
for fname=fieldnames(dataout.root)'
    dataout.(fname{1})=dataout.root.(fname{1});
end
dataout=rmfield(dataout,'root');
nfields=numel(fieldnames(dataout));
dataout=orderfields(dataout,[n_nonroot_fields+1:nfields 1:n_nonroot_fields]);

% END OF MAIN FUNCTION %%

% FUNCTIONS BELOW SHARE SCOPE WITH MAIN FUNCTIONS  
    function read_data_fields()
        data_defs=define_data_types();
        for count_data=1:numel(data_defs)
            tmpout=struct();
            switch data_defs(count_data).type
                case 'array'
                    % First chack if there is any of this data in the raw data
                    if ~any(data.headers==hex2dec(data_defs(count_data).header)), continue, end;
                    if ~init_array_indexing(data_defs(count_data).in_structure), continue, end; 
                    field_meta=data_defs(count_data).fields(1);
                    tmpout.(field_meta.name)=get_array_data(data_defs(count_data).header,field_meta, data_defs(count_data).in_structure);
                case 'fields'
                    fdata=data.headers==hex2dec(data_defs(count_data).header);
                    if ~any(fdata), continue, end
                    n_per_ens=accumarray(data.ensid(fdata)',ones(sum(fdata),1),[numel(ens.pos),1],@sum,0,false);
                    for count_field=1:numel(data_defs(count_data).fields)
                        field_meta=data_defs(count_data).fields(count_field);
                        tmpdat=reshape(parse_blocks(bsxfun(@plus,(0:field_meta.size-1)'* sizeof(field_meta.type),data.offset(fdata)+field_meta.offset-1),field_meta.type), field_meta.size, sum(fdata));
                        if any(n_per_ens>1)
                            out_dat=mat2cell(tmpdat,field_meta.size,n_per_ens);
                        else
                            nv=nullval(field_meta.type);
                            out_dat=nv(ones(field_meta.size,numel(ens.pos)));
                            out_dat(:,data.ensid(fdata))=tmpdat;
                        end
                        tmpout.(field_meta.name)=out_dat;
                    end
                otherwise
                    warning('readADCP:UnknownDataType',[data_defs(count_data).name,' data of unknown type: ', data_defs(count_data).type,', skipping.']);
                    continue
            end
            if isa(data_defs(count_data).postaction,'function_handle')
                tmpout=data_defs(count_data).postaction(tmpout);
            end
            assign_output(tmpout,data_defs(count_data))
        end
        
    end

    function assign_output(out_data,data_meta)
        if ~isfield(dataout,data_meta.in_structure)
            dataout.(data_meta.in_structure)=out_data;
        else
            for fieldname=fieldnames(out_data)'
                dataout.(data_meta.in_structure).(fieldname{1})=out_data.(fieldname{1});
            end
        end
    end



    function tf=init_array_indexing(struct_name)
        % Here we create indices that map the data as laid out
        % in the PD0 files to the shape of the output matrices
        tf=false;
        if isfield(array,struct_name), return, end
        if ~isfield(dataout.(struct_name),'nbins') || ~isfield(dataout.(struct_name),'usedbeams')
            warning('readADCP:NoVariableLeader','Cannot read array data without leader information')
            return
        end
        array.(struct_name).nbins=double(dataout.(struct_name).nbins); % fix for scalar nbins?
        array.(struct_name).nbins(array.(struct_name).nbins==nullval(class(dataout.(struct_name).nbins)))=0;
        array.(struct_name).usedbeams=double(dataout.(struct_name).usedbeams);
        array.(struct_name).usedbeams(array.(struct_name).nbins==nullval(class(dataout.(struct_name).usedbeams)))=0;
        ndat=array.(struct_name).nbins.*array.(struct_name).usedbeams;
        
        % Search for ensembles with data
        fgood=find(ndat~=0);
        ndat_good=ndat;
        ndat_good(ndat_good==0)=[];
        
        % Ensemble indices
        ensidx=zeros(1,sum(ndat_good));
        ensidx(cumsum(ndat_good(1:end-1))+1)=1;
        ensidx_fgood=cumsum(ensidx)+1;
        ensidx=fgood(ensidx_fgood);
        
        % Beam indices
        cnvels=cumsum([0  ndat_good]);
        beamidx=mod(cumsum(ones(1,sum(ndat_good)))-cnvels(ensidx_fgood)-1,double(dataout.(struct_name).usedbeams(ensidx)))+1; % This is slow
        
        % Cell indices
        cncells=cumsum([0 array.(struct_name).nbins(fgood)]);
        cellidx=cumsum(beamidx==1)-cncells(ensidx_fgood);
        
        % Output matrix linear indices
        array.(struct_name).idces=sub2ind([double(max(array.(struct_name).nbins)),numel(ens.pos),double(max(array.(struct_name).usedbeams))],cellidx,ensidx,beamidx);
        
        % Store results
        array.(struct_name).ensidx=ensidx;
        array.(struct_name).ndat=ndat;
        
        % Set output to true, since it seems everything went well
        tf=true;
    end

    function data_idx=find_last_datafield(header)
        data_idx=nan(1,numel(ens.pos));
        ftype=find(data.headers==hex2dec(header));
        data_idx(data.ensid(ftype))=ftype; % Keeps last, if multiple data blocks of same type. nan means data not found in ensemble
    end

    function val=get_array_data(header, field_meta, struct_name)
        nval=nullval(field_meta.type);
        nbins_max=max(array.(struct_name).nbins);
        nbeams_max=max(array.(struct_name).usedbeams);
        val=nval(ones(nbins_max,numel(ens.pos),nbeams_max)); % Initialize
        if isempty(val), return, end;
        ndat=array.(struct_name).ndat;        % locate array data
        pos=find_last_datafield(header); % Locate all  data
        if (all(isnan(pos))), return, end;
        
        % get indices 
        ensidx=array.(struct_name).ensidx;
        idces=array.(struct_name).idces;
        
        % compute offest in array field to each datum
        offset_in_array_field=cumsum([0 sizeof(field_meta.type)*ones(1,sum(ndat)-1)]); % pretend all velocities are behind each other
        cnvels=cumsum([0 ndat(1:end-1)]); % compute number of data in previous ensembles
        offset_in_array_field=offset_in_array_field-cnvels(ensidx)*sizeof(field_meta.type); % subtract them
        
        % account for possible ensembles without data (does this ever happen? assume yes for now)
        f_bad_ens=find(~isfinite(pos));
        f_bad_idx=ismember(ensidx,f_bad_ens); % Get index of all ensembles with missing data
        ensidx(f_bad_idx)=[];
        offset_in_array_field(f_bad_idx)=[];
        idces(f_bad_idx)=[];

        % compute position in raw data
        pos=data.offset(pos(ensidx))+2+offset_in_array_field; % This line is slow, any way to improve?
        val(idces)=parse_blocks(pos,field_meta.type);
    end

    function data=locate_data_fields()
        n_data=sum(ens.ndat);
        data.ensid=zeros(1,n_data);
        data.ensid(cumsum(ens.ndat)+1)=1;
        data.ensid=cumsum(data.ensid)+1;
        data.ensid(end)=[];
        pos_data_offset=(1:n_data);
        tmp_n_data_in_ens=[0 cumsum(ens.ndat)];
        pos_data_offset=pos_data_offset-tmp_n_data_in_ens(data.ensid);
        data.offset=double(parse_blocks(ens.pos(data.ensid)+6+(pos_data_offset-1)*2,'uint16'))+ens.pos(data.ensid);
        data.headers=parse_blocks(data.offset,'uint16');
        f_nmea=data.headers==8226;                                         % Replalce generic nmea header with specific nmea header
        data.headers(f_nmea)=parse_blocks(data.offset(f_nmea)+2,'uint16');
    end

    function ens_pos=find_valid_ensembles()
        ens_pos=[];
        if numel(raw_dat)<4, return, end
        
        % Search for valid ensemble headers
        ens_pos=find([all([raw_dat(1:end-3)==127 raw_dat(2:end-2)==127],2); false])'; % find valid ensemble headers
        ens_bytes=parse_blocks(ens_pos+2,'uint16');
        if isempty(ens_pos), return, end
        
        % Remove if header position + ensbytes exceeds size of buffer
        fbad=ens_pos+double(ens_bytes)-1>numel(raw_dat)-2;
        ens_pos(fbad)=[];
        ens_bytes(fbad)=[];
        if isempty(ens_pos), return, end
        
        
        % Remove if checksum in file does not match with computed one
        csum=[0 cumsum(double(raw_dat))'];
        csum=uint16(mod(csum(ens_pos+double(ens_bytes))-csum(ens_pos),65536));
        csumr=typecast(reshape(raw_dat(bsxfun(@plus,ens_pos+double(ens_bytes),[0;1])),[],1),'uint16')';
        fbad=csumr~=csum;
        ens_pos(fbad)=[];
        ens_bytes(fbad)=[];
        if isempty(ens_pos), return, end
        
        % Remove if ens_pos + ens_bytes +2 is higher than start of next
        % ensemble (Ensembles shouldn't overlap, but don't need to be
        % contiguous), and yes, this happens even after checksum checks.
        while numel(ens_pos)>1 % loop while there's still more than 2 ensembles
            fbad=find((ens_pos(1:end-1)+double(ens_bytes(1:end-1))+2)>ens_pos(2:end))+1; % Search for overlapping ensembles
            if isempty(fbad), break, end % if no bad ensemble is found, exit loop
            fbad=fbad([true; diff(fbad)~=1]); % remove only first of series of consecutive overlapping ensembles
            ens_pos(fbad)=[]; % remove
            ens_bytes(fbad)=[]; % remove
        end
    end

    function out=parse_blocks(pos,type)
        if isempty(pos), out = []; return, end
        
        if strcmp(type,'uint8')
            out=raw_dat(pos);
        elseif strcmp(type,'char')
            out=char(raw_dat(pos));
        else
            nb=sizeof(type);
            out=reshape(typecast(reshape(raw_dat(permute(bsxfun(@plus,pos,shiftdim(0:nb-1,-1)),[3 1 2])),[],1),type),size(pos));
        end
    end
end

% FUNCTIONS BELOW DO NOT SHARE VARIABLES WITH MAIN FUNCTION

function sysstr=sysinterp(sysid)
    sysstrtemp='';
    switch sysid(1:3)
        case '000'
            sysstrtemp=[sysstrtemp,'75 KHz System'];
        case '100'
            sysstrtemp=[sysstrtemp,'150 KHz System'];
        case '010'
            sysstrtemp=[sysstrtemp,'300 KHz System'];
        case '110'
            sysstrtemp=[sysstrtemp,'600 KHz System'];
        case '001'
            sysstrtemp=[sysstrtemp,'1200 KHz System'];
        case '101'
            sysstrtemp=[sysstrtemp,'2400 KHz System'];
    end
    if sysid(4)=='0'
        sysstrtemp=[sysstrtemp,', Concave beam pattern'];
    else
        sysstrtemp=[sysstrtemp,', Convex beam pattern'];
    end
    switch sysid(5:6)
        case '00'
            sysstrtemp=[sysstrtemp,', Sensor configuration #1'];
        case '10'
            sysstrtemp=[sysstrtemp,', Sensor configuration #2'];
        case '11'
            sysstrtemp=[sysstrtemp,', Sensor configuration #3'];
    end
    if sysid(7)=='0'
        sysstrtemp=[sysstrtemp,', XDCR HD Not attached'];
    else
        sysstrtemp=[sysstrtemp,', XDCR HD Attached'];
    end
    if sysid(8)=='0'
        sysstrtemp=[sysstrtemp,', Up facing beam'];
    else
        sysstrtemp=[sysstrtemp,', Down facing beam'];
    end
    switch sysid(9:10)
        case '00'
            sysstrtemp=[sysstrtemp,', Beam angle: 15E'];
        case '10'
            sysstrtemp=[sysstrtemp,', Beam angle: 20E'];
        case '11'
            sysstrtemp=[sysstrtemp,', Beam angle: 30E'];
        case '01'
            sysstrtemp=[sysstrtemp,', Unknown beam angle'];
    end
    switch sysid(13:16)
        case '0010'
            sysstrtemp=[sysstrtemp,', 4 Beam Janus configuration'];
        case '1010'
            sysstrtemp=[sysstrtemp,', 5 Beam Janus configuration, 3 Demods'];
        case '1111'
            sysstrtemp=[sysstrtemp,', 5 Beam Janus configuration, 2 Demods'];
    end
    sysstr=sysstrtemp;
end

% Function to interprete coordinate information (bits are reversed wrt manual of ADCP)
function corstr=corinterpr(corid)
    corstrtemp='';
    switch corid(4:5)
        case '00'
            corstrtemp=[corstrtemp,'Beam Coordinates Used'];
        case '10'
            corstrtemp=[corstrtemp,'Instrument Coordinates Used'];
        case '01'
            corstrtemp=[corstrtemp,'Ship Coordinates Used'];
        case '11'
            corstrtemp=[corstrtemp,'Earth Coordinates Used'];
    end
    if corid(1)=='1'
        corstrtemp=[corstrtemp,', Bin mapping used'];
    else
        corstrtemp=[corstrtemp,', Bin mapping not used'];
    end
    if corid(2)=='1'
        corstrtemp=[corstrtemp,', 3 beams solution allowed'];
    else
        corstrtemp=[corstrtemp,', 3 beams solution not allowed'];
    end
    if corid(3)=='1'
        corstrtemp=[corstrtemp,', Tilts used (only for ship or eart coordinates)'];
    else
        corstrtemp=[corstrtemp,', Tilts not used'];
    end
    corstr=corstrtemp;
end

function nb=sizeof(type)
    if strcmp(type,'char')
        nb=1;
        return
    end
    nb=cast(1,type); %#ok<NASGU>
    nb=whos('nb');
    nb=nb.bytes;
end

function val=nullval(type)
    if strcmp(type(1:3),'int'), val=intmin(type);
    elseif strcmp(type(1:4),'uint'), val = intmax(type);
    elseif any(strcmp(type,{'single','double'})), val=nan;
    elseif strcmp(type,'char'), val=' ';
    else val=0;
    end
end