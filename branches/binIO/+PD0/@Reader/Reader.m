classdef Reader < handle
    properties(Access=protected)
        raw_data
        ensemble
        data
        array
        dynamic_data_type_traits
        fixed_leader_traits
        variable_leader_traits
    end
    properties(Dependent, GetAccess=public, SetAccess=private)
        n_ensembles
        dynamic_traits
    end
    methods
        function val=get.n_ensembles(obj)
            val=numel(obj.ensemble.pos);
        end
        function load_data_from_file(obj,searchstr)
            files=dir(searchstr); % read all filenames
            assert(~isempty(files),'readDeployment:NoFileFound','Could not find any file'); % check there is at least one file with the given deployment name
            files([files.isdir])=[];
            files=cellfun(@fullfile, {files.folder}, {files.name},'UniformOutput', false);
            nfiles=numel(files);
            tmp_dat=cell(nfiles,1);
            filesize=zeros(nfiles,1);
            for cntfile=1:nfiles                                                       % Loop for files
                [fid,openmessage]=fopen(files{cntfile},'r');                           % Open file in read-only binary mode
                if fid==-1                                                             % If opening file fails
                    warning('readADCP:WrongFile',[openmessage,': ',files{cntfile},...
                        ', Skipping this file.'])                                      % Show warning
                    continue                                                           % Continue to next file
                end
                [tmp_dat{cntfile}, filesize(cntfile)]=fread(fid,'*uint8');             % read all data in file
                fclose(fid);                                                           % close file
            end
            obj.raw_data=cat(1,tmp_dat{:});
            obj.ensemble.pos=obj.find_valid_ensembles();
            assert(~isempty(obj.ensemble.pos),'readADCP:NoEnsemble','Could not find any valid ensemble')      % Generate error
            obj.ensemble.ndat=double(obj.raw_data(obj.ensemble.pos+5)');
            obj.data=obj.locate_data_blocks();
        end
        function obj=Reader(varargin)
            [obj.fixed_leader_traits, obj.variable_leader_traits, obj.dynamic_data_type_traits]=PD0.define_data_types();
        end
        function data=get_data(obj)
            data=PD0.Data();
            data.fixed_leader=obj.read_block(obj.fixed_leader_traits);
            data.variable_leader=obj.read_block(obj.variable_leader_traits);
            data_block_q=obj.dynamic_data_type_traits;
            while ~isempty(data_block_q)
                current_data_block_traits=data_block_q(1);
                switch class(current_data_block_traits)
                    case 'PD0.ArrayBlock_Traits'
                        vl=current_data_block_traits.nbins_block;          % Get name of the leader
                        if ~isprop(data,vl.name)                           % If leader has not been read yet
                            data_block_q(isequal(data_block_q,vl))=[];     % Put it on top of the que
                            data_block_q=[vl; data_block_q];               %#ok<AGROW> % Continue to item on top of que
                            continue
                        end
                        fl=current_data_block_traits.nbins_field;
                        if ~isprop(data.(vl.name),fl.name)
                            warning('cannot read data due to missing leader information')
                            data_block_q(1)=[];
                            continue
                        end
                        tmp_dat=obj.read_array(current_data_block_traits, data);
                    case 'PD0.DataBlock_Traits'
                        tmp_dat=obj.read_block(current_data_block_traits);
                end
                if ~isempty(tmp_dat)
                    propmeta=data.addprop(current_data_block_traits.name);
                    % set property characteristics
                    data.(current_data_block_traits.name)=tmp_dat;
                end
                data_block_q(1)=[];
            end
        end
    end
    
    methods (Access=protected)
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
        function array_block=read_array(obj,traits,data)
            % First chack if there is any of this data in the raw data
            if ~any(obj.data.headers==hex2dec(traits.header))
                return
            end
            if ~obj.init_array_indexing(traits,data)
                return
            end
            array_block=PD0.DataBlock(traits);
            array_block.(traits.fields.name)=obj.get_array_data(traits.header,traits.fields);
            if ~isempty(array_block) && ~isempty(traits.postaction)
                array_block=traits.postaction(array_block);
            end
        end
        function data_block=read_block(obj,traits)
            data_block=PD0.DataBlock.empty;
            fdata=obj.data.headers==hex2dec(traits.header);
            if ~any(fdata)
                return
            end
            data_block=PD0.DataBlock(traits);
            n_per_ens=accumarray(obj.data.ensid(fdata)',ones(sum(fdata),1),[numel(obj.ensemble.pos),1],@sum,0,false);
            for count_field=1:numel(traits.fields)
                field_meta=traits.fields(count_field);
                tmpdat=reshape(obj.parse_blocks(bsxfun(@plus,(0:field_meta.size-1)'* PD0.Reader.sizeof(field_meta.type),obj.data.offset(fdata)+field_meta.offset-1),field_meta.type), field_meta.size, sum(fdata));
                if any(n_per_ens>1)
                    out_dat=mat2cell(tmpdat,field_meta.size,n_per_ens);
                else
                    nv=PD0.Reader.nullval(field_meta.type);
                    out_dat=nv(ones(field_meta.size,numel(obj.ensemble.pos)));
                    out_dat(:,obj.data.ensid(fdata))=tmpdat;
                end
                data_block.(field_meta.name)=out_dat;
            end
            if ~isempty(data_block) && ~isempty(traits.postaction)
                data_block=traits.postaction(data_block);
            end
        end
        function ens_pos=find_valid_ensembles(obj)
            ens_pos=[];
            if numel(obj.raw_data)<4, return, end
            
            % Search for valid ensemble headers
            ens_pos=find([all([obj.raw_data(1:end-3)==127 obj.raw_data(2:end-2)==127],2); false])'; % find valid ensemble headers
            ens_bytes=obj.parse_blocks(ens_pos+2,'uint16');
            if isempty(ens_pos), return, end
            
            % Remove if header position + ensbytes exceeds size of buffer
            fbad=ens_pos+double(ens_bytes)-1>numel(obj.raw_data)-2;
            ens_pos(fbad)=[];
            ens_bytes(fbad)=[];
            if isempty(ens_pos), return, end
            
            
            % Remove if checksum in file does not match with computed one
            csum=[0 cumsum(double(obj.raw_data))'];
            csum=uint16(mod(csum(ens_pos+double(ens_bytes))-csum(ens_pos),65536));
            csumr=typecast(reshape(obj.raw_data(bsxfun(@plus,ens_pos+double(ens_bytes),[0;1])),[],1),'uint16')';
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
        function data=locate_data_blocks(obj)
            n_data=sum(obj.ensemble.ndat);
            data.ensid=zeros(1,n_data);
            data.ensid(cumsum(obj.ensemble.ndat)+1)=1;
            data.ensid=cumsum(data.ensid)+1;
            data.ensid(end)=[];
            pos_data_offset=(1:n_data);
            tmp_n_data_in_ens=[0 cumsum(obj.ensemble.ndat)];
            pos_data_offset=pos_data_offset-tmp_n_data_in_ens(data.ensid);
            data.offset=double(obj.parse_blocks(obj.ensemble.pos(data.ensid)+6+(pos_data_offset-1)*2,'uint16'))+obj.ensemble.pos(data.ensid);
            data.headers=obj.parse_blocks(data.offset,'uint16');
            f_nmea=data.headers==8226;                                         % Replalce generic nmea header with specific nmea header
            data.headers(f_nmea)=obj.parse_blocks(data.offset(f_nmea)+2,'uint16');
        end
        function out=parse_blocks(obj,pos,type)
            if isempty(pos), out = []; return, end
            
            if strcmp(type,'uint8')
                out=obj.raw_data(pos);
            elseif strcmp(type,'char')
                out=char(obj.raw_data(pos));
            else
                nb=PD0.Reader.sizeof(type);
                out=reshape(typecast(reshape(obj.raw_data(permute(bsxfun(@plus,pos,shiftdim(0:nb-1,-1)),[3 1 2])),[],1),type),size(pos));
            end
        end
        function tf=init_array_indexing(obj,traits,data)
            % Here we create indices that map the data as laid out
            % in the PD0 files to the shape of the output matrices
            tf=false;
            block=traits.nbins_block.name;
            field=traits.nbins_field.name;
            if isfield(obj.array,block)
                tf=true;
                return
            end
            
            obj.array.(block).nbins=double(data.(block).(field)); % fix for scalar nbins?
            obj.array.(block).nbins(obj.array.(block).nbins==PD0.Reader.nullval(class(data.(block).(field))))=0;
            obj.array.(block).usedbeams=traits.nbeams;
            nbeams=ones(1,obj.n_ensembles).*traits.nbeams;
            ndat=obj.array.(block).nbins.*obj.array.(block).usedbeams;
            
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
            % expand nbeams
            beamidx=mod(cumsum(ones(1,sum(ndat_good)))-cnvels(ensidx_fgood)-1,double(nbeams(ensidx)))+1; % This is slow
            
            % Cell indices
            cncells=cumsum([0 obj.array.(block).nbins(fgood)]);
            cellidx=cumsum(beamidx==1)-cncells(ensidx_fgood);
            
            % Output matrix linear indices
            obj.array.(block).idces=sub2ind([double(max(obj.array.(block).nbins)),obj.n_ensembles,traits.nbeams],cellidx,ensidx,beamidx);
            
            % Store results
            obj.array.(block).ensidx=ensidx;
            obj.array.(block).ndat=ndat;
            
            % Set output to true, since it seems everything went well
            tf=true;
        end
    end
    methods (Static, Access=protected)
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
            if strcmp(type(1:3),'int')
                val=intmin(type);
            elseif strcmp(type(1:4),'uint')
                val = intmax(type);
            elseif any(strcmp(type,{'single','double'}))
                val=nan;
            elseif strcmp(type,'char')
                val=' ';
            else
                val=0;
            end
        end
    end
end