classdef rdireader < reader & handle
    properties
        filename;
        fid;
    end
    methods
        function read(obj)
            obj.fid=fopen(obj.filename,'r','l'); % reading, in little-endian
            fpos=0;                                                                %Set file position pointer to beginning of file
            enscnt=0;
            while 1                                                                %Loop
                headpos=obj.search_head(fpos);                                      %Search for header
                if headpos==-1, break, end                                         %If no header was found, exit the loop
                if ~obj.checksum(headpos)                                          %If ensemble is not valid
                   warning('rdireader:CorruptEnsemble',['Invalid ensemble: ',num2str(enscnt+1)]) %Generate warning
                   fpos=headpos+1;                                                 %Point fpos to position after last header
                   continue                                                        %Restart loop
                end
                fseek(obj.fid,headpos+2,-1);                                           %Move to position to read ensemble size
                EnsBytes=fread(obj.fid,1,'uint16=>double');                            %read number of bytes in ensemble
                fpos=headpos+EnsBytes+2;                                           %Point to byte behind last read ensemble
                enscnt=enscnt+1;                                                   %Increase ensemble counter with one
                EnsStart(enscnt)=headpos;%#ok<AGROW>                               %Store starting position of ensemble %#ok<AGROW> 
                [NDataTypes{enscnt},DataOffset{enscnt},DataHeader{enscnt}]=...
                    obj.read_head(headpos);%#ok<AGROW>                              %Read info about where to find which data %#ok<AGROW> 
            end
            if isempty(EnsStart)                                                   %If no valid ensembles are found in current file
                warning('readadcp2:NoEnsembleInFile',...
                    ['No valid ensembles in file: ',obj.filename])               %Generate a warning
                fclose(obj.fid);                                                       %Close file
            else
                disp(['Found ',num2str(enscnt),...
                    ' valid ensembles in file: ',obj.filename]);                 %Display amount of valid ensembles in current file
            end

            pos_fl=EnsStart(1)+DataOffset{1}(DataHeader{1}(:,1)==0,1);
            fseek(obj.fid,pos_fl+8,-1);
            nbeams=fread(obj.fid,1,'uint8');
            nbins=fread(obj.fid,1,'uint8');
            nens=enscnt;
            obj.adcp=adcp(nbins,nens,nbeams);
            obj.read_fixed_leader(pos_fl);                         %Call function to read the fixed leader


            fclose(obj.fid);
            obj.status=reader_status.Valid;
        end
    end
    
    methods(Access=protected)
        headpos=search_head(obj,fpos)
        isvalid=checksum(obj,fpos)
        [ndat, offset,header]=read_head(obj,fpos)
        read_fixed_leader(obj,fpos)
    end
end