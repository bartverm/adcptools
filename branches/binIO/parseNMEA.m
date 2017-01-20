function data_out=parseNMEA(instr, varargin)


[standard, proprietary]=defineNMEA(); % load patterns to read

%% check input
P=inputParser; % create inpute parser
P.addRequired('instr',@(x) (iscellstr(x) | ischar(x))); % input should be cell of strings of character array
P.addOptional('string_id', false, @(x) (isscalar(x) && islogical(x)));
P.addParameter('Proprietary', true, @(x) (isscalar(x) && islogical(x)));
P.parse(instr, varargin{:});% parse input
string_id=P.Results.string_id;
read_prop=P.Results.Proprietary;

if ischar(instr) % if input is char
    instr=cellstr(instr); % convert to cell of char
end
n_in=numel(instr); % number of input strings

%% process standard NMEA strings
[strings,~]=regexp(cellstr(instr),...
    strcat(standard(1).fields.pattern),'names','split');                   % split nmea strings in talker id, msg id, message and checksum

fgood=find(cellfun(@numel,strings)==1);                                    % find non empty messages
strings=[strings{fgood}];                                                  % Keep only good messages
fgood_cs=arrayfun(@(x) csum([x.talker, x.id, x.content],...
    x.checksum),strings);                                                  % Check the checksum
fgood=fgood(fgood_cs);                                                     % remove ids of messages with bad checksum (needed for string_id)
strings=strings(fgood_cs);                                                 % Keep only messages with good checksum

ids_in=unique({strings(:).id});                                            % Get all unique message ids in input strings
for c_ids=1:numel(ids_in)                                                  % loop over the available message ids
    c_pat=standard(strcmpi({standard.name}, ids_in{c_ids}));               % get definition of current pattern
    if isempty(c_pat)                                                      % if no definition of current message is found
        warning(['don''t know how to read ', ids_in{c_ids}]);              % warn user 
        continue                                                           % continue with next message
    end
    fstrings=find(strcmpi({strings.id},c_pat.name));                       % find input strings with current message id
    [tmpdat,~]=regexp({strings(fstrings).content},...
        strcat(c_pat.fields.pattern),'names','split');                     % parse the current messages
    fempty=~cellfun(@isscalar, tmpdat);                                    % find strings that failed to parse (yes, this happens even after checking checksum)
    fstrings(fempty)=[];                                                   % remove these strings from further processing
    tmpdat=[tmpdat{~fempty}];                                              % all output structures can be the same and so can be concatenated
    for c_field=c_pat.fields                                               % loop over fields in current message
        tmpfield=textscan(strjoin({tmpdat.(c_field.name)},'\n'),...
            c_field.conversion);                           % convert field from char to right format
        if ~string_id
            data_out.(c_pat.name).(c_field.name)=...
                repmat(c_field.nullval,c_field.size,n_in);                     % initialize output field
            data_out.(c_pat.name).(c_field.name)(:,fgood(fstrings))=...
                horzcat(tmpfield{:})';                                         % store output
        else
            data_out.(c_pat.name).(c_field.name)=...
                horzcat(tmpfield{:})'; 
        end
    end
    if isa(c_pat.postaction,'function_handle')                             % if the message has a postaction function
        data_out.(c_pat.name)=c_pat.postaction(data_out.(c_pat.name));     % run the postaction
    end
    if string_id
        data_out.(c_pat.name).string_id=reshape(fgood(fstrings),1,[]);
    end
        
end

%% process proprietary NMEA messages

if read_prop
    for cnt_prop=1:numel(proprietary)
        c_prop=proprietary(cnt_prop);
        [strings,~]=regexp(cellstr(instr),...
            strcat(c_prop.head, c_prop.fields.pattern, c_prop.tail),...
            'names','split');                                                 % search for and read the current patterns
        fgood=find(cellfun(@numel,strings)==1);                                    % find non empty messages
        strings=[strings{fgood}];                                                  % Keep only good messages
        if c_prop.tail_is_checksum
            fgood_cs=cellfun(@(x) csum(x(2:end-2),x(end-1:end)),instr(fgood));                                                  % Check the checksum
            fgood=fgood(fgood_cs);                                                     % remove ids of messages with bad checksum (needed for string_id)
            strings=strings(fgood_cs);                                                 % Keep only messages with good checksum
        end
        for cnt_field=1:numel(c_prop.fields)                                 % loop over fields in current message
            c_field=c_prop.fields(cnt_field);
            tmpfield=textscan(strjoin({strings.(c_field.name)},'|'),...
                c_field.conversion,'delimiter','|');                           % convert field from char to right format
            if ~string_id
                data_out.(c_prop.name).(c_field.name)=...
                    repmat(c_field.nullval,c_field.size,n_in);                     % initialize output field
                data_out.(c_prop.name).(c_field.name)(:,fgood(fstrings))=...
                    horzcat(tmpfield{:})';                                         % store output
            else
                data_out.(c_prop.name).(c_field.name)=...
                    horzcat(tmpfield{:})'; 
            end
        end
        if isa(c_prop.postaction,'function_handle')                             % if the message has a postaction function
            data_out.(c_prop.name)=c_prop.postaction(data_out.(c_prop.name));     % run the postaction
        end
        if string_id
            data_out.(c_prop.name).string_id=reshape(fgood,1,[]);
        end        
    end
end

end


function tf=csum(string,sum)

msgu = uint8(string);    % convert characters in string to double values
cs_calc=zeros(1,'uint8');
for count = 1:length(msgu)       % checksum calculation ignores $ at start
    cs_calc = bitxor(cs_calc,msgu(count));  % checksum calculation
end
cs_calc = dec2hex(cs_calc,2);
tf=strcmp(cs_calc,sum);
end