classdef Message < matlab.mixin.Heterogeneous & handle
% Generic NMEA messages class. This is an abstract class that is used to
% define NMEA messages and cannot be used directly, apart from the static
% method get_all_messages.
%
% Message properties (read only):
%   regular_expression - for a generic NMEA string
%   fields - returns the fields in the current message
%   format - returns the format for use with textscan
%   name - name of the message, same as message identifier in the message
%
% Message methods:
%   parse - parse a string containing nmea messages
%   checksum_is_valid - check validity of NMEA checksum
%
% Message methods (static):
%   get_all_messages - returns all available NMEA messages
%   parse_all - parse nmea data
%
% Message abstract methods:
%
%   name_static - defining the msg_id of the message
%   fields_static - returning the fields in the message
%  
%   see also:nmea, nmea.Field
%   Copyright 2020, Bart Vermeulen

%     This file is part of the NMEA toolbox.
% 
%     NMEA toolbox is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     NMEA toolbox is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with NMEA toolbox.  If not, see <https://www.gnu.org/licenses/>.

    properties(SetAccess = protected)
        talker_id_pattern (1,:) char = '[a-zA-Z]{2}'
        msg_id_pattern (1,:) char
        needs_valid_checksum (1,1) logical = true;
    end
    properties(Constant)
        start_pattern = "\$"
        start_field_pattern = "[a-zA-Z]+"
        content_pattern = "(?:,[^,\*\$\n]*)+"
        csum_pattern = "(\*[a-fA-F0-9]{2})?"
        generic_pattern = ...
                nmea.Message.start_pattern+...
                nmea.Message.start_field_pattern+...
                nmea.Message.content_pattern+...
                nmea.Message.csum_pattern;
    end
    properties(Dependent)
        message_pattern
    end
    properties(SetAccess=protected)       
% nmea.Message/fields read only property
%
%   returns an array of nmea.Field defining the data in the nmea message
%
%   see also: nmea, nmea.Message, nmea.Field
        fields
        
% nmea.Message/format read only property
%
%   returns a concatenation of the format strings defined in the fields, to
%   be used with textscan to extract data from NMEA messages.
%
%   see also: nmea, nmea.Message, nmea.Field, textscan
        format

% nmea.Message/name read only property
%
%   returns the name of the current message. Matches with the message
%   identiefier in the NMEA message. Three characater row vector.
%
%   see also: nmea, nmea.Message
        name

    end
    methods
        function val=get.format(obj)
            val=join([obj.fields.format]);
        end

        function val = get.message_pattern(obj)
            val = obj.start_pattern +...
                "(?<talker_id>" + obj.talker_id_pattern + ")" + ...
                obj.msg_id_pattern + ","+...
                join([obj.fields.named_pattern],",")+...
                obj.csum_pattern;
        end
    end
    methods (Access=protected)
        function out = post_process(obj,in)
            idx_dat=1;
            if isscalar(obj)
                for cf=1:numel(obj.fields)
                    out.(obj.fields(cf).name)=obj.fields(cf).post_process(in(idx_dat:idx_dat+obj.fields(cf).n_formats-1));
                    idx_dat=idx_dat+obj.fields(cf).n_formats;
                end
            else
                error('running post_process on object array not supported')
            end
        end
    end
    methods (Sealed)
        function val=parse(obj,str)
% function to parse NMEA data
%
%   dat=parse(obj,str) will parse the string data in str, using the
%   Messgages defined in the nmea.Message array obj. Returns a structure
%   containing one field with the name of the message, that holds a
%   structure with all the data read by the message.
%   Next to the message data each of the message structures will contain
%   the following additional fields:
%
%   talker_id - holding the talker ID of the nmea messages
%   char_pos - holding the starting position of the message in str. This
%   differs when str is a array of characters of a string array or a cell
%   of character arrays. See help of MATLAB builtin regexp function for an
%   explanation of its format.
%
%   see also: nmea, nmea.Message, get_all_messages, regexp
            
            val = struct;
            % generic NMEA part
            if ischar(str) && ~isrow(str)
                str = reshape(str,1,[]);
            end
            [matched,pos]=regexp(str,nmea.Message.generic_pattern,...
                'match','start'); % find NMEA strings
            if isempty(matched)
                return
            end
            csum_valid=nmea.Message.checksum_is_valid(matched); % compute checksum

            % process individual messages
            for co=1:numel(obj)
                tmp_dat=regexp(matched,obj(co).message_pattern,"names"); % find strings with right message id
                cur_msgs = ~cellfun(@isempty,tmp_dat);
                if obj(co).needs_valid_checksum
                    cur_msgs = cur_msgs & csum_valid;
                end
                tmp_dat=[tmp_dat{cur_msgs}];
                if isempty(tmp_dat)
                    continue
                end
                dat = struct;
                for cntf = 1:numel(obj(co).fields)
                    cf = obj(co).fields(cntf);
                    if isa(cf,'nmea.SkipField')
                        continue
                    end
                    if cf.format == ""
                        dat.(cf.name) = {tmp_dat.(cf.name)};
                    else
                        dat.(cf.name)=textscan(...
                            [strjoin({tmp_dat.(cf.name)},','),','],... % trailing comma for empty fields
                            join(cf.format),'Delimiter',',');
                        dat.(cf.name) = cf.post_process(dat.(cf.name));
                    end
                end
                dat.char_pos=reshape(pos(cur_msgs),[],1);
                val.(obj(co).name)=dat;
            end
        end 
    end
    methods (Static, Sealed)
        function val=checksum_is_valid(str)
% Checks validity of NMEA message checksum
%
%   val=checksum_is_valid(str) returns whether checksum of message in str
%   is correct. If a cell of strings or string array is given returns
%   validity of each message.
%
%   see also: nmea, nmea.Message
            if ischar(str)
                if numel(str)<4
                    val = false;
                    return
                end
                msgu = uint8(str(2:end-3));    % convert characters in string to numbers
                cs_calc=zeros(1,'uint8');
                for count = 1:length(msgu)       
                    cs_calc = bitxor(cs_calc,msgu(count));
                end
                try
                    val=hex2dec(str(end-1:end))==cs_calc;
                catch
                    val = false;
                end
            elseif iscellstr(str) || isstring(str)
                val=cellfun(@nmea.Message.checksum_is_valid,str);
            end
        end
        function msgs=get_all_messages()
% Returns array with all available NMEA message objects
%
%   msgs=get_all_messages() Returns array with all available NMEA messages.
%   see under 'see also' below, for supported messages.
%
%   see also: nmea, Message, GGAMessage, VTGMessage, HDTMessage,
%   ZDAMessage, GMPMessage, TROMessage, LINMessage, SPDMessage, ROTMessage,
%   INFMessage
            msgs=[nmea.GGAMessage;...
              nmea.VTGMessage;...
              nmea.HDTMessage;...
              nmea.ZDAMessage;...
              nmea.GMPMessage;...
              nmea.TROMessage;...
              nmea.LINMessage;...
              nmea.SPDMessage;...
              nmea.ROTMessage;...
              nmea.INFMessage;...
              nmea.RDENSMessage;...
              nmea.PSATHPRMessage;...
              nmea.RMCMessage;...
              nmea.DBTMessage;...
              nmea.DBKMessage;...
              nmea.DBSMessage];
        end
        function val=parse_all(str)
% Parse string data using all available messages
%
%   dat=nmea.Message.parse_all(str) reads all known messages in the given
%   string or cell of character arrays and returns the data available per
%   message.
%
%   see also: nmea, Message, parse
           msgs=nmea.Message.get_all_messages();
           val=msgs.parse(str);
        end
    end
end