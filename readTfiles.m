function out = readTfiles(adcp,filenames)
%
%

%    Copyright 2013 Bart Vermeulen
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




%% parse input
inp=inputParser;     % Create an object of the InputParser class
inp.addRequired('adcp',@isstruct);                                       % Add required input inadcp
inp.addRequired('filenames',@(x) (iscellstr(x) | ischar(x)));           % Add the required variable 'mmeafilename' and check for right format
inp.parse(adcp,filenames);                                            % Parse input
if ischar(filenames)
    filenames=cellstr(filenames);                                    % Change character array to cell
end



clear inp

%% read data
% Initialize some variables
nFiles=numel(filenames);
out=struct();

% File header
out.note1=cell(1,nFiles);
out.note2=cell(1,nFiles);
out.depthCellLength=nan(1,nFiles);
out.blanking=nan(1,nFiles);
out.adcpDepth=nan(1,nFiles);
out.nbins=nan(1,nFiles);
out.pingsperens=nan(1,nFiles);
out.timeperens=nan(1,nFiles);
out.profilingMode=nan(1,nFiles);


out = initvars(out, 100);


csect=1;
lastens=csect;
for cf=1:numel(filenames)
    fid=fopen(filenames{cf});
    assert(fid~=-1);
    out.note1{cf}=fgetl(fid);
    out.note2{cf}=fgetl(fid);
    [dat,cnt]=fscanf(fid,'%d',7);
    if (7 == cnt)
        out.depthCellLength(cf)=dat(1);
        out.blanking(cf)=dat(2);
        out.adcpDepth(cf)=dat(3);
        out.nbins(cf)=dat(4);
        out.pingsperens(cf)=dat(5);
        out.timeperens(cf)=dat(6);
        out.profilingMode(cf)=dat(7);
        while ~(feof(fid))
            if csect > numel(out.pitch)
                out = initvars(out, csect * 2,nrows);
            end
            dat=fscanf(fid,'%d',9);
            if feof(fid), break, end
            out.timeV(csect,:)=dat(1:6);
            out.timeV(csect,6)=out.timeV(csect,6)+dat(7)/100;
            out.segNum(csect)=dat(8);
            out.ensInSeg(csect)=dat(9);
            dat=fscanf(fid,'%f',4);
            out.pitch(csect)=dat(1);
            out.roll(csect)=dat(2);
            out.corrHead(csect)=dat(3);
            out.adcpTemp(csect)=dat(4);
            dat=fscanf(fid,'%f',12);
            out.btvel(:,csect)=dat(1:4)';
            out.gpsvel(:,csect)=dat(5:8)';
            out.depth(:,csect)=dat(9:12)';
            dat=fscanf(fid,'%f',5);
            out.elapsDist(csect)=dat(1);
            out.elapsTime(csect)=dat(2);
            out.distNorth(csect)=dat(3);
            out.distEast(csect)=dat(4);
            out.distGood(csect)=dat(5);
            dat=fscanf(fid,'%f',5);
            out.lat(csect)=dat(1);
            out.long(csect)=dat(2);
            dat=fscanf(fid,'%f',9);
            out.qmid(csect)=dat(1);
            out.qtop(csect)=dat(2);
            out.qbot(csect)=dat(3);
            out.qstart(csect)=dat(4);
            out.distStart(csect)=dat(5);
            out.qend(csect)=dat(6);
            out.distEnd(csect)=dat(7);
            out.startDepthMid(csect)=dat(8);
            out.endDepthMid(csect)=dat(9);
            dat=textscan(fid,'%d %s %s %s %f %f',1);
            nrows=dat{1};
            out.lengthUnit(csect)=dat{2};
            out.velRef(csect)=dat(3);
            out.intensityUnit(csect)=dat(4);
            out.intensityScFact(csect)=dat{5};
            out.soundAbsorbtion(csect)=dat{6};
            dat=textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f',nrows);
            if isfield(out,'veldepth') && nrows > size(out.veldepth,1)
                pad=nan(nrows-size(out.veldepth,1),size(out.veldepth,2));
                out.veldepth=cat(1,out.veldepth,pad);
                out.velmag=cat(1,out.velmag,pad);
                out.veldir=cat(1,out.veldir,pad);
                out.pgood=cat(1,out.pgood,pad);
                out.Q=cat(1,out.Q,pad);
                pad=nan([size(pad) 4]);
                out.vel=cat(1,out.vel,pad);
                out.echo=cat(1,out.echo,pad);
            end % isfield out
            out.veldepth(1:nrows,csect)=dat{1};
            out.velmag(1:nrows,csect)=dat{2};
            out.veldir(1:nrows,csect)=dat{3};
            out.vel(1:nrows,csect,:)=cat(3,dat{4},dat{5},dat{6},dat{7});
            out.echo(1:nrows,csect,:)=cat(3,dat{8},dat{9},dat{10},dat{11});
            out.pgood(1:nrows,csect)=dat{12};
            out.Q(1:nrows,csect)=dat{13};
            csect=csect+1;
        end % while ~feof
    end % if cnt == 7
    fclose(fid);
    out.filenumber(lastens:csect-1)=cf;
    lastens=csect;
end

out = trimvars(out, lastens);

%% Try to align data if number of ensembles don't match based on time vector
if (isfield(out,'ensInSeg'))
    nens=numel(adcp.ensnum);
    if csect-1~=nens && all(out.ensInSeg==1)
        time1=datenum(adcp.timeV);
        timeV=out.timeV;
        if timeV(1,1)<50
            timeV(:,1)=timeV(:,1)+2000; 
        else 
            timeV(:,1)=timeV(:,1)+1900; 
        end
        time2=datenum(timeV);
        [~, i1, i2]=intersect(time1,time2);
        if isempty(i1)
            warning('readTfiles:noTimeOverlap','Could not find time overlap between transect file and raw ADCP file')
            return;
        end
        out2=out;
        out.timeV=nan(nens,6);
        out.timeV(i1,:)=out2.timeV(i2,:);
        allfields=fields(out);
        dim1flds=[11:16 20:35 39 40 48];
        for cfld=dim1flds
            out.(allfields{cfld})=nan(1,nens);
            out.(allfields{cfld})(i1)=out2.(allfields{cfld})(i2);
        end
        cell1fld=36:38;
        for cfld=cell1fld
            out.(allfields{cfld})=cell(1,nens);
            out.(allfields{cfld})(i1)=out2.(allfields{cfld})(i2);
        end
        dim2flds=[17:19 41:43 46 47];
        for cfld=dim2flds
            out.(allfields{cfld})=nan(size(out2.(allfields{cfld}),1),nens);
            out.(allfields{cfld})(:,i1)=out2.(allfields{cfld})(:,i2);
        end
        dim3flds=[44 45];
        for cfld=dim3flds
            out.(allfields{cfld})=nan(size(out2.(allfields{cfld}),1),nens,size(out2.(allfields{cfld}),3));
            out.(allfields{cfld})(:,i1,:)=out2.(allfields{cfld})(:,i2,:);
        end
    end
end % if ~ isemty out.ensInSeg
end

function out=initvars(out,init_sect,nrows)

if nargin < 3
    nrows = 1;
end
% initialize variables
out.timeV(init_sect,6)=nan;
out.segNum(init_sect)=nan;
out.ensInSeg(init_sect)=nan;
out.pitch(init_sect)=nan;
out.roll(init_sect)=nan;
out.corrHead(init_sect)=nan;
out.adcpTemp(init_sect)=nan;
out.btvel(4,init_sect)=nan;
out.gpsvel(4,init_sect)=nan;
out.depth(4,init_sect)=nan;
out.elapsDist(init_sect)=nan;
out.elapsTime(init_sect)=nan;
out.distNorth(init_sect)=nan;
out.distEast(init_sect)=nan;
out.distGood(init_sect)=nan;
out.lat(init_sect)=nan;
out.long(init_sect)=nan;
out.qmid(init_sect)=nan;
out.qtop(init_sect)=nan;
out.qbot(init_sect)=nan;
out.qstart(init_sect)=nan;
out.distStart(init_sect)=nan;
out.qend(init_sect)=nan;
out.distEnd(init_sect)=nan;
out.startDepthMid(init_sect)=nan;
out.endDepthMid(init_sect)=nan;
out.lengthUnit(init_sect)={[]};
out.velRef(init_sect)={[]};
out.intensityUnit(init_sect)={[]};
out.intensityScFact(init_sect)=nan;
out.soundAbsorbtion(init_sect)=nan;

out.veldepth(nrows,init_sect)=nan;
out.velmag(nrows,init_sect)=nan;
out.veldir(nrows,init_sect)=nan;
out.vel(nrows,init_sect,4)=nan;
out.echo(nrows,init_sect,4)=nan;
out.pgood(nrows,init_sect)=nan;
out.Q(nrows,init_sect)=nan;
end

function out=trimvars(out,lastens)

% initialize variables
out.timeV(lastens:end,:)=[];
out.segNum(lastens:end)=[];
out.ensInSeg(lastens:end)=[];
out.pitch(lastens:end)=[];
out.roll(lastens:end)=[];
out.corrHead(lastens:end)=[];
out.adcpTemp(lastens:end)=[];
out.btvel(:,lastens:end)=[];
out.gpsvel(:,lastens:end)=[];
out.depth(:,lastens:end)=[];
out.elapsDist(lastens:end)=[];
out.elapsTime(lastens:end)=[];
out.distNorth(lastens:end)=[];
out.distEast(lastens:end)=[];
out.distGood(lastens:end)=[];
out.lat(lastens:end)=[];
out.long(lastens:end)=[];
out.qmid(lastens:end)=[];
out.qtop(lastens:end)=[];
out.qbot(lastens:end)=[];
out.qstart(lastens:end)=[];
out.distStart(lastens:end)=[];
out.qend(lastens:end)=[];
out.distEnd(lastens:end)=[];
out.startDepthMid(lastens:end)=[];
out.endDepthMid(lastens:end)=[];
out.lengthUnit(lastens:end)=[];
out.velRef(lastens:end)=[];
out.intensityUnit(lastens:end)=[];
out.intensityScFact(lastens:end)=[];
out.soundAbsorbtion(lastens:end)=[];

out.veldepth(:,lastens:end)=[];
out.velmag(:,lastens:end)=[];
out.veldir(:,lastens:end)=[];
out.vel(:,lastens:end,:)=[];
out.echo(:,lastens:end,:)=[];
out.pgood(:,lastens:end)=[];
out.Q(:,lastens:end)=[];
end