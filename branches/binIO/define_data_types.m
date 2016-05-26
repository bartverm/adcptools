function datadefs=define_data_types()

%%% FIXED LEADER %%%
fl.name='Fixed Leader';
fl.header=0;
fl.type='fields';
fl.fields(1 ).name='firmver';           fl.fields(1 ).size=1;    fl.fields(1 ).type='uint8';   fl.fields(1 ).offset=3;
fl.fields(2 ).name='firmrev';           fl.fields(2 ).size=1;    fl.fields(2 ).type='uint8';   fl.fields(2 ).offset=4;
fl.fields(3 ).name='sysconf';           fl.fields(3 ).size=1;    fl.fields(3 ).type='uint16';   fl.fields(3 ).offset=5;
fl.fields(4 ).name='SymData';           fl.fields(4 ).size=1;    fl.fields(4 ).type='uint8';   fl.fields(4 ).offset=7;
fl.fields(5 ).name='LagLength';         fl.fields(5 ).size=1;    fl.fields(5 ).type='uint8';   fl.fields(5 ).offset=8;
fl.fields(6 ).name='usedbeams';         fl.fields(6 ).size=1;    fl.fields(6 ).type='uint8';   fl.fields(6 ).offset=9;
fl.fields(7 ).name='nbins';             fl.fields(7 ).size=1;    fl.fields(7 ).type='uint8';   fl.fields(7 ).offset=10;
fl.fields(8 ).name='pingperens';        fl.fields(8 ).size=1;    fl.fields(8 ).type='uint16';  fl.fields(8 ).offset=11;
fl.fields(9 ).name='binsize';           fl.fields(9 ).size=1;    fl.fields(9 ).type='uint16';  fl.fields(9 ).offset=13;
fl.fields(10).name='blnk';              fl.fields(10).size=1;    fl.fields(10).type='uint16';  fl.fields(10).offset=15;
fl.fields(11).name='mode';              fl.fields(11).size=1;    fl.fields(11).type='uint8';   fl.fields(11).offset=17;
fl.fields(11).name='minthrsh';          fl.fields(11).size=1;    fl.fields(11).type='uint8';   fl.fields(11).offset=18;
fl.fields(12).name='ncodrep';           fl.fields(12).size=1;    fl.fields(12).type='uint8';   fl.fields(12).offset=19;
fl.fields(13).name='minpercgood';       fl.fields(13).size=1;    fl.fields(13).type='uint8';   fl.fields(13).offset=20;
fl.fields(14).name='maxerrvel';         fl.fields(14).size=1;    fl.fields(14).type='uint16';  fl.fields(14).offset=21;
fl.fields(15).name='Tbetweenpng_min';   fl.fields(15).size=1;    fl.fields(15).type='uint8';   fl.fields(15).offset=23;
fl.fields(16).name='Tbetweenpng_sec';   fl.fields(16).size=1;    fl.fields(16).type='uint8';   fl.fields(16).offset=24;
fl.fields(17).name='Tbetweenpng_cs';    fl.fields(17).size=1;    fl.fields(17).type='uint8';   fl.fields(17).offset=25;
fl.fields(18).name='corinfo';           fl.fields(18).size=1;    fl.fields(18).type='uint8';   fl.fields(18).offset=26;
fl.fields(19).name='headalign';         fl.fields(19).size=1;    fl.fields(19).type='uint16';  fl.fields(19).offset=27;
fl.fields(20).name='headbias';          fl.fields(20).size=1;    fl.fields(20).type='uint16';  fl.fields(20).offset=29;
fl.fields(21).name='sensource';         fl.fields(21).size=1;    fl.fields(21).type='uint8';   fl.fields(21).offset=31;
fl.fields(22).name='senavail';          fl.fields(22).size=1;    fl.fields(22).type='uint8';   fl.fields(22).offset=32;
fl.fields(23).name='distmidbin1';       fl.fields(23).size=1;    fl.fields(23).type='uint16';  fl.fields(23).offset=33;
fl.fields(24).name='lngthtranspulse';   fl.fields(24).size=1;    fl.fields(24).type='uint16';  fl.fields(24).offset=35;
fl.fields(25).name='watrefbins';        fl.fields(25).size=2;    fl.fields(25).type='uint8';   fl.fields(25).offset=37;
fl.fields(26).name='mintarget';         fl.fields(26).size=1;    fl.fields(26).type='uint8';   fl.fields(26).offset=39;
fl.fields(27).name='lowlattrig';        fl.fields(27).size=1;    fl.fields(27).type='uint8';   fl.fields(27).offset=40;
fl.fields(28).name='distpulse';         fl.fields(28).size=1;    fl.fields(28).type='uint16';  fl.fields(28).offset=41;
fl.fields(29).name='cpuserial';         fl.fields(29).size=8;    fl.fields(29).type='uint8';   fl.fields(29).offset=43;
fl.fields(30).name='bandwidth';         fl.fields(30).size=1;    fl.fields(30).type='uint16';  fl.fields(30).offset=51;
fl.fields(31).name='syspower';          fl.fields(31).size=1;    fl.fields(31).type='uint8';   fl.fields(31).offset=53;
fl.fields(32).name='basefreqid';        fl.fields(32).size=1;    fl.fields(32).type='uint8';   fl.fields(32).offset=54;
fl.fields(33).name='serial';            fl.fields(33).size=4;    fl.fields(33).type='uint8';   fl.fields(33).offset=55;
fl.fields(34).name='beamangle';         fl.fields(34).size=1;    fl.fields(34).type='uint8';   fl.fields(34).offset=59;
datadefs(1)=fl;

%%% VARIABLE LEADER %%%
vl.name='Variable Leader';
vl.header=128;
vl.type='fields';
vl.fields(1 ).name='ensnum_lsb';        vl.fields(1 ).size=1;    vl.fields(1 ).type='uint16';  vl.fields(1 ).offset=3;
vl.fields(2 ).name='time';              vl.fields(2 ).size=7;    vl.fields(2 ).type='uint8';   vl.fields(2 ).offset=5;
vl.fields(3 ).name='ensnum_msb';        vl.fields(3 ).size=1;    vl.fields(3 ).type='uint8';   vl.fields(3 ).offset=12;
vl.fields(4 ).name='BITcheck';          vl.fields(4 ).size=2;    vl.fields(4 ).type='uint8';   vl.fields(4 ).offset=13;
vl.fields(5 ).name='speedsound';        vl.fields(5 ).size=1;    vl.fields(5 ).type='uint16';  vl.fields(5 ).offset=15;
vl.fields(6 ).name='depthtransd';       vl.fields(6 ).size=1;    vl.fields(6 ).type='uint16';  vl.fields(6 ).offset=17;
vl.fields(7 ).name='heading';           vl.fields(7 ).size=1;    vl.fields(7 ).type='uint16';  vl.fields(7 ).offset=19;
vl.fields(8 ).name='pitch';             vl.fields(8 ).size=1;    vl.fields(8 ).type='int16';   vl.fields(8 ).offset=21;
vl.fields(9 ).name='roll';              vl.fields(9 ).size=1;    vl.fields(9 ).type='int16';   vl.fields(9 ).offset=23;
vl.fields(10).name='salinity';          vl.fields(10).size=1;    vl.fields(10).type='uint16';  vl.fields(10).offset=25;
vl.fields(11).name='temperature';       vl.fields(11).size=1;    vl.fields(11).type='int16';   vl.fields(11).offset=27;
vl.fields(12).name='prepingT_min';      vl.fields(12).size=1;    vl.fields(12).type='uint8';   vl.fields(12).offset=29;
vl.fields(13).name='prepingT_sec';      vl.fields(13).size=1;    vl.fields(13).type='uint8';   vl.fields(13).offset=30;
vl.fields(14).name='prepingT_cs';       vl.fields(14).size=1;    vl.fields(14).type='uint8';   vl.fields(14).offset=31;
vl.fields(15).name='headstd';           vl.fields(15).size=1;    vl.fields(15).type='uint8';   vl.fields(15).offset=32;
vl.fields(16).name='pitchstd';          vl.fields(16).size=1;    vl.fields(16).type='uint8';   vl.fields(16).offset=33;
vl.fields(17).name='rollstd';           vl.fields(17).size=1;    vl.fields(17).type='uint8';   vl.fields(17).offset=34;
vl.fields(18).name='ADC';               vl.fields(18).size=8;    vl.fields(18).type='uint8';   vl.fields(18).offset=35;
vl.fields(19).name='errorstat';         vl.fields(19).size=4;    vl.fields(19).type='uint8';   vl.fields(19).offset=43;
vl.fields(20).name='pressure';          vl.fields(20).size=1;    vl.fields(20).type='uint32';  vl.fields(20).offset=49;
vl.fields(21).name='pressurevar';       vl.fields(21).size=1;    vl.fields(21).type='uint32';  vl.fields(21).offset=53;
vl.fields(22).name='time_y2k';          vl.fields(22).size=8;    vl.fields(22).type='uint8';   vl.fields(22).offset=58;
datadefs(2)=vl;

%%% VELOCITY DATA %%%
vel.name='Velocity';
vel.header=256;
vel.type='array';
vel.fields(1).name='VEL';
vel.fields(1).type='int16';
datadefs(3)=vel;

%%% ECHO DATA %%%
echo.name='Echo Intensity';
echo.header=768;
echo.type='array';
echo.fields(1).name='ECHO';
echo.fields(1).type='uint8';
datadefs(4)=echo;

%%% CORRELATION DATA %%%
corr.name='Correlation';
corr.header=512;
corr.type='array';
corr.fields(1).name='CORR';
corr.fields(1).type='uint8';
datadefs(5)=corr;

%%% PERCENTAGE GOOD DATA %%%
perc.name='Percentage good';
perc.header=1024;
perc.type='array';
perc.fields(1).name='PERC';
perc.fields(1).type='uint8';
datadefs(6)=perc;


end
