function datadefs=define_data_types()
% defines PD0 data types

% Remarks:
%   -Make sure the nbins and usedbeams are always in the same sturcture as the array fields that need these variables
%   -Seems the name field is unused now


%%% FIXED LEADER %%% (tested)
fl.name='fixed_leader';
fl.header='0';
fl.type='fields';
fl.in_structure='root';
fl.postaction=[];
fl.fields(1 ).name='firmver';           fl.fields(1 ).size=1;    fl.fields(1 ).type='uint8';   fl.fields(1 ).offset=3;
fl.fields(2 ).name='firmrev';           fl.fields(2 ).size=1;    fl.fields(2 ).type='uint8';   fl.fields(2 ).offset=4;
fl.fields(3 ).name='sysconf';           fl.fields(3 ).size=1;    fl.fields(3 ).type='uint16';  fl.fields(3 ).offset=5;
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

%%% VARIABLE LEADER %%% (tested)
vl.name='variable_leader';
vl.header='0080';
vl.type='fields';
vl.in_structure='root';
vl.postaction=[];
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
datadefs(end+1)=vl;

%%% VELOCITY DATA %%% (tested)
vel.name='velocity';
vel.header='0100';
vel.type='array';
vel.in_structure='root';
vel.postaction=[];
vel.fields(1).name='VEL'; vel.fields(1).type='int16';
datadefs(end+1)=vel;

%%% ECHO DATA %%% (tested)
echo.name='echo';
echo.header='0200';
echo.type='array';
echo.in_structure='root';
echo.postaction=[];
echo.fields(1).name='ECHO'; echo.fields(1).type='uint8';
datadefs(end+1)=echo;

%%% CORRELATION DATA %%% (tested)
corr.name='correlation';
corr.header='0300';
corr.type='array';
corr.in_structure='root';
corr.postaction=[];
corr.fields(1).name='CORR'; corr.fields(1).type='uint8';
datadefs(end+1)=corr;

%%% PERCENTAGE GOOD DATA %%% (tested)
perc.name='percentage';
perc.header='0400';
perc.type='array';
perc.in_structure='root';
perc.postaction=[];
perc.fields(1).name='PERC'; perc.fields(1).type='uint8';
datadefs(end+1)=perc;

%% BOTTOM TRACKING DATA %%% (tested)
bt.name='bottom_tracking';
bt.header='0600';
bt.type='fields';
bt.in_structure='root';
bt.postaction=[];
bt.fields(1 ).name='btpingperens';   bt.fields(1 ).size=1;   bt.fields(1 ).type='uint16';  bt.fields(1 ).offset=3;
bt.fields(2 ).name='reacqdelay';     bt.fields(2 ).size=1;   bt.fields(2 ).type='uint16';  bt.fields(2 ).offset=5;
bt.fields(3 ).name='mincormag';      bt.fields(3 ).size=1;   bt.fields(3 ).type='uint8';   bt.fields(3 ).offset=7;
bt.fields(4 ).name='minevampl';      bt.fields(4 ).size=1;   bt.fields(4 ).type='uint8';   bt.fields(4 ).offset=8;
bt.fields(5 ).name='btminpergood';   bt.fields(5 ).size=1;   bt.fields(5 ).type='uint8';   bt.fields(5 ).offset=9;
bt.fields(6 ).name='btmode';         bt.fields(6 ).size=1;   bt.fields(6 ).type='uint8';   bt.fields(6 ).offset=10;
bt.fields(7 ).name='btmaxerrv';      bt.fields(7 ).size=1;   bt.fields(7 ).type='uint16';  bt.fields(7 ).offset=11;
bt.fields(8 ).name='btrange';        bt.fields(8 ).size=4;   bt.fields(8 ).type='uint16';  bt.fields(8 ).offset=17;
bt.fields(9 ).name='btvel';          bt.fields(9 ).size=4;   bt.fields(9 ).type='int16';   bt.fields(9 ).offset=25;
bt.fields(10).name='btcor';          bt.fields(10).size=4;   bt.fields(10).type='uint8';   bt.fields(10).offset=33;
bt.fields(11).name='btevampl';       bt.fields(11).size=1;   bt.fields(11).type='uint8';   bt.fields(11).offset=37;
bt.fields(12).name='btpercgood';     bt.fields(12).size=1;   bt.fields(12).type='uint8';   bt.fields(12).offset=41;
bt.fields(13).name='minlyrsize';     bt.fields(13).size=1;   bt.fields(13).type='uint16';  bt.fields(13).offset=45;
bt.fields(14).name='nearbnd';        bt.fields(14).size=1;   bt.fields(14).type='uint16';  bt.fields(14).offset=47;
bt.fields(15).name='farbnd';         bt.fields(15).size=1;   bt.fields(15).type='uint16';  bt.fields(15).offset=49;
bt.fields(16).name='reflyrvel';      bt.fields(16).size=4;   bt.fields(16).type='int16';   bt.fields(16).offset=51;
bt.fields(17).name='reflyrcor';      bt.fields(17).size=4;   bt.fields(17).type='uint8';   bt.fields(17).offset=59;
bt.fields(18).name='reflyrint';      bt.fields(18).size=4;   bt.fields(18).type='uint8';   bt.fields(18).offset=63;
bt.fields(19).name='reflyrpergood';  bt.fields(19).size=4;   bt.fields(19).type='uint8';   bt.fields(19).offset=67;
bt.fields(20).name='maxdepth';       bt.fields(20).size=1;   bt.fields(20).type='uint16';  bt.fields(20).offset=71;
bt.fields(21).name='rssiamp';        bt.fields(21).size=4;   bt.fields(21).type='uint8';   bt.fields(21).offset=73;
bt.fields(22).name='gain';           bt.fields(22).size=1;   bt.fields(22).type='uint8';   bt.fields(22).offset=77;
bt.fields(23).name='btrange_msb';    bt.fields(23).size=4;   bt.fields(23).type='uint8';   bt.fields(23).offset=78;
datadefs(end+1)=bt;

%% SURFACE LAYER VELOCITY LEADER %%
sll.name='surface_layer_leader';
sll.header='0010';
sll.type='fields';
sll.in_structure='surface_layer';
sll.postaction=@(x) set_usedbeams_to_n(x,4,size(x.nbins,2));
sll.fields(1).name='nbins';    sll.fields(1).size=1; sll.fields(1).type='uint8';  sll.fields(1).offset=3;
sll.fields(2).name='binsize';  sll.fields(2).size=1; sll.fields(2).type='uint16'; sll.fields(2).offset=4;
sll.fields(3).name='bin1dist'; sll.fields(3).size=1; sll.fields(3).type='uint16'; sll.fields(3).offset=6;
datadefs(end+1)=sll;

%% SURFACE LAYER VELOCITY %%
slv.name='surface_layer_velocity';
slv.header='0110';
slv.type='array';
slv.in_structure='surface_layer';
slv.postaction=[];
slv.fields(1).name='vel'; slv.fields(1).type='int16';
datadefs(end+1)=slv;

%% SURFACE LAYER CORRELATION %%
slc.name='surface_layer_correlation';
slc.header='0210';
slc.type='array';
slc.in_structure='surface_layer';
slc.postaction=[];
slc.fields(1).name='corr'; slc.fields(1).type='uint8';
datadefs(end+1)=slc;

%% SURFACE LAYER ECHO %%
sle.name='surface_layer_echo';
sle.header='0310';
sle.type='array';
sle.in_structure='surface_layer';
sle.postaction=[];
sle.fields(1).name='echo'; sle.fields(1).type='uint8';
datadefs(end+1)=sle;

%% SURFACE LAYER PERCENTAGE GOOD %%
slp.name='surface_layer_percentage';
slp.header='0410';
slp.type='array';
slp.in_structure='surface_layer';
slp.postaction=[];
slp.fields(1).name='perc'; slp.fields(1).type='uint8';
datadefs(end+1)=slp;

%% SURFACE LAYER STATUS %%
sls.name='surface_layer_status';
sls.header='0510';
sls.type='array';
sls.in_structure='surface_layer';
sls.postaction=[];
sls.fields(1).name='status'; sls.fields(1).type='uint8';
datadefs(end+1)=sls;

%% AUTOMATIC MODE 3 %%
am3.name='mode3_settings';
am3.header='4401';
am3.type='fields';
am3.in_structure='mode3_settings';
am3.postaction=[];
am3.fields( 1).name='n_beams';           am3.fields( 1).size=1; am3.fields( 1).type='uint8' ; am3.fields( 1).offset=3;
am3.fields( 2).name='beam1_setup';       am3.fields( 2).size=1; am3.fields( 2).type='uint8' ; am3.fields( 2).offset=4;
am3.fields( 3).name='beam1_depth';       am3.fields( 3).size=1; am3.fields( 3).type='uint16'; am3.fields( 3).offset=5;
am3.fields( 4).name='beam1_ping_count';  am3.fields( 4).size=1; am3.fields( 4).type='uint8' ; am3.fields( 4).offset=7;
am3.fields( 5).name='beam1_ping_type';   am3.fields( 5).size=1; am3.fields( 5).type='uint8' ; am3.fields( 5).offset=8;
am3.fields( 6).name='beam1_nbins';       am3.fields( 6).size=1; am3.fields( 6).type='uint16'; am3.fields( 6).offset=9;
am3.fields( 7).name='beam1_bin_size';    am3.fields( 7).size=1; am3.fields( 7).type='uint16'; am3.fields( 7).offset=11;
am3.fields( 8).name='beam1_bin1dist';    am3.fields( 8).size=1; am3.fields( 8).type='uint16'; am3.fields( 8).offset=13;
am3.fields( 9).name='beam1_code_resps';  am3.fields( 9).size=1; am3.fields( 9).type='uint8' ; am3.fields( 9).offset=15;
am3.fields(10).name='beam1_xmit_length'; am3.fields(10).size=1; am3.fields(10).type='uint16'; am3.fields(10).offset=16;
am3.fields(11).name='beam1_lag_length';  am3.fields(11).size=1; am3.fields(11).type='uint16'; am3.fields(11).offset=18;
am3.fields(12).name='beam1_xmit_bwidth'; am3.fields(12).size=1; am3.fields(12).type='uint8' ; am3.fields(12).offset=20;
am3.fields(13).name='beam1_rcv_bwidth';  am3.fields(13).size=1; am3.fields(13).type='uint8' ; am3.fields(13).offset=21;
am3.fields(14).name='beam1_min_ping_t';  am3.fields(14).size=1; am3.fields(14).type='uint16'; am3.fields(14).offset=22;
am3.fields(15).name='beam2_setup';       am3.fields(15).size=1; am3.fields(15).type='uint8' ; am3.fields(15).offset=24;
am3.fields(16).name='beam2_depth';       am3.fields(16).size=1; am3.fields(16).type='uint16'; am3.fields(16).offset=25;
am3.fields(17).name='beam2_ping_count';  am3.fields(17).size=1; am3.fields(17).type='uint8' ; am3.fields(17).offset=27;
am3.fields(18).name='beam2_ping_type';   am3.fields(18).size=1; am3.fields(18).type='uint8' ; am3.fields(18).offset=28;
am3.fields(19).name='beam2_nbins';       am3.fields(19).size=1; am3.fields(19).type='uint16'; am3.fields(19).offset=29;
am3.fields(20).name='beam2_bin_size';    am3.fields(20).size=1; am3.fields(20).type='uint16'; am3.fields(20).offset=31;
am3.fields(21).name='beam2_bin1dist';    am3.fields(21).size=1; am3.fields(21).type='uint16'; am3.fields(21).offset=33;
am3.fields(22).name='beam2_code_resps';  am3.fields(22).size=1; am3.fields(22).type='uint8' ; am3.fields(22).offset=35;
am3.fields(23).name='beam2_xmit_length'; am3.fields(23).size=1; am3.fields(23).type='uint16'; am3.fields(23).offset=36;
am3.fields(24).name='beam2_lag_length';  am3.fields(24).size=1; am3.fields(24).type='uint16'; am3.fields(24).offset=38;
am3.fields(25).name='beam2_xmit_bwidth'; am3.fields(25).size=1; am3.fields(25).type='uint8' ; am3.fields(25).offset=40;
am3.fields(26).name='beam2_rcv_bwidth';  am3.fields(26).size=1; am3.fields(26).type='uint8' ; am3.fields(26).offset=41;
am3.fields(27).name='beam2_min_ping_t';  am3.fields(27).size=1; am3.fields(27).type='uint16'; am3.fields(27).offset=42;
am3.fields(28).name='beam3_setup';       am3.fields(28).size=1; am3.fields(28).type='uint8' ; am3.fields(28).offset=44;
am3.fields(29).name='beam3_depth';       am3.fields(29).size=1; am3.fields(29).type='uint16'; am3.fields(29).offset=45;
am3.fields(30).name='beam3_ping_count';  am3.fields(30).size=1; am3.fields(30).type='uint8' ; am3.fields(30).offset=47;
am3.fields(31).name='beam3_ping_type';   am3.fields(31).size=1; am3.fields(31).type='uint8' ; am3.fields(31).offset=48;
am3.fields(32).name='beam3_nbins';       am3.fields(32).size=1; am3.fields(32).type='uint16'; am3.fields(32).offset=49;
am3.fields(33).name='beam3_bin_size';    am3.fields(33).size=1; am3.fields(33).type='uint16'; am3.fields(33).offset=51;
am3.fields(34).name='beam3_bin1dist';    am3.fields(34).size=1; am3.fields(34).type='uint16'; am3.fields(34).offset=53;
am3.fields(35).name='beam3_code_resps';  am3.fields(35).size=1; am3.fields(35).type='uint8' ; am3.fields(35).offset=55;
am3.fields(36).name='beam3_xmit_length'; am3.fields(36).size=1; am3.fields(36).type='uint16'; am3.fields(36).offset=56;
am3.fields(37).name='beam3_lag_length';  am3.fields(37).size=1; am3.fields(37).type='uint16'; am3.fields(37).offset=58;
am3.fields(38).name='beam3_xmit_bwidth'; am3.fields(38).size=1; am3.fields(38).type='uint8' ; am3.fields(38).offset=60;
am3.fields(39).name='beam3_rcv_bwidth';  am3.fields(39).size=1; am3.fields(39).type='uint8' ; am3.fields(39).offset=61;
am3.fields(40).name='beam3_min_ping_t';  am3.fields(40).size=1; am3.fields(40).type='uint16'; am3.fields(40).offset=62;
am3.fields(41).name='beam4_setup';       am3.fields(41).size=1; am3.fields(41).type='uint8' ; am3.fields(41).offset=64;
am3.fields(42).name='beam4_depth';       am3.fields(42).size=1; am3.fields(42).type='uint16'; am3.fields(42).offset=65;
am3.fields(43).name='beam4_ping_count';  am3.fields(43).size=1; am3.fields(43).type='uint8' ; am3.fields(43).offset=67;
am3.fields(44).name='beam4_ping_type';   am3.fields(44).size=1; am3.fields(44).type='uint8' ; am3.fields(44).offset=68;
am3.fields(45).name='beam4_nbins';       am3.fields(45).size=1; am3.fields(45).type='uint16'; am3.fields(45).offset=69;
am3.fields(46).name='beam4_bin_size';    am3.fields(46).size=1; am3.fields(46).type='uint16'; am3.fields(46).offset=71;
am3.fields(47).name='beam4_bin1dist';    am3.fields(47).size=1; am3.fields(47).type='uint16'; am3.fields(47).offset=73;
am3.fields(48).name='beam4_code_resps';  am3.fields(48).size=1; am3.fields(48).type='uint8' ; am3.fields(48).offset=75;
am3.fields(49).name='beam4_xmit_length'; am3.fields(49).size=1; am3.fields(49).type='uint16'; am3.fields(49).offset=76;
am3.fields(50).name='beam4_lag_length';  am3.fields(50).size=1; am3.fields(50).type='uint16'; am3.fields(50).offset=78;
am3.fields(51).name='beam4_xmit_bwidth'; am3.fields(51).size=1; am3.fields(51).type='uint8' ; am3.fields(51).offset=80;
am3.fields(52).name='beam4_rcv_bwidth';  am3.fields(52).size=1; am3.fields(52).type='uint8' ; am3.fields(52).offset=81;
am3.fields(53).name='beam4_min_ping_t';  am3.fields(53).size=1; am3.fields(53).type='uint16'; am3.fields(53).offset=82;
am3.fields(54).name='reserved';          am3.fields(54).size=1; am3.fields(54).type='uint8' ; am3.fields(54).offset=84;
datadefs(end+1)=am3;

%% FIRMWARE STATUS %%
fs.name='firmware_status';
fs.header='4400';
fs.type='fields';
fs.in_structure='firmware_status';
fs.postaction=[];
fs.fields(1).name='version_alpha';  fs.fields(1).size=1;  fs.fields(1).type='char';  fs.fields(1).offset=3;
fs.fields(2).name='version_branch'; fs.fields(2).size=14; fs.fields(2).type='char';  fs.fields(2).offset=4;
fs.fields(3).name='test_data';      fs.fields(3).size=2;  fs.fields(3).type='char';  fs.fields(3).offset=18;
fs.fields(4).name='test_switches';  fs.fields(4).size=2;  fs.fields(4).type='char';  fs.fields(4).offset=20;
fs.fields(5).name='reserved';       fs.fields(5).size=1;  fs.fields(5).type='uint8'; fs.fields(5).offset=22;
datadefs(end+1)=fs;

%% VERTICAL BEAM RANGE %%
vr.name='vertical_beam_range';
vr.header='4100';
vr.type='fields';
vr.in_structure='vertical_beam';
vr.postaction=[];
vr.fields(1).name='eval_amp'; vr.fields(1).size=1; vr.fields(1).type='uint8' ; vr.fields(1).offset=3;
vr.fields(2).name='rssi_amp'; vr.fields(2).size=1; vr.fields(2).type='uint8' ; vr.fields(2).offset=4;
vr.fields(3).name='range';    vr.fields(3).size=1; vr.fields(3).type='uint32'; vr.fields(3).offset=5;
vr.fields(4).name='status';   vr.fields(4).size=1; vr.fields(4).type='uint8' ; vr.fields(4).offset=9;
datadefs(end+1)=vr;

%% VERTICAL BEAM LEADER %%
vbl.name='vertical_beam_leader';
vbl.header='0F01';
vbl.type='fields';
vbl.in_structure='vertical_beam';
vbl.postaction=@(x) set_usedbeams_to_n(x,1,size(nbins,2));
vbl.fields(1).name='nbins';       vbl.fields(1).size=1; vbl.fields(1).type='uint16'; vbl.fields(1).offset=3;
vbl.fields(2).name='npings';      vbl.fields(2).size=1; vbl.fields(2).type='uint16'; vbl.fields(2).offset=5;
vbl.fields(3).name='binsize';     vbl.fields(3).size=1; vbl.fields(3).type='uint16'; vbl.fields(3).offset=7;
vbl.fields(4).name='bin1dist';    vbl.fields(4).size=1; vbl.fields(4).type='uint16'; vbl.fields(4).offset=9;
vbl.fields(5).name='xmit_length'; vbl.fields(5).size=1; vbl.fields(5).type='uint16'; vbl.fields(5).offset=13;
vbl.fields(6).name='lag_length';  vbl.fields(6).size=1; vbl.fields(6).type='uint16'; vbl.fields(6).offset=15;
vbl.fields(7).name='ncod_el';     vbl.fields(7).size=1; vbl.fields(7).type='uint16'; vbl.fields(7).offset=17;
datadefs(end+1)=vbl;

%% VERTICAL BEAM VELOCITY %%
vbv.name='vertical_beam_velocity';
vbv.header='0A00';
vbv.type='array';
vbv.in_structure='vertical_beam';
vbv.postaction=[];
vbv.fields(1).name='vel'; vbv.fields(1).type='int16';
datadefs(end+1)=vbv;

%% VERTICAL BEAM CORRELATION %%
vbc.name='vertical_beam_correlation';
vbc.header='0B00';
vbc.type='array';
vbc.in_structure='vertical_beam';
vbc.postaction=[];
vbc.fields(1).name='corr'; vbc.fields(1).type='uint8';
datadefs(end+1)=vbc;

%% VERTICAL BEAM ECHO %%
vbe.name='vertical_beam_echo';
vbe.header='0C00';
vbe.type='array';
vbe.in_structure='vertical_beam';
vbe.postaction=[];
vbe.fields(1).name='echo'; vbe.fields(1).type='uint8';
datadefs(end+1)=vbe;

%% VERTICAL BEAM PERCENTAGE GOOD %%
vbp.name='vertical_beam_percentage';
vbp.header='0D00';
vbp.type='array';
vbp.in_structure='vertical_beam';
vbp.postaction=[];
vbp.fields(1).name='perc'; vbp.fields(1).type='uint8';
datadefs(end+1)=vbp;

%% VERTICAL BEAM STATUS %%
vbs.name='vertical_beam_status';
vbs.header='0E00';
vbs.type='array';
vbs.in_structure='vertical_beam';
vbs.postaction=[];
vbs.fields(1).name='status'; vbs.fields(1).type='uint8';
datadefs(end+1)=vbs;

%% BEAM CORRECTION MATRIX %%
bcm.name='beam_correction_matrix';
bcm.header='3200';
bcm.type='fields';
bcm.in_structure='beam_correction_matrix';
bcm.postaction=@reshape_corr_matrix;
bcm.fields(1).name='cal_matrix'; bcm.fields(1).size=16; bcm.fields(1).type='int16'; bcm.fields(1).offset=3;
datadefs(end+1)=bcm;


%% NMEA WRII NEW GGA %% (tested)

gga_wr2_2.name='gga_wr2_2';
gga_wr2_2.header='0068';
gga_wr2_2.type='fields';
gga_wr2_2.in_structure='gga_wr2_2';
gga_wr2_2.postaction=[];
gga_wr2_2.fields(1 ).name='dt';             gga_wr2_2.fields(1 ).size=1;   gga_wr2_2.fields(1 ).type='double';   gga_wr2_2.fields(1 ).offset=7 ;
gga_wr2_2.fields(2 ).name='header';         gga_wr2_2.fields(2 ).size=7;   gga_wr2_2.fields(2 ).type='char';     gga_wr2_2.fields(2 ).offset=15;
gga_wr2_2.fields(3 ).name='utc';            gga_wr2_2.fields(3 ).size=10;  gga_wr2_2.fields(3 ).type='char';     gga_wr2_2.fields(3 ).offset=22;
gga_wr2_2.fields(4 ).name='latitude';       gga_wr2_2.fields(4 ).size=1;   gga_wr2_2.fields(4 ).type='double';   gga_wr2_2.fields(4 ).offset=32;
gga_wr2_2.fields(5 ).name='lat_ns';         gga_wr2_2.fields(5 ).size=1;   gga_wr2_2.fields(5 ).type='char';     gga_wr2_2.fields(5 ).offset=40;
gga_wr2_2.fields(6 ).name='longitude';      gga_wr2_2.fields(6 ).size=1;   gga_wr2_2.fields(6 ).type='double';   gga_wr2_2.fields(6 ).offset=41;
gga_wr2_2.fields(7 ).name='lon_ew';         gga_wr2_2.fields(7 ).size=1;   gga_wr2_2.fields(7 ).type='char';     gga_wr2_2.fields(7 ).offset=49;
gga_wr2_2.fields(8 ).name='quality';        gga_wr2_2.fields(8 ).size=1;   gga_wr2_2.fields(8 ).type='uint8';    gga_wr2_2.fields(8 ).offset=50;
gga_wr2_2.fields(9 ).name='num_sat';        gga_wr2_2.fields(9 ).size=1;   gga_wr2_2.fields(9 ).type='uint8';    gga_wr2_2.fields(9 ).offset=51;
gga_wr2_2.fields(10).name='hdop';           gga_wr2_2.fields(10).size=1;   gga_wr2_2.fields(10).type='single';   gga_wr2_2.fields(10).offset=52;
gga_wr2_2.fields(11).name='alt';            gga_wr2_2.fields(11).size=1;   gga_wr2_2.fields(11).type='single';   gga_wr2_2.fields(11).offset=56;
gga_wr2_2.fields(12).name='alt_unit';       gga_wr2_2.fields(12).size=1;   gga_wr2_2.fields(12).type='char';     gga_wr2_2.fields(12).offset=60;
gga_wr2_2.fields(13).name='geoid';          gga_wr2_2.fields(13).size=1;   gga_wr2_2.fields(13).type='single';   gga_wr2_2.fields(13).offset=61;
gga_wr2_2.fields(14).name='geoid_unit';     gga_wr2_2.fields(14).size=1;   gga_wr2_2.fields(14).type='char';     gga_wr2_2.fields(14).offset=65;
gga_wr2_2.fields(15).name='age_dgps';       gga_wr2_2.fields(15).size=1;   gga_wr2_2.fields(15).type='single';   gga_wr2_2.fields(15).offset=66;
gga_wr2_2.fields(16).name='ref_station_id'; gga_wr2_2.fields(16).size=1;   gga_wr2_2.fields(16).type='uint16';   gga_wr2_2.fields(16).offset=70;
datadefs(end+1)=gga_wr2_2;

%% NMEA WRII NEW VTG %% (tested)
vtg_wr2_2.name='vtg_wr2_2';
vtg_wr2_2.header='0069';
vtg_wr2_2.type='fields';
vtg_wr2_2.in_structure='vtg_wr2_2';
vtg_wr2_2.postaction=[];
vtg_wr2_2.fields(1 ).name='dt';                    vtg_wr2_2.fields(1 ).size=1; vtg_wr2_2.fields(1 ).type='double'; vtg_wr2_2.fields(1 ).offset=7 ;
vtg_wr2_2.fields(2 ).name='header';                vtg_wr2_2.fields(2 ).size=7; vtg_wr2_2.fields(2 ).type='char';   vtg_wr2_2.fields(2 ).offset=15;
vtg_wr2_2.fields(3 ).name='track_dir_true';        vtg_wr2_2.fields(3 ).size=1; vtg_wr2_2.fields(3 ).type='single'; vtg_wr2_2.fields(3 ).offset=22;
vtg_wr2_2.fields(4 ).name='true_indicator';        vtg_wr2_2.fields(4 ).size=1; vtg_wr2_2.fields(4 ).type='char';   vtg_wr2_2.fields(4 ).offset=26;
vtg_wr2_2.fields(5 ).name='track_dir_magn';        vtg_wr2_2.fields(5 ).size=1; vtg_wr2_2.fields(5 ).type='single'; vtg_wr2_2.fields(5 ).offset=27;
vtg_wr2_2.fields(6 ).name='magn_indicator';        vtg_wr2_2.fields(6 ).size=1; vtg_wr2_2.fields(6 ).type='char';   vtg_wr2_2.fields(6 ).offset=31;
vtg_wr2_2.fields(7 ).name='speed_over_ground_kts'; vtg_wr2_2.fields(7 ).size=1; vtg_wr2_2.fields(7 ).type='single'; vtg_wr2_2.fields(7 ).offset=32;
vtg_wr2_2.fields(8 ).name='kts_indicator';         vtg_wr2_2.fields(8 ).size=1; vtg_wr2_2.fields(8 ).type='char';   vtg_wr2_2.fields(8 ).offset=36;
vtg_wr2_2.fields(9 ).name='speed_over_ground_kmh'; vtg_wr2_2.fields(9 ).size=1; vtg_wr2_2.fields(9 ).type='single'; vtg_wr2_2.fields(9 ).offset=37;
vtg_wr2_2.fields(10).name='kmh_indicator';         vtg_wr2_2.fields(10).size=1; vtg_wr2_2.fields(10).type='char';   vtg_wr2_2.fields(10).offset=41;
vtg_wr2_2.fields(11).name='mode_indicator';        vtg_wr2_2.fields(11).size=1; vtg_wr2_2.fields(11).type='char';   vtg_wr2_2.fields(11).offset=42;
datadefs(end+1)=vtg_wr2_2;

%% NMEA WRII NEW HDT %% (untested)
hdt_wr2_2.name='hdt_wr2_2';
hdt_wr2_2.header='006B';
hdt_wr2_2.type='fields';
hdt_wr2_2.in_structure='hdt_wr2_2';
hdt_wr2_2.postaction=[];
hdt_wr2_2.fields(1 ).name='dt';                    hdt_wr2_2.fields(1 ).size=1; hdt_wr2_2.fields(1 ).type='double'; hdt_wr2_2.fields(1 ).offset=7 ;
hdt_wr2_2.fields(2 ).name='header';                hdt_wr2_2.fields(2 ).size=7; hdt_wr2_2.fields(2 ).type='char';   hdt_wr2_2.fields(2 ).offset=15;
hdt_wr2_2.fields(3 ).name='heading';               hdt_wr2_2.fields(3 ).size=1; hdt_wr2_2.fields(3 ).type='double'; hdt_wr2_2.fields(3 ).offset=22;
hdt_wr2_2.fields(4 ).name='true_indicator';        hdt_wr2_2.fields(4 ).size=1; hdt_wr2_2.fields(4 ).type='char';   hdt_wr2_2.fields(4 ).offset=30;
datadefs(end+1)=hdt_wr2_2;

%% NMEA WRII NEW DBT %% (untested)
dbt_wr2_2.name='dbt_wr2_2';
dbt_wr2_2.header='006A';
dbt_wr2_2.type='fields';
dbt_wr2_2.in_structure='dbt_wr2_2';
dbt_wr2_2.postaction=[];
dbt_wr2_2.fields(1 ).name='dt';            dbt_wr2_2.fields(1 ).size=1; dbt_wr2_2.fields(1 ).type='double'; dbt_wr2_2.fields(1 ).offset=7 ;
dbt_wr2_2.fields(2 ).name='header';        dbt_wr2_2.fields(2 ).size=7; dbt_wr2_2.fields(2 ).type='char';   dbt_wr2_2.fields(2 ).offset=15;
dbt_wr2_2.fields(3 ).name='water_depth_ft';dbt_wr2_2.fields(3 ).size=1; dbt_wr2_2.fields(3 ).type='single'; dbt_wr2_2.fields(3 ).offset=22;
dbt_wr2_2.fields(4 ).name='ft_indicator';  dbt_wr2_2.fields(4 ).size=1; dbt_wr2_2.fields(4 ).type='char';   dbt_wr2_2.fields(4 ).offset=26;
dbt_wr2_2.fields(5 ).name='water_depth_m'; dbt_wr2_2.fields(5 ).size=1; dbt_wr2_2.fields(5 ).type='single'; dbt_wr2_2.fields(5 ).offset=27;
dbt_wr2_2.fields(6 ).name='m_indicator';   dbt_wr2_2.fields(6 ).size=1; dbt_wr2_2.fields(6 ).type='char';   dbt_wr2_2.fields(6 ).offset=31;
dbt_wr2_2.fields(7 ).name='water_depth_f'; dbt_wr2_2.fields(7 ).size=1; dbt_wr2_2.fields(7 ).type='single'; dbt_wr2_2.fields(7 ).offset=32;
dbt_wr2_2.fields(8 ).name='f_indicator';   dbt_wr2_2.fields(8 ).size=1; dbt_wr2_2.fields(8 ).type='char';   dbt_wr2_2.fields(8 ).offset=36;
datadefs(end+1)=dbt_wr2_2;

%% NMEA WRII OLD GGA %% (untested)

gga_wr2.name='gga_wr2';
gga_wr2.header='0064';
gga_wr2.type='fields';
gga_wr2.in_structure='gga_wr2';
gga_wr2.postaction=[];
gga_wr2.fields(1 ).name='dt';             gga_wr2.fields(1 ).size=1;   gga_wr2.fields(1 ).type='double';   gga_wr2.fields(1 ).offset=7 ;
gga_wr2.fields(2 ).name='header';         gga_wr2.fields(2 ).size=10;  gga_wr2.fields(2 ).type='char';     gga_wr2.fields(2 ).offset=15;
gga_wr2.fields(3 ).name='utc';            gga_wr2.fields(3 ).size=10;  gga_wr2.fields(3 ).type='char';     gga_wr2.fields(3 ).offset=25;
gga_wr2.fields(4 ).name='latitude';       gga_wr2.fields(4 ).size=1;   gga_wr2.fields(4 ).type='double';   gga_wr2.fields(4 ).offset=35;
gga_wr2.fields(5 ).name='lat_ns';         gga_wr2.fields(5 ).size=1;   gga_wr2.fields(5 ).type='char';     gga_wr2.fields(5 ).offset=43;
gga_wr2.fields(6 ).name='longitude';      gga_wr2.fields(6 ).size=1;   gga_wr2.fields(6 ).type='double';   gga_wr2.fields(6 ).offset=44;
gga_wr2.fields(7 ).name='lon_ew';         gga_wr2.fields(7 ).size=1;   gga_wr2.fields(7 ).type='char';     gga_wr2.fields(7 ).offset=52;
gga_wr2.fields(8 ).name='quality';        gga_wr2.fields(8 ).size=1;   gga_wr2.fields(8 ).type='uint8';    gga_wr2.fields(8 ).offset=53;
gga_wr2.fields(9 ).name='num_sat';        gga_wr2.fields(9 ).size=1;   gga_wr2.fields(9 ).type='uint8';    gga_wr2.fields(9 ).offset=54;
gga_wr2.fields(10).name='hdop';           gga_wr2.fields(10).size=1;   gga_wr2.fields(10).type='single';   gga_wr2.fields(10).offset=55;
gga_wr2.fields(11).name='alt';            gga_wr2.fields(11).size=1;   gga_wr2.fields(11).type='single';   gga_wr2.fields(11).offset=59;
gga_wr2.fields(12).name='alt_unit';       gga_wr2.fields(12).size=1;   gga_wr2.fields(12).type='char';     gga_wr2.fields(12).offset=63;
gga_wr2.fields(13).name='geoid';          gga_wr2.fields(13).size=1;   gga_wr2.fields(13).type='single';   gga_wr2.fields(13).offset=64;
gga_wr2.fields(14).name='geoid_unit';     gga_wr2.fields(14).size=1;   gga_wr2.fields(14).type='char';     gga_wr2.fields(14).offset=68;
gga_wr2.fields(15).name='age_dgps';       gga_wr2.fields(15).size=1;   gga_wr2.fields(15).type='single';   gga_wr2.fields(15).offset=69;
gga_wr2.fields(16).name='ref_station_id'; gga_wr2.fields(16).size=1;   gga_wr2.fields(16).type='uint16';   gga_wr2.fields(16).offset=73;
datadefs(end+1)=gga_wr2;

%% NMEA WRII OLD VTG %% (untested)
vtg_wr2.name='vtg_wr2';
vtg_wr2.header='0065';
vtg_wr2.type='fields';
vtg_wr2.in_structure='vtg_wr2';
vtg_wr2.postaction=[];
vtg_wr2.fields(1 ).name='dt';                    vtg_wr2.fields(1 ).size=1; vtg_wr2.fields(1 ).type='double'; vtg_wr2.fields(1 ).offset=7 ;
vtg_wr2.fields(2 ).name='header';                vtg_wr2.fields(2 ).size=10;vtg_wr2.fields(2 ).type='char';   vtg_wr2.fields(2 ).offset=15;
vtg_wr2.fields(3 ).name='track_dir_true';        vtg_wr2.fields(3 ).size=1; vtg_wr2.fields(3 ).type='single'; vtg_wr2.fields(3 ).offset=25;
vtg_wr2.fields(4 ).name='true_indicator';        vtg_wr2.fields(4 ).size=1; vtg_wr2.fields(4 ).type='char';   vtg_wr2.fields(4 ).offset=29;
vtg_wr2.fields(5 ).name='track_dir_magn';        vtg_wr2.fields(5 ).size=1; vtg_wr2.fields(5 ).type='single'; vtg_wr2.fields(5 ).offset=30;
vtg_wr2.fields(6 ).name='magn_indicator';        vtg_wr2.fields(6 ).size=1; vtg_wr2.fields(6 ).type='char';   vtg_wr2.fields(6 ).offset=34;
vtg_wr2.fields(7 ).name='speed_over_ground_kts'; vtg_wr2.fields(7 ).size=1; vtg_wr2.fields(7 ).type='single'; vtg_wr2.fields(7 ).offset=35;
vtg_wr2.fields(8 ).name='kts_indicator';         vtg_wr2.fields(8 ).size=1; vtg_wr2.fields(8 ).type='char';   vtg_wr2.fields(8 ).offset=39;
vtg_wr2.fields(9 ).name='speed_over_ground_kmh'; vtg_wr2.fields(9 ).size=1; vtg_wr2.fields(9 ).type='single'; vtg_wr2.fields(9 ).offset=40;
vtg_wr2.fields(10).name='kmh_indicator';         vtg_wr2.fields(10).size=1; vtg_wr2.fields(10).type='char';   vtg_wr2.fields(10).offset=44;
vtg_wr2.fields(11).name='mode_indicator';        vtg_wr2.fields(11).size=1; vtg_wr2.fields(11).type='char';   vtg_wr2.fields(11).offset=45;
datadefs(end+1)=vtg_wr2;

%% NMEA WRII OLD HDT %% (untested)
hdt_wr2.name='hdt_wr2';
hdt_wr2.header='0067';
hdt_wr2.type='fields';
hdt_wr2.in_structure='hdt_wr2';
hdt_wr2.postaction=[];
hdt_wr2.fields(1 ).name='dt';                    hdt_wr2.fields(1 ).size=1; hdt_wr2.fields(1 ).type='double'; hdt_wr2.fields(1 ).offset=7 ;
hdt_wr2.fields(2 ).name='header';                hdt_wr2.fields(2 ).size=10;hdt_wr2.fields(2 ).type='char';   hdt_wr2.fields(2 ).offset=15;
hdt_wr2.fields(3 ).name='heading';               hdt_wr2.fields(3 ).size=1; hdt_wr2.fields(3 ).type='double'; hdt_wr2.fields(3 ).offset=25;
hdt_wr2.fields(4 ).name='true_indicator';        hdt_wr2.fields(4 ).size=1; hdt_wr2.fields(4 ).type='char';   hdt_wr2.fields(4 ).offset=33;
datadefs(end+1)=hdt_wr2;

%% NMEA WRII OLD DBT %% (untested)
dbt_wr2.name='dbt_wr2';
dbt_wr2.header='0066';
dbt_wr2.type='fields';
dbt_wr2.in_structure='dbt_wr2';
dbt_wr2.postaction=[];
dbt_wr2.fields(1 ).name='dt';            dbt_wr2.fields(1 ).size=1; dbt_wr2.fields(1 ).type='double'; dbt_wr2.fields(1 ).offset=7 ;
dbt_wr2.fields(2 ).name='header';        dbt_wr2.fields(2 ).size=10;dbt_wr2.fields(2 ).type='char';   dbt_wr2.fields(2 ).offset=15;
dbt_wr2.fields(3 ).name='water_depth_ft';dbt_wr2.fields(3 ).size=1; dbt_wr2.fields(3 ).type='single'; dbt_wr2.fields(3 ).offset=25;
dbt_wr2.fields(4 ).name='ft_indicator';  dbt_wr2.fields(4 ).size=1; dbt_wr2.fields(4 ).type='char';   dbt_wr2.fields(4 ).offset=29;
dbt_wr2.fields(5 ).name='water_depth_m'; dbt_wr2.fields(5 ).size=1; dbt_wr2.fields(5 ).type='single'; dbt_wr2.fields(5 ).offset=30;
dbt_wr2.fields(6 ).name='m_indicator';   dbt_wr2.fields(6 ).size=1; dbt_wr2.fields(6 ).type='char';   dbt_wr2.fields(6 ).offset=34;
dbt_wr2.fields(7 ).name='water_depth_f'; dbt_wr2.fields(7 ).size=1; dbt_wr2.fields(7 ).type='single'; dbt_wr2.fields(7 ).offset=35;
dbt_wr2.fields(8 ).name='f_indicator';   dbt_wr2.fields(8 ).size=1; dbt_wr2.fields(8 ).type='char';   dbt_wr2.fields(8 ).offset=39;
datadefs(end+1)=dbt_wr2;

%% NMEA WR GGA %%
gga_wr.name='gga_wr';
gga_wr.header='2101';
gga_wr.type='fields';
gga_wr.in_structure='gga_wr';
gga_wr.postaction=@(x) parseNMEA(x.message');
gga_wr.fields(1).name='message';        gga_wr.fields(1).size=94; gga_wr.fields(1).type='char'; gga_wr.fields(1).offset=3;
datadefs(end+1)=gga_wr;

%% NMEA WR DBT %%
dbt_wr.name='dbt_wr';
dbt_wr.header='2100';
dbt_wr.type='fields';
dbt_wr.in_structure='dbt_wr';
dbt_wr.postaction=@(x) parseNMEA(x.message');
dbt_wr.fields(1).name='message';        dbt_wr.fields(1).size=38; dbt_wr.fields(1).type='char'; dbt_wr.fields(1).offset=3;
datadefs(end+1)=dbt_wr;


end


function x=set_usedbeams_to_n(x,n,nens)
  x.usedbeams=uint8(ones(1,nens)*n);
end

function x=reshape_corr_matrix(x)
  x.cal_matrix=permute(reshape(x.cal_matrix,4,4,[]),[2 1 3]);
end

