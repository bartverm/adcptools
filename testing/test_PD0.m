clearvars
clear classes

addpath ../branches/oo/
addpath ../trunk/

%% Construct empty
a=rdi.PD0
assert(isa(a,'rdi.PD0'));

%% Copy construct
b=rdi.PD0(a)
assert(a==b);
clear a b

%% Construct with buffer
fid=fopen('MASTR000.000','r','l');
dat=fread(fid,'*uint8');
fclose(fid);
a=rdi.PD0(dat);
assert(isa(a,'rdi.PD0'));

%% Construct with file
a=rdi.PD0('MASTR000.000');
assert(isa(a,'rdi.PD0'));


%% Compare parsing speed
tic, dat=readADCP('MASTR000.000'); t1=toc;
tic, a=rdi.PD0('MASTR000.000'), t2=toc;

