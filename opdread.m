% opdread - Read a bio-rad OPD file (from iCycler)
% If a filename is given, use it in local directory if it exists
% If filename is empty or not provided, check for a file in current directory
% If there isn't any file in current directory then use newest in server directory
% If 'copyfromserver' is set, and the file in server directory is newer or the local file doesn't exist, then
%  copy from there
function x=opdread(file,varargin)
defaults=struct('copyfromserver',true,'remotedir',[]);
args=processargs(defaults,varargin);
if nargin<1 || isempty(file)
  file='Data *';
end

% Check for local data file
localfile=dir(file);
if length(localfile)<1
  fprintf('Data file "%s" not found in local directory, checking server...\n',file);
elseif length(localfile)>1
  error('Multiple data files in local directory matching "%s"\n',file);
else
  file=localfile.name;
  fprintf('Found local file: %s\n', file);
end

if exist(args.remotedir,'dir') && (isempty(localfile) || args.copyfromserver)
  remotefilename=[args.remotedir,'/',file];
  remotefile=dir(remotefilename);
  if length(remotefile)<1
    fprintf('No remote files found matching %s\n',remotefilename);
  else
    [~,ord]=sort(datenum({remotefile.date}));
    remotefile=remotefile(ord(end));
    fprintf('Remote file: %s\n',remotefilename);
    if ~isempty(localfile) && remotefile.bytes==localfile.bytes
      fprintf('Remote file has same size as localfile, not copying\n');
    elseif ~isempty(localfile) && remotefile.datenum<localfile.datenum
      fprintf('Local file is newer, not copying\n');
    else
      if isempty(localfile)
        file=[args.remotedir,'/',remotefile.name];
      end
      if args.copyfromserver
        cmd=sprintf('cp -p "%s" .',remotefilename);
        fprintf('Executing <%s> ...',cmd);
        [s,r]=system(cmd);
        if s~=0
          error(r);
        else
          fprintf('done\n');
        end
      end
    end
  end
end

if isempty(localfile) && ~exist(args.remotedir)
  if isempty(args.remotedir)
    error('Local file "%s" does not exist and remotedir not specified.\n',file);
  else
    error('Local file "%s" does not exist and remote directory (%s) not mounted.\n',file, args.remotedir);
  end
end

[fd,msg]=fopen(file,'r','ieee-le');
if fd==-1
  error('Unable to open %s: %s',file,msg);
end
dfile=dir(file);
x=struct('filename',dfile.name);
gotWSRT=0;
while 1
  [code,cnt]=fread(fd,4,'uint8=>char');
  if cnt<4
    break;
  end
  code=code';
  if strcmp(code,'WSRT')
    % Called by readWIRT() 
    fread(fd,44,'uint8');
    gotWSRT=1;
    continue;
  elseif gotWSRT
    % Spectral data - hacked to find this block immediately after
    % last WSRT block
    fseek(fd,-4,'cof');  % Back up to reread the code
    x.data=readSpectral(fd,x.TCOD.nsamples);
    gotWSRT=0;
  elseif exist(['read',code])==0
    error(['Invalid section code: 0x%02x%02x%02x%02x at offset ' ...
           '0x%x'],code, ftell(fd)-4);
  elseif isfield(x,code)
    x.(code)(end+1)=eval(['read',code,'(fd,x)']);
  else
    x.(code)=eval(['read',code,'(fd,x)']);
  end
end
fclose(fd);
% Restructure data
x=opdanalyze(x);

% Read n filler bytes from fd and note if they are non-zero
function x=readfill(fd,n)
pad=fread(fd,n,'*uint8');
for i=1:n
  if pad(i)~=0
    error('opdread: filler byte(s) at offset 0x%x is non-zero: 0x%x\n', ...
            ftell(fd)-n+i-1,pad(i));
    return;
  end
end

function x=readString(fd,n)
x=deblank(fread(fd,n,'uint8=>char')');

function x=readTCOD(fd,mx)
x.version=readString(fd,5);
x.user=readString(fd,25);
x.notes=readString(fd,800);
x.filterSet=readFilterSet(fd);
x.numdyes=fread(fd,1,'*uint');
readfill(fd,2);
x.nsamples=fread(fd,1,'*ushort');
readfill(fd,6);
x.protocol=readString(fd,260);
x.plateLayout=readString(fd,260);

function x=readFilterSet(fd)
x.id=fread(fd,1,'*uint');
x.name=readString(fd,20);
for i=1:6
  desc=readString(fd,20);
  index=fread(fd,1,'*ushort');
  x.excitation(i)=struct('desc',desc,'index',index);
end
for i=1:6
  desc=readString(fd,20);
  index=fread(fd,1,'*ushort');
  x.emission(i)=struct('desc',desc,'index',index);
end
readfill(fd,2);

function x=readTCPS(fd,mx)
x.version1=readString(fd,5);
readfill(fd,2);
x.v1=fread(fd,1,'*ushort');
x.user1=readString(fd,825);
x.v2=fread(fd,1,'*ushort');
readfill(fd,6);
x.numdyes=fread(fd,1,'*ushort');
readfill(fd,6);
x.count=fread(fd,1,'*ushort');
readfill(fd,6);
x.version2=readString(fd,5);
readfill(fd,1000);
x.v4=fread(fd,2,'*ushort');
x.v5=fread(fd,3,'*float');
readfill(fd,5);
timestamp=fread(fd,1,'uint');
x.timestamp=timestamp/3600/24+datenum(1970,1,1);  % In UTC
readfill(fd,14);
x.v9=fread(fd,1,'*ushort');
readfill(fd,2);
x.volume=fread(fd,1,'*float');
readfill(fd,8);
x.version3=readString(fd,5);
x.v10=fread(fd,1,'*ushort');
readfill(fd,8);
x.version4=readString(fd,8);
readfill(fd,2);
x.user2=readString(fd,825);
x.v11=fread(fd,1,'*ushort');
readfill(fd,4);
x.nsteps=fread(fd,1,'*ushort');
readfill(fd,30);
x.version5=readString(fd,5);
x.user3=readString(fd,25);
readfill(fd,1105);
for i=1:x.numdyes
  x.layerdata(:,i)=fread(fd,3,'*uint');
end

function x=readWIRT(fd,mx)
readfill(fd,4);
x.platepos=fread(fd,1,'*ushort');
x.numlayers=fread(fd,1,'*ushort');
x.wsrt_offset=fread(fd,1,'*uint');
curpos=ftell(fd);
fseek(fd,x.wsrt_offset,'bof');
for i=1:x.numlayers
  if ~strcmp(readString(fd,4),'WSRT')
    error('Out of sync trying to read WSRT at offset 0x%x\n', ...
          ftell(fd)-4);
  end
  x.wellsettings(i)=readWSRT(fd,mx);
end
fseek(fd,curpos,'bof');

% Real implementation of readWSRT
function x=readWSRT(fd,mx)
x.name=readString(fd,20);
x.v1=fread(fd,1,'*ushort');
x.conc=fread(fd,1,'double');
x.v2=fread(fd,1,'*ushort');
readfill(fd,2);
x.v2b=fread(fd,1,'*ushort');
x.v3=fread(fd,1,'*ushort');
x.dye=fread(fd,1,'*ushort');
readfill(fd,2);
x.index=fread(fd,1,'*ushort');


function x=readSpectral(fd,nsamples)
if nsamples==0
  x=[];
end
for i=1:nsamples
  x(i)=readSample(fd);
end

function x=readSample(fd)
x.bias=fread(fd,1,'*ushort');
x.exposure=10*2^fread(fd,1,'ushort');
x.v1=fread(fd,1,'*float');
x.i3=fread(fd,1,'*ushort');
x.count=fread(fd,1,'*ushort');
x.measnum=fread(fd,1,'*ushort');
x.i4=fread(fd,8,'*ushort');
x.stage=fread(fd,1,'*ushort');
x.step=fread(fd,1,'*ushort');
x.rep=fread(fd,1,'*ushort');
x.temperature1=fread(fd,1,'*float');
x.dwellTime=fread(fd,1,'*ushort');
x.dwellTimeRemaining=fread(fd,1,'*ushort');
x.temperature2=fread(fd,4,'*float');
x.filter1=fread(fd,1,'*ushort');
x.filter2=fread(fd,1,'*ushort');
x.measindex=fread(fd,1,'*ushort');
readfill(fd,2);
x.time=fread(fd,1,'double');
x.validFlag=fread(fd,1,'*ushort');
readfill(fd,6);
for j=1:x.count
  x.well(j,1)=fread(fd,1,'*ushort');
  x.v(:,j)=fread(fd,2,'*float');
  x.clipped(j,1)=fread(fd,1,'*ushort');
end

function x=readRMEL(fd,mx)
x.name=readString(fd,20);
x.i1=fread(fd,1,'*uint');
x.i2=fread(fd,1,'*ushort');
readfill(fd,12);
for i=1:x.i1
  if ~strcmp(readString(fd,4),'RMEF')
    error('Out of sync trying to read RMEF at offset 0x%x\n', ...
          ftell(fd)-4);
  end
  x.RMEF(i)=readRMEF(fd,mx);
end

function x=readRMEF(fd,mx)
x.excitation=readString(fd,20);
x.emission=readString(fd,20);
x.response=fread(fd,1,'*double');
readfill(fd,16);

function x=readWFFR(fd,mx)
x.v1=fread(fd,5,'*ushort');
readfill(fd,38);

function x=readWFFP(fd,mx)
x.i1=fread(fd,1,'*ushort');
readfill(fd,2);
x.dye=readString(fd,20);
x.excitation=readString(fd,20);
x.emission=readString(fd,20);
x.cnt=fread(fd,1,'*ushort');
readfill(fd,40);
x.val=fread(fd,x.cnt,'*float');

function x=readPIFB(fd,mx)
for i=1:5
  x.version{i}=readString(fd,25);
end
readfill(fd,100);
x.pc=readString(fd,15);
readfill(fd,100);
x.windowsversion=fread(fd,8,'*ushort');
x.servicepack=readString(fd,14);
readfill(fd,218);
x.v=fread(fd,6,'*uint');
readfill(fd,100);
x.location=readString(fd,256);
readfill(fd,400);
for i=1:9
  x.labels{i}=readString(fd,25);
end
for i=1:mx.TCPS.nsteps
  x.step(i)=readStep(fd);
end
x.i2=fread(fd,2,'*ushort');
x.v1=fread(fd,1,'double');
x.v2=fread(fd,1,'double');
x.i3=fread(fd,2,'*ushort');

function x=readStep(fd)
x.cycle=fread(fd,1,'*uint');
x.repeat=fread(fd,1,'*ushort');
x.step=fread(fd,1,'*ushort');
x.v2=fread(fd,1,'*float');
x.duration=fread(fd,1,'*ushort');
x.temperature=fread(fd,1,'*float');
x.gradient=fread(fd,1,'*float');
x.i3=fread(fd,5,'*ushort');

function x=readUDLP(fd,mx)
x.layer=fread(fd,1,'*uint');
x.cnt=fread(fd,2,'*ushort');
for i=1:length(x.cnt)
  x.enabledwells(i,:)=fread(fd,x.cnt(i),'*ushort');
end
x.i1=fread(fd,6,'*uint');
x.i2=fread(fd,1,'*ushort');
x.i2=fread(fd,1,'*ushort');
x.v1=fread(fd,1,'double');
x.v2=fread(fd,1,'double');

function x=readTCCP(fd,mx)
x.i1=fread(fd,8,'*ushort');
readfill(fd,6);
x.i2=fread(fd,4,'*ushort');
readfill(fd,26);
x.i2=fread(fd,3,'*ushort');
readfill(fd,70);

function x=readPSPR(fd,mx)
x.c1=readString(fd,1);
x.numwirt=fread(fd,1,'*uint');
for i=1:x.numwirt
  if ~strcmp(readString(fd,4),'WIRT')
    error('Out of sync trying to read WIRT at offset 0x%x\n',ftell(fd)-4);
  end
  x.WIRT(i)=readWIRT(fd,mx);
end
while (1)
  code=readString(fd,4);
  if strcmp(code,'WSRT')
    % Called by readWIRT() 
    fread(fd,44,'uint8');
  else
    % Spectral data - hacked to find this block immediately after
    % last WSRT block
    fseek(fd,-4,'cof');  % Back up to reread the code
    break
  end
end
readfill(fd,102);
