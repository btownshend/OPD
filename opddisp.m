% opddisp - display an OPD structure
function opddisp(x,prefix,maxdisp)
global lastprefix
if nargin<2
  prefix='';
end
if nargin<3
  maxdisp=3
end
if length(prefix)==0
  lastprefix='';
end
if iscell(x)
  for i=1:min(length(x),maxdisp)
    opddisp(x{i},sprintf('%s[%d]',prefix,i),maxdisp);
  end
  if length(x)>maxdisp
    fprintf('%s[%d..%d] ...\n',pf(prefix),maxdisp+1,length(x));
  end
elseif length(x)>1 && ~ischar(x)
  for i=1:min(length(x),maxdisp)
    opddisp(x(i),sprintf('%s[%d]',prefix,i),maxdisp);
  end
  if length(x)>maxdisp
    fprintf('%s[%d..%d] ...\n',pf(prefix),maxdisp+1,length(x));
  end
elseif isstruct(x)
  f=fieldnames(x);
  for i=1:length(f)
    opddisp(x.(f{i}),[prefix,'.',f{i}],maxdisp);
  end
else
  fprintf('%s:\t',pf(prefix));
  if isfloat(x)
    fprintf('%f ',x);
  elseif isinteger(x)
    fprintf('%d ',x);
  elseif ischar(x)
    fprintf('"%s" ',x);
  else
    fprintf('UNKNOWN TYPE: %s',class(x));
  end
  fprintf('\n');
end

function p=pf(prefix)
global lastprefix
if length(prefix)>0 && prefix(1)=='.'
  prefix=prefix(2:end);
end
nmatch=0;
for i=1:min(length(prefix),length(lastprefix))
  if prefix(i)~=lastprefix(i)
    break;
  end
  nmatch=i;
end
while nmatch>0 && nmatch~=length(prefix) && prefix(nmatch+1)~='.' && prefix(nmatch+1)~='['
  nmatch=nmatch-1;
end
if nmatch==0
  p=prefix;
else
  p=[blanks(nmatch),prefix(nmatch+1:end)];
end
lastprefix=prefix;
