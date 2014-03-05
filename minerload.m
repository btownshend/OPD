function [r,z]=minerload(job,s)
if nargin<2
  if job(1)>='0' && job(1)<='9'
    s=urlread(sprintf('http://ewindup.info/miner/Results/Miner_%s_Analyzed_Data.txt',job));
  else
    fd=fopen(job,'r');
    if fd==-1
      error('Unable to open %s for reading\n', job);
    end
    s=fread(fd,inf,'*char')';
    fclose(fd);
  end
end
line=1;
z={};
while ~isempty(s)
  [tok,s]=strtok(s,10);
  tok=tok(tok~=13);   % Remove cr
  cnt=1;
  while ~isempty(tok)
    [elem,tok]=strtok(tok,9);
    z{line,cnt}=elem;
    cnt=cnt+1;
  end
  line=line+1;
end
r=struct();
for i=1:size(z,1)
  field=strrep(z{i,1},':','');
  field(~isstrprop(field,'alphanum'))='_';
  v=z(i,2:end);
  if ~strcmp(field,'SampleNames')
    % Convert to numeric
    v=strrep(v,'Error!','NaN');
    vv=[];
    for k=1:length(v)
      vv(k)=str2num(v{k});
    end
    r.(field)=vv;
  else
    r.(field)=v;
  end
  if strcmp(field,'CT')
    break;
  end
end
r.job=job;
