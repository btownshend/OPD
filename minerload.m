function [r,z]=minerload(job,s)
if nargin<2
  if strncmp(job,'file://',7)
    s=urlread(job);
  else
    s=urlread(sprintf('http://ewindup.info/miner/Results/Miner_%s_Analyzed_Data.txt',job));
  end
end
line=1;
z={};
while ~isempty(s)
  [tok,s]=strtok(s,10);
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
