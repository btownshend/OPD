% Print out info appropriate for pasting into Miner
function minerdump(v,samps)
if nargin<2
  sampsel=1:length(v.WIRT);
else
  sampsel=find(ismember([v.WIRT.platepos],samps));
end

fprintf('Cycle\t');
for ii=1:length(sampsel)
  i=sampsel(ii);
  pos=v.WIRT(i).platepos;
  wellnames{i}=sprintf('%c%d',floor(double(pos)/12)+'A',mod(pos,12)+1);
  if i<length(sampsel)
    fprintf('%s\t',wellnames{i});
  else
    fprintf('%s\n',wellnames{i});
  end      
end

sc=v.avg.scaled(:,:,sampsel);
for i=1:size(sc,1)
  fprintf('%f\t',v.avg.cycle(i));
  for j=1:size(sc,3)-1
    fprintf('%f\t',sc(i,1,j));    
  end
  fprintf('%f\n', sc(i,1,end));
end
  