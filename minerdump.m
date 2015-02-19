% Print out info appropriate for pasting into Miner
% Fixlib- true to remove declines in trace after a peak, such as occurs with over-amplified library -- BUT, doesn't seem to help MINER in fitting these... (2/2/15)
function minerdump(v,samps,fixlib)
if nargin<3
  fixlib=false;
end

if nargin<2 || isempty(samps)
  sampsel=1:length(v.WIRT);
else
  sampsel=find(ismember([v.WIRT.platepos],samps));
end

fprintf('Cycle\t');
for ii=1:length(sampsel)
  i=sampsel(ii);
  pos=v.WIRT(i).platepos;
  wellnames{i}=sprintf('%c%d',floor(double(pos)/12)+'A',mod(pos,12)+1);
  if ii<length(sampsel)
    fprintf('%s\t',wellnames{i});
  else
    fprintf('%s\n',wellnames{i});
  end      
end

sc=v.avg.scaled(:,:,sampsel);
for i=1:size(sc,1)
  fprintf('%f\t',v.avg.cycle(i)+0.63);
  for j=1:size(sc,3)
    if fixlib
      t=max(sc(1:i,1,j));
    else
      t=sc(i,1,j);
    end
    if j<size(sc,3)
      fprintf('%f\t',t);    
    else
      fprintf('%f\n', t);
    end
  end
end
  