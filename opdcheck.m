% Check out OPD data to see if it makes sense
function opdcheck(v,samps,varargin)
defaults=struct('firststage',false);
args=processargs(defaults,varargin);
wellnms=wellnames(v);
if nargin<2
  sampsel=1:length(v.WIRT);
elseif iscell(samps)
  w=[];
  for i=1:length(samps)
    w(end+1)=find(strcmp(samps{i},wellnms));
  end
  sampsel=w;
else
  sampsel=ismember([v.WIRT.platepos],samps);
end
  
pnum=1;
for stage=1:length(v.stageavg)
  avg=v.stageavg{stage};
  delta=(avg.temp(end)-avg.temp(1))>1;
  for dye=1:size(avg.scaled,2)
    s=squeeze(avg.scaled(:,dye,sampsel));
    if ~args.firststage && (length(v.stageavg)>1 || size(avg.scaled,2)>1)
      subplot(length(v.stageavg),size(avg.scaled,2),pnum);
      pnum=pnum+1;
    end
    if 0 && stage==1 && size(s,1)>10
      % Subtract out baseline
      for k=1:size(s,2)
        s(:,k)=s(:,k)-mean(s(2:7,k));
      end
    end
    if delta
      [b,a]=butter(12,0.3);
      for i=1:size(s,2)
        % Normalize to 1 total fluorescence
        s(:,i)=s(:,i)/s(1,i);
      end
      dF=-filtfilt(b,a,diff(s));
      plot((avg.temp(2:end)+avg.temp(1:end-1))/2,dF,'-');
      c=axis;
      axis([c(1),c(2),0,c(4)]);
      ylabel('-d(Fluorescence)/d(Cycle)');
      xlabel('Temp (C)');
    else
      plot(avg.cycle,s,'-');
      ylabel('Fluorescence');
      xlabel('Cycle');
    end
    title(sprintf('%s - %s',v.filename,v.WFFP(dye).dye));
    if stage==1
      if nargin>=2 || length(wellnms)<20
        legend(wellnms{sampsel},'Location','EastOutside');
      end
    end
  end
  if args.firststage
    break;
  end
end
