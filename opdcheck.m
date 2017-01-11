% Check out OPD data to see if it makes sense
function opdcheck(v,samps,varargin)
% showall - show even wells that appear to be empty/unused (default: false)
defaults=struct('firststage',false,'basecycles',[],'thresh',[],'clf',true,'showall',false);
args=processargs(defaults,varargin);
wellnms=wellnames(v);
if nargin<2 || isempty(samps)
  sampsel=1:length(v.WIRT);
  if ~args.showall
    levdiff=squeeze(diff(v.avg.scaled([1,end],1,:)));
    mindiff=max(levdiff)*0.1;
    sampsel=intersect(sampsel,find(levdiff>mindiff));
    if length(sampsel)<length(v.WIRT)
      fprintf('Displaying only %d/%d wells that change by more than %.0f over course of run\n', length(sampsel), length(v.WIRT), mindiff);
    end
  end
elseif iscell(samps)
  w=[];
  for i=1:length(samps)
    f=find(strcmp(samps{i},wellnms));
    if ~isempty(f)
      w(end+1)=f;
    end
  end
  sampsel=w;
else
  sampsel=ismember([v.WIRT.platepos],samps);
end
  
pnum=1;
if args.clf
  clf;
end
for stage=1:length(v.stageavg)
  avg=v.stageavg{stage};
  delta=(avg.temp(end)-avg.temp(1))>1;
  for dye=1:size(avg.scaled,2)
    s=squeeze(avg.scaled(:,dye,sampsel));
    if ~args.firststage && (length(v.stageavg)>1 || size(avg.scaled,2)>1)
      subplot(length(v.stageavg),size(avg.scaled,2),pnum);
      pnum=pnum+1;
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
      if ~isempty(args.basecycles) && stage==1 && size(s,1)>max(args.basecycles)
        % Subtract out baseline
        for k=1:size(s,2)
          s(:,k)=s(:,k)-mean(s(args.basecycles,k));
        end
      end
      plot(avg.cycle,s,'-');
      ylabel('Fluorescence');
      xlabel('Cycle');
      if ~isempty(args.basecycles) && stage==1
        hold on;
        c=axis;
        plot(min(args.basecycles-.5)*[1,1],c(3:4),':k');
        plot(max(args.basecycles+.5)*[1,1],c(3:4),':k');
        axis(c);
      end
      if ~isempty(args.thresh)
        hold on;
        c=axis;
        plot(c(1:2),args.thresh*[1,1],':k');
        axis(c);
      end
    end
    title(sprintf('%s - %s',v.filename,v.WFFP(dye).dye));
    if stage==1
      if nargin>=2 || length(wellnms)<20
        legend(wellnms{sampsel},'Location','EastOutside');
        for i=1:size(s,2)
          text(avg.cycle(end)+0.5,s(end,i),wellnms{sampsel(i)});
        end
        c=axis;
        c(2)=c(2)+1.5;
        axis(c);
      end
    end
  end
  if args.firststage
    break;
  end
end
