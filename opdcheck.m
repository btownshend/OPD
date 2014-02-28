% Check out OPD data to see if it makes sense
function opdcheck(v,samps)
clf;
for i=1:length(v.WIRT)
  pos=v.WIRT(i).platepos;
  wellnames{i}=sprintf('%c%d',floor(double(pos)/12)+'A',mod(pos,12)+1);
end
if nargin<2
  sampsel=1:length(v.WIRT)
else
  sampsel=ismember([v.WIRT.platepos],samps);
end
  
pnum=1;
for stage=1:length(v.stageavg)
  avg=v.stageavg{stage};
  delta=(avg.temp(end)-avg.temp(1))>1;
  for dye=1:size(avg.scaled,2)
    s=squeeze(avg.scaled(:,dye,sampsel));
    subplot(length(v.stageavg),size(avg.scaled,2),pnum);
    pnum=pnum+1;
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
      if nargin>=2 || length(wellnames)<20
        legend(wellnames{sampsel},'Location','EastOutside');
      end
    end
  end
end
