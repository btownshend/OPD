% Plot melt curves
function [ut,sf,sfg]=opdmelt(v,wells,varargin)
defaults=struct('gelplot',false,'wellnames',[]);
args=processargs(defaults,varargin);
if nargin<2 || isempty(wells)
  wells=1:size(v.all.scaled,3);
end
for i=1:length(v.WIRT)
  pos=v.WIRT(i).platepos;
  wellnames{i}=sprintf('%c%d',floor(double(pos)/12)+'A',mod(pos,12)+1);
end
if iscell(wells)
  w=[];
  for i=1:length(wells)
    nw=find(strcmp(wells{i},wellnames));
    if isempty(nw)
      fprintf('Well %s not present in OPD file -- ignoring\n');
    else
      w(end+1)=nw;
    end
  end
  wells=w;
end
if isempty(args.wellnames)
  args.wellnames=wellnames(wells);
elseif ~iscell(args.wellnames)
  error('Wellnames arg must be a cell array of length %d',length(wells));
elseif length(args.wellnames)~=length(wells)
  error('Have %d wells, but %d well names',length(wells), length(args.wellnames));
end
% Select melt cycles (assume it is the one with 100 or 200 repeats
stage=v.PIFB.step(ismember([v.PIFB.step.repeat],[100,150,200])).cycle;
sel=(v.all.stage==stage);
if sum(sel)<100
  fprintf('No melt data\n');
  ut=[];
  return;
end
temp=[v.data(sel).temperature1];
trange=max(temp)-min(temp);
npoints=sum(sel);
%fprintf('Taking melt data from cycle %d: npoints=%d, temp=%.1f-%.1f\n',stage,npoints,min(temp),max(temp));
%setfig('response');clf;
%freqz(filtb,filta);
rtemp=round(temp);
ut=unique(rtemp);
[filtb,filta]=ellip(4,0.1,60,double(trange/length(ut)));
%filtb=fir1(10,double(trange/length(ut)/4));filta=[1 0];
%setfig('melt');clf;
for dye=1:size(v.all.scaled,2)
  s=squeeze(v.all.scaled(sel,dye,wells));
  if size(v.all.scaled,2)>1
    subplot(size(v.all.scaled,2),1,dye);
  end
  cycle=v.all.cycle(sel);
  sfg=[];
  sf=[];
  sfgi=[];
  alltemp=min(ut):.1:max(ut);
  for k=1:size(s,2)
    for t=1:length(ut)
      cur=rtemp==ut(t);
      sf(t,k)=mean(s(cur,k));
    end
  end
  % Compute slope of curve over 58-68C range (fluor efficiency loss)
  effsel=ut>=60 & ut<=68;
  effsel=[1,2,length(ut)-1,length(ut)];
  for k=1:size(s,2)
    pfit=polyfit(ut(effsel),sf(effsel,k)',1);
    loss=polyval(pfit,ut)';
    sf(:,k)=sf(:,k)./loss;
    sf(sf(:,k)<0,k)=0;
    sf(sf(:,k)>1,k)=1;
    sfg(:,k)=gradient(filtfilt(filtb,filta,sf(:,k)));
    sfgi(:,k)=spline(ut,sfg(:,k),alltemp);
  end
  %plot(ut,-sf,'-');
  if args.gelplot
    xx=sfgi;
    xx(end+1,:)=nan;
    xx(:,end+1)=nan;
    pcolor(xx);
    shading flat;
    set(gca,'XTick',(1:length(args.wellnames))+0.5);
    set(gca,'XTickLabel',args.wellnames);
    set(gca,'XTickLabelRotation',45);
    tvals=ceil(alltemp(1)/5)*5:5:alltemp(end);
    tpos=interp1(alltemp,1:size(sfgi,1),tvals);
    set(gca,'YTick',tpos);
    set(gca,'YTickLabel',arrayfun(@(z) sprintf('%.0f',z),tvals,'UniformOutput',false));
    ylabel('Temp (C)');
  else
    plot(alltemp,-sfgi,'-');
    xlabel('Temperature');
    ylabel('Fraction melting/deg');
    %fprintf('Length(wells)=%d\n',length(wells));
    if length(wellnames(wells))<=20
      legend(wellnames{wells},'Location','NorthWest');
    end
  end
  title(sprintf('%s - %s',v.filename,v.WFFP(dye).dye));
end

