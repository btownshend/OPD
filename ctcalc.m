% Calculate the Ct's
function opd=ctcalc(opd,varargin)
defaults=struct('basecycles',2:8,'thresh',nan,'doplot',false,'debug',false,'samps',[],'maxcterror',0.4,'fulow',[],'fuhigh',[]);
args=processargs(defaults,varargin);
wellnms=wellnames(opd);
if isempty(args.samps)
  fu=squeeze(opd.avg.scaled);
else
  w=[];
  for i=1:length(args.samps)
    f=find(strcmp(args.samps{i},wellnms));
    if ~isempty(f)
      w(end+1)=f;
    end
  end
  sampsel=w;
  fu=squeeze(opd.avg.scaled(:,:,sampsel));
end

if isempty(args.fulow) || isnan(args.fulow)
  if isempty(args.thresh) || isnan(args.thresh)
    args.fulow=4*mean(std(fu(args.basecycles,:),1));
  else
    args.fulow=args.thresh/2;
  end
end
if isempty(args.thresh) || isnan(args.thresh)
  args.thresh=args.fulow*2;
end
if isempty(args.fuhigh) || isnan(args.fuhigh)
  args.fuhigh=args.thresh*4;
end

if (args.debug)
  fprintf('Thresh=%.1f, FU(low)=%.1f, FU(high)=%.1f\n', args.thresh, args.fulow, args.fuhigh);
end

% Loop over each trace
fuexp=nan(size(fu));
ct=nan(1,size(fu,2));
for i=1:size(fu,2)
  basecycles=args.basecycles;
  for baseredo=1:3
    baseline=prctile(fu(basecycles,i),25);
    fu(:,i)=fu(:,i)-baseline;
    estart=find(fu(:,i)>args.fulow & (1:size(fu,1))'>max(basecycles),1);
    if ~isempty(estart)
      elast=find(fu(estart:end,i)>args.fuhigh,1);
      if isempty(elast)
        elast=size(fu,1);
      else
        elast=elast+estart-2;
      end
      sel=false(size(fu,1),1);
      sel(estart:elast)=true;
      sel(fu(:,i)<0)=false;
      cntr=1:size(fu,1);
      fit=polyfit(log10(fu(sel,i)),cntr(sel)',1);
      ctpred=polyval(fit,log10(fu(sel,i)));
      rmserror=sqrt(mean((cntr(sel)'-ctpred).^2));
      ct(i)=polyval(fit,log10(args.thresh));
      fuexp(sel,i)=fu(sel,i);
      fuexp(~sel,i)=nan;
    else
      estart=nan;
      elast=nan;
      fit=[nan,nan];
      rmserror=nan;
    end
    if args.debug
      fprintf('Sample %s: baseline=%.1f [%d-%d], max=%.1f, estart=%d, elast=%d, fit=(%f,%f), ct=%.1f rmserr=%.1f\n', wellnms{i}, baseline, min(basecycles), max(basecycles), max(fu(:,i)), estart, elast,fit,ct(i),rmserror);
    end
    if baseredo==1 && isfinite(estart)
      % Reset the base cycles to just before the area of interest
      basecycles=max(min(args.basecycles),estart-6):(estart-1);
    end
  end
  eff=10^(1/fit(1));
  if rmserror>args.maxcterror  || eff<1.3 || eff>2.5
    fprintf('Bad fit: Sample %s: baseline=%.1f [%d-%d], max=%.1f, estart=%d, elast=%d, eff=%.2f, (acceptable is 1.3-2.5), ct=%.1f rmserr=%.1f (max=%.1f)\n', wellnms{i}, baseline, min(basecycles), max(basecycles), max(fu(:,i)), estart, elast,eff,ct(i),rmserror,args.maxcterror);
    ct(i)=nan;
  end
  opd.fit(i,:)=fit;
end
opd.ct=ct;

% Create a full grid version of ct
if isempty(args.samps)
  w=wellnames(opd);
else
  w=args.samps;
end
opd.ctgrid=nan(8,12);
for j='A':'H'
  for i=1:12
    ind=find(strcmp(w,sprintf('%c%d',j,i)));
    if ~isempty(ind)
      opd.ctgrid(j-'A'+1,i)=ct(ind);
    end
  end
end

if args.doplot
  setfig(opd.filename); clf;
  subplot(211);
  plot(fu);
  hold on;
  c=axis;
  plot([c(1),c(2)],args.fulow*[1,1],':');
  plot([c(1),c(2)],args.fuhigh*[1,1],':');
  plot(args.basecycles(1)*[1,1],[c(3),c(4)],':');
  plot(args.basecycles(end)*[1,1],[c(3),c(4)],':');
  if ~isempty(args.samps) && length(args.samps)<20
    legend(args.samps,'Location','EastOutside');
  end
  title(opd.filename);
  subplot(212);
  semilogy(fu,'y');
  hold on;
  semilogy(fu(:,~isfinite(ct)),'m');
  semilogy(fuexp(:,isfinite(ct)),'.-');
  plot(ct,0*ct+args.thresh,'o');
  plot([c(1),c(2)],args.fulow*[1,1],':');
  plot([c(1),c(2)],args.fuhigh*[1,1],':');
  c=axis; c(3)=args.fulow/10; axis(c);
  w=wellnames(opd);
  if ~isempty(args.samps)
    w=args.samps;
  end
  % if length(w)==length(ct)
  %   for i=1:length(ct)
  %     fprintf('%s\t%.2f\n',w{i},ct(i));
  %   end
  % else
  %   for i=1:length(ct)
  %     fprintf('%.2f\n',ct(i));
  %   end
  % end  

  fprintf('\t');
  for i=1:12
    fprintf('%d\t',i);
  end
  fprintf('\n');
  for j='A':'H'
    fprintf('%c\t',j);
    for i=1:12
      ind=find(strcmp(w,sprintf('%c%d',j,i)));
      if isempty(ind)
        fprintf('\t');
      else
        fprintf('%.2f\t',ct(ind));
      end
    end
    fprintf('\n');  
  end
end
