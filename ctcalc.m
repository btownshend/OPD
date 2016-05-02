% Calculate the Ct's
function opd=ctcalc(opd,varargin)
defaults=struct('basecycles',2:8,'thresh',nan,'doplot',false,'debug',false,'samps',[],'maxcterror',0.4,'fulow',[],'fuhigh',[],'showall',false);
args=processargs(defaults,varargin);
wellnms=wellnames(opd);
if isempty(args.samps)
  fu=squeeze(opd.avg.scaled);
  if ~args.showall
    levdiff=squeeze(max(opd.avg.scaled(:,1,:))-min(opd.avg.scaled(:,1,:)));
    mindiff=max(levdiff)*0.2;
    sampsel=find(levdiff>mindiff);
    if length(sampsel)<size(fu,2)
      fprintf('Displaying only %d/%d wells that change by more than %.0f over course of run\n', length(sampsel), size(fu,2), mindiff);
      fu=fu(:,sampsel);
      wellnms=wellnms(sampsel);
    end
  end
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
  wellnms=args.samps;
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
  [baseline,poly,lastbaseline,basestd]=findbase(fu(:,i),'basecycles',args.basecycles,'debug',args.debug,'doplot',args.doplot&&size(fu,2)==1);
  fu(:,i)=fu(:,i)-baseline;
  if args.fulow<basestd*3
    fprintf('ctcalc: Sample %s: Warning: fulow (%.0f) < baseline stdev*3 (=%.1f)\n', wellnms{i}, args.fulow, 3*basestd);
  end
  estart=find(fu(:,i)>args.fulow & (1:size(fu,1))'>lastbaseline,1);
  if ~isempty(estart)
    elast=find(fu(estart:end,i)>min(args.fuhigh,max(fu(:,i))/2),1);
    if isempty(elast)
      elast=size(fu,1);
    else
      elast=elast+estart-2;
    end
    sel=false(size(fu,1),1);
    sel(estart:elast)=true;
    sel(fu(:,i)<0)=false;
    if sum(sel)>=2
      cntr=1:size(fu,1);
      fit=polyfit(log10(fu(sel,i)),cntr(sel)',1);
      ctpred=polyval(fit,log10(fu(sel,i)));
      rmserror=sqrt(mean((cntr(sel)'-ctpred).^2));
      ct(i)=polyval(fit,log10(args.thresh));
      fuexp(sel,i)=fu(sel,i);
      fuexp(~sel,i)=nan;
    else
      fprintf('No exponential range found for sample %s: baseline=%.2fx+%.1f [%d-%d], max=%.1f, estart=%d\n', wellnms{i}, poly, min(args.basecycles), lastbaseline, max(fu(:,i)), estart);
      estart=nan;
      elast=nan;
      fit=[nan,nan];
      rmserror=nan;
    end
  else
    fprintf('No exponential range found for sample %s: baseline=%.2fx+%.1f [%d-%d], max=%.1f\n', wellnms{i}, poly, min(args.basecycles), lastbaseline, max(fu(:,i)));
    estart=nan;
    elast=nan;
    fit=[nan,nan];
    rmserror=nan;
  end
  if args.debug
    fprintf('Sample %s: baseline=%.2fx+%.1f [%d-%d], max=%.1f, estart=%d, elast=%d, fit=(%f,%f), ct=%.1f rmserr=%.1f\n', wellnms{i}, poly, min(args.basecycles),lastbaseline, max(fu(:,i)), estart, elast,fit,ct(i),rmserror);
  end
  eff=10^(1/fit(1));
  if rmserror>args.maxcterror  || eff<1.3 || eff>2.5
    fprintf('Bad fit: Sample %s: baseline=%.2fx+%.1f [%d-%d], max=%.1f, estart=%d, elast=%d, eff=%.2f, (acceptable is 1.3-2.5), ct=%.1f rmserr=%.1f (max=%.1f)\n', wellnms{i}, poly, min(args.basecycles), lastbaseline, max(fu(:,i)), estart, elast,eff,ct(i),rmserror,args.maxcterror);
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
    ind=find(strcmp(wellnms,sprintf('%c%d',j,i)));
    if ~isempty(ind)
      opd.ctgrid(j-'A'+1,i)=ct(ind);
    end
  end
end

if args.doplot
  setfig(opd.filename); clf;
  subplot(311);
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

  subplot(312);
  semilogy(fu,'y');
  hold on;
  semilogy(fu(:,~isfinite(ct)),'m');
  semilogy(fuexp(:,isfinite(ct)),'.-');
  plot(ct,0*ct+args.thresh,'o');
  plot([c(1),c(2)],args.fulow*[1,1],':');
  plot([c(1),c(2)],args.fuhigh*[1,1],':');
  %ylim=get(gca,'YLim');
  ylim=[args.fulow/10,max(args.fuhigh,max(fu(:)))*1.1]; 
  set(gca,'YLim',ylim);
  
  subplot(313);
  failed=~isfinite(ct);
  plot(fu(:,~isfinite(ct)));
  hold on;
  plot([c(1),c(2)],args.fulow*[1,1],':k');
  plot([c(1),c(2)],args.fuhigh*[1,1],':k');
  %c=axis; c(3)=args.fulow/10; axis(c);
  % if length(wellnms)==length(ct)
  %   for i=1:length(ct)
  %     fprintf('%s\t%.2f\n',wellnms{i},ct(i));
  %   end
  % else
  %   for i=1:length(ct)
  %     fprintf('%.2f\n',ct(i));
  %   end
  % end  

  fprintf('  ');
  for i=1:12
    fprintf('    %2d ',i);
  end
  fprintf('\n');
  for j='A':'H'
    fprintf('%c ',j);
    for i=1:12
      ind=find(strcmp(wellnms,sprintf('%c%d',j,i)));
      if isempty(ind)
        fprintf('       ');
      else
        fprintf(' %5.2f ',ct(ind));
      end
    end
    fprintf('\n');  
  end
end
