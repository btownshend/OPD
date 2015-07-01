% Calculate the Ct's
function opd=ctcalc(opd,varargin)
defaults=struct('basecycles',2:8,'thresh',nan,'doplot',false,'debug',false,'samps',[]);
args=processargs(defaults,varargin);
if isempty(args.samps)
  fu=squeeze(opd.avg.scaled);
else
  wellnms=wellnames(opd);
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

fulow=4*mean(std(fu(args.basecycles,:),1));
if isnan(args.thresh)
  args.thresh=fulow*2;
else
  fulow=args.thresh/2;
end
fuhigh=args.thresh*4;
if (args.debug)
  fprintf('Thresh=%.1f, FU(low)=%.1f, FU(high)=%.1f\n', args.thresh, fulow, fuhigh);
end

% Loop over each trace
fuexp=nan(size(fu));
ct=nan(1,size(fu,2));
emin=12; emax=12;
for i=1:size(fu,2)
  basecycles=args.basecycles;
  for baseredo=1:2
    baseline=prctile(fu(basecycles,i),25);
    fu(:,i)=fu(:,i)-baseline;
    estart=find(fu(:,i)>fulow & (1:size(fu,1))'>max(basecycles),1);
    if ~isempty(estart)
      elast=find(fu(estart:end,i)>fuhigh,1);
      if isempty(elast)
        elast=size(fu,1);
      else
        elast=elast+estart-1;
      end
      sel=false(size(fu,1),1);
      sel(estart:elast)=true;
      sel(fu(:,i)<0)=false;
      cntr=1:size(fu,1);
      fit=polyfit(log10(fu(sel,i)),cntr(sel)',1);
      ct(i)=polyval(fit,log10(args.thresh));
      fuexp(sel,i)=fu(sel,i);
      emin=min(estart,emin);
      emax=max(elast,emax);
    else
      estart=nan;
      elast=nan;
      fit=[nan,nan];
    end
    if args.debug
      fprintf('Sample %d: baseline=%.1f [%d-%d], max=%.1f, estart=%d, elast=%d, fit=(%f,%f), ct=%.1f\n', i, baseline, min(basecycles), max(basecycles), max(fu(:,i)), estart, elast,fit,ct(i));
    end
    if baseredo==1 && isfinite(estart)
      basecycles=max(min(args.basecycles),estart-8):(estart-2);
    end
  end
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
  setfig('ctplot'); clf;
  subplot(211);
  plot(fu);
  hold on;
  c=axis;
  plot([c(1),c(2)],fulow*[1,1],':');
  plot([c(1),c(2)],fuhigh*[1,1],':');
  plot(args.basecycles(1)*[1,1],[c(3),c(4)],':');
  plot(args.basecycles(end)*[1,1],[c(3),c(4)],':');
  if ~isempty(args.samps) && length(args.samps)<20
    legend(args.samps,'Location','EastOutside');
  end
  subplot(212);
  semilogy(fu,'y');
  hold on;
  semilogy(fuexp,'.-');
  hold on;
  plot(ct,0*ct+args.thresh,'o');
  if any(isfinite(fuexp(:)))
    % axis([emin,emax,nanmin(fuexp(:)),nanmax(fuexp(:))]);
  end
  c=axis; c(3)=10; axis(c);
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
