% Calculate the Ct's
function opd=ctcalc(opd,varargin)
defaults=struct('basecycles',2:8,'thresh',nan,'doplot',false,'debug',false,'samps',[]);
args=processargs(defaults,varargin);
if isempty(args.samps)
  fu=squeeze(opd.avg.scaled);
else
  fu=squeeze(opd.avg.scaled(:,:,args.samps));
end

fulow=4*mean(std(fu(args.basecycles,:),1));
if isnan(args.thresh)
  args.thresh=fulow*2;
end
fuhigh=args.thresh*4;

% Loop over each trace
fuexp=nan(size(fu));
ct=nan(1,size(fu,2));
emin=12; emax=12;
for i=1:size(fu,2)
  baseline=mean(fu(args.basecycles,i));
  fu(:,i)=fu(:,i)-baseline;
  estart=find(fu(:,i)>fulow,1);
  if ~isempty(estart)
    elast=max(find(fu(:,i)<fuhigh));
    fit=polyfit(log10(fu(estart:elast,i)),(estart:elast)',1);
    ct(i)=polyval(fit,log10(args.thresh));
    fuexp(estart:elast,i)=fu(estart:elast,i);
    emin=min(estart,emin);
    emax=max(elast,emax);
  else
    estart=nan;
    elast=nan;
    fit=[nan,nan];
  end
  if args.debug
    fprintf('Sample %d: baseline=%.1f, max=%.1f, estart=%d, elast=%d, fit=(%f,%f), ct=%.1f\n', i, baseline, max(fu(:,i)), estart, elast,fit,ct(i));
  end
end
opd.ct=ct;


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
  subplot(212);
  semilogy(fu,'y');
  hold on;
  semilogy(fuexp);
  hold on;
  plot(ct,0*ct+args.thresh,'o');
  if any(isfinite(fuexp(:)))
    % axis([emin,emax,nanmin(fuexp(:)),nanmax(fuexp(:))]);
  end
  w=wellnames(opd);
  if ~isempty(args.samps)
    w=w(args.samps);
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
