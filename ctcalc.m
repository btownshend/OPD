% Calculate the Ct's
function opd=ctcalc(opd,varargin)
defaults=struct('basecycles',2:8,'thresh',100);
args=processargs(defaults,varargin);
debug=0;
if nargin<2
  args.basecycles=2:8;
end
fu=squeeze(opd.avg.scaled);
% Loop over each trace
fuexp=nan(size(fu));
ct=nan(1,size(fu,2));
emin=12; emax=12;
for i=1:size(fu,2)
  baseline=mean(fu(args.basecycles,i));
  fu(:,i)=fu(:,i)-baseline;
  estart=find(fu(:,i)>baseline/5,1);
  if ~isempty(estart)
    elast=max(find(fu(:,i)<baseline*0.8));
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
  if debug
    fprintf('Sample %d: baseline=%.1f, max=%.1f, estart=%d, elast=%d, fit=(%f,%f), ct=%.1f\n', i, baseline, max(fu(:,i)), estart, elast,fit,ct(i));
  end
end
opd.ct=ct;
