% Plot the Ct's of the OPD data in z.avg.scaled
function ct=ctplot(z,basecycles)
fu=squeeze(z.avg.scaled);
if nargin<2
  basecycles=2:max(3,min(8,size(fu,1)));
end
thresh=100;
if size(fu,1)<max(basecycles)
  error('Not enough cycles -- need %d, have %d\n', max(basecycles), size(fu,1));
end
% Loop over each trace
fuexp=nan(size(fu));
ct=nan(1,size(fu,2));
emin=12; emax=12;
for i=1:size(fu,2)
  baseline=mean(fu(basecycles,i));
  fu(:,i)=fu(:,i)-baseline;
  estart=find(fu(:,i)>baseline/5,1);
  if ~isempty(estart)
    elast=max(find(fu(:,i)<baseline*0.8));
    fit=polyfit(log10(fu(estart:elast,i)),(estart:elast)',1);
    ct(i)=polyval(fit,log10(thresh));
    fuexp(estart:elast,i)=fu(estart:elast,i);
    emin=min(estart,emin);
    emax=max(elast,emax);
  else
    estart=nan;
    elast=nan;
    fit=[nan,nan];
  end
  fprintf('Sample %d: baseline=%.1f, max=%.1f, estart=%d, elast=%d, fit=(%f,%f), ct=%.1f\n', i, baseline, max(fu(:,i)), estart, elast,fit,ct(i));
end
setfig('ctplot'); clf;
subplot(211);
plot(fu);
subplot(212);
semilogy(fuexp);
hold on;
plot(ct,0*ct+thresh,'o');
if any(isfinite(fuexp(:)))
  axis([emin,emax,nanmin(fuexp(:)),nanmax(fuexp(:))]);
end
w=wellnames(z);
if length(w)==length(ct)
  for i=1:length(ct)
    fprintf('%s\t%.2f\n',w{i},ct(i));
  end
else
  for i=1:length(ct)
    fprintf('%.2f\n',ct(i));
  end
end  

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