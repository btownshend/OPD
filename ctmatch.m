% Align the fluorescence plots using variable (LO,CT,HI) for each plot to minimize MSE
function opd=ctmatch(opd,varargin)
defaults=struct('doplot',false,'debug',false,'samps',[]);
args=processargs(defaults,varargin);
if isempty(args.samps)
  fu=squeeze(opd.avg.scaled);
else
  fu=squeeze(opd.avg.scaled(:,:,args.samps));
end


ref=fu(:,1);
ref=alignfu(ref,[ref(1),ref(end),0,1]);
setfig('ctmatch.1'); clf;
plot(ref);
hold on;
plot(alignfu(ref,[0,1,0,1]));
plot(alignfu(ref,[0,2,0,1]));
plot(alignfu(ref,[0,1,1,1]));
plot(alignfu(ref,[0,1,-1,1]));
match=matchfu(ref,ref)
maxerr=0.2;

fit2=matchfu(fu(:,2),ref);
fprintf('Sample %d: lo=%.0f, hi=%.0f, ct=%.2f\n', 2,fit2);
fitted2=alignfu(fu(:,2),fit2);
setfig('ctmatch.1.1');clf;
plot(ref);
hold on;
plot(fitted2);

for j=1:2
  % Loop over each trace
  fit=nan(size(fu,2),4);
  err=nan(size(fu,2),1);
  for i=1:size(fu,2)
    fit(i,:)=matchfu(fu(:,i),ref);
    fitted(:,i)=alignfu(fu(:,i),fit(i,:));
    err(i)=norm(fitted(:,i)-ref');
    if args.debug
      fprintf('Sample %d: lo=%.0f, hi=%.0f, ct=%.2f, hill=%.2f err=%.2f\n', i,fit(i,:),err(i));
    end
  end
  setfig(sprintf('ctmatch.%d',j+1)); clf;
  sel=err<maxerr;
  subplot(211);
  plot(fitted(:,sel));
  subplot(212);
  plot(err);
  ref=median(fitted(:,sel)');
end
return

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
    w=w{args.samps};
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

% Fit to the reference with fit parameters lo,hi,ct
function fit=matchfu(fu,ref)
xinit=[min(fu),max(fu),0.1,1];
opts=optimset('MaxFunEvals',10000,'MaxIter',10000,'TolX',1e-10,'TolFun',1e-10,'Display','notify');
fit=fminsearch(@(z) norm(alignfu(fu,z)-ref), xinit,opts);

function scaledfu=alignfu(fu,fit)
if fit(1)<0
  fit(1)=0;
end
maxshift=10;
if fit(3)<-maxshift
  fit(3)=-maxshift;
end
if fit(3)>maxshift
  fit(3)=maxshift;
end
fit(4)=min(1.1,max(0.9,fit(4)));
scaledfu=interp1(1:length(fu),(fu-fit(1))/(fit(2)-fit(1)),(1:length(fu))*fit(4)-fit(3),'pchip','extrap');
%scaledfu(1:ceil(fit(3)))=(fu(1)-fit(1))/(fit(2)-fit(1));
%scaledfu(end-ceil(-fit(3)+1):end)=(fu(end)-fit(1))/(fit(2)-fit(1));