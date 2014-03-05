% Calculate the Ct's
function opd=ctcalc2(opd,varargin)
defaults=struct('doplot',false,'debug',false,'samps',[]);
args=processargs(defaults,varargin);
if isempty(args.samps)
  fu=squeeze(opd.avg.scaled);
else
  fu=squeeze(opd.avg.scaled(:,:,args.samps));
end


% Loop over each trace
fuexp=nan(size(fu));
ct=nan(1,size(fu,2));
cpsdm=nan(1,size(fu,2));
emin=12; emax=12;
fuexp=nan(size(fu));
cycle=opd.avg.cycle;
for i=1:size(fu,2)
  [cf,g]=L4P(cycle,fu(:,i),[9000,-2.8,21,2400],[0,-4,0,0],[inf,0,inf,inf]);
  a=cf.A-cf.D;
  b=cf.B;
  x0=cf.C;
  y0=cf.D;
  fuexp(:,i)=y0+a./(1+(cycle./x0).^b);
  if true
    a=6257.1;
    b=-2.78764;
    x0=20.5513;
    y0=2382.25;
    rnoise=33.6852;
    dynthresh=2559.2;
  end
  cpsdm(i)=x0*((sqrt(3*b^2*(b^2-1))-2*(1-b^2))/(b^2+3*b+2))^(1/b);
  endexp(i)=y0+a./(1+(cpsdm(i)./x0).^b);
  cpspe(i)=x0*((a-rnoise)/rnoise)^(1/b);
  cdt(i)=L4Pinv(cf,dynthresh);
  ct(i)=.4411*cpspe(i)+.6417*cpsdm(i);
  ctf=y0+a./(1+(cpspe(i)./x0).^b);
  if args.debug
    fprintf('Sample %d: a=%.1f, b=%.1f, x0=%.1f, y0=%.1f Noise(SPE)=%.4f, EndofExpPhase(SDM)=%.2f, CP(SPE)=%.5f (%.1f), CP(SDM)=%.5f, f(%.4f)=%.4f\n', i, a,b,x0,y0,rnoise,endexp(i),cpspe(i),ctf,cpsdm(i),cdt(i),dynthresh);
  end
end
opd.ct=ct;

if args.doplot
  setfig('ctplot'); clf;
  subplot(311);
  plot(cycle,fu,'.');
  hold on;
  plot(cycle,fuexp,'g');
  c=axis;
  plot(cpspe(1)*[1,1],c(3:4),':');
  plot(cpsdm(1)*[1,1],c(3:4),':');
  title('Raw data');
  subplot(312);
  semilogy(fu-y0,'.');
  hold on;
  semilogy(fuexp-y0,'g');
  c=axis;
  plot(cpspe(1)*[1,1],c(3:4),':');
  plot(cpsdm(1)*[1,1],c(3:4),':');
  subplot(313);
  plot(cycle,fu-fuexp,'g');
  title('Residual');
  
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
keyboard
