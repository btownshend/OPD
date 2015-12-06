% Find the linear baseline of a trace
function [bestbaseline,bestpoly,bestlast,beststd]=findbase(x,varargin)
defaults=struct('basecycles',1:5,'doplot',false,'debug',false,'maxslope',0,'minslope',-1,'maxstd',2);
args=processargs(defaults,varargin);
x=x(:);

beststd=inf;
bestlast=nan;
includeNext=false;
for last=max(args.basecycles):length(x)-1
  poly=polyfit((min(args.basecycles):last)',x(min(args.basecycles):last),1);
  if poly(1)>args.maxslope || poly(1)<args.minslope
    if poly(1)>args.maxslope
      xx=x-(1:length(x))'*args.maxslope;
    else
      xx=x-(1:length(x))'*args.minslope;
    end
    intcpt=mean(xx(min(args.basecycles):last));
    if args.debug
      fprintf('Initial poly for range [%d,%d] was %.2fx + %.2f, reducing slope to %.2fx + %.2f\n', min(args.basecycles),last,poly,args.maxslope,intcpt);
    end
    poly=[args.maxslope,intcpt];
  end
  baseline=polyval(poly,1:length(x))';
  delta=x-baseline;
  dstd=sqrt(mean(delta(min(args.basecycles):last).^2));
  if all(delta(last+1:min(last+10,end))>0)
    if args.debug
      fprintf('At last=%d, poly=[%.2f,%.2f], std=%.2f\n', last, poly,dstd);
    end
    if dstd<=max(beststd,args.maxstd)
      bestpoly=poly;
      bestbaseline=baseline;
      beststd=dstd;
      bestlast=last;
    end
  end
end
if isinf(beststd)
  bval=prctile(x,10);
  bestbaseline=ones(size(x))*bval;
  fprintf('findbase: Warning: Unable to identify slope of base, using %.0f as baseline\n',bval);
  bestpoly=[1,bval];
  bestlast=max(find(x<bestbaseline));
end
if bestpoly(1)>args.maxslope || bestpoly(1)<args.minslope
  fprintf('findbase: Warning: Baseline slope is %.2f; out of range [%.2f, %.2f]\n', bestpoly(1),args.minslope,args.maxslope);
end

if args.doplot
  setfig('findbase');clf;
  subplot(211);
  plot(x);
  hold on;
  rng=min(args.basecycles):min(bestlast+5,length(x));
  plot(rng,bestbaseline(rng),'r');
  plot(bestlast,bestbaseline(bestlast),'ro');
  c=axis;
  sel=1:min(bestlast+5,length(x));
  c(3)=min([bestbaseline(sel);x(sel)])-20;
  c(4)=max([bestbaseline(sel);x(sel)])+20;
  axis(c);
  subplot(212);
  xc=x-bestbaseline;
  semilogy(xc);
  hold on;
  semilogy(xc(1:bestlast),'r');
  c=axis;
  c(3)=1;
  c(4)=max(xc);
  axis(c);
end
