% Plot melt curves
function [s,st,sf]=opdmelt2(v,wells)
if nargin<2
  wells=1:size(v.all.scaled,3);
end
clf;
for i=1:length(v.WIRT)
  pos=v.WIRT(i).platepos;
  wellnames{i}=sprintf('%c%d',floor(double(pos)/12)+'A',mod(pos,12)+1);
end
% Select melt cycles (assume it is stage 7)
stage=v.PIFB.step([v.PIFB.step.repeat]==100).cycle;
sel=(v.all.stage==stage);
temp=[v.data(sel).temperature1];
trange=max(temp)-min(temp);
npoints=sum(sel);
fprintf('Taking melt data from cycle %d: npoints=%d, temp=%.1f-%.1f\n',stage,npoints,min(temp),max(temp));
setfig('response');clf;
setfig('melt');clf;
ut=unique(temp);
alltemp=min(ut):.05:max(ut);
for dye=1:size(v.all.scaled,2)
  s=squeeze(v.all.scaled(sel,dye,wells));
  st=[];
  for t=1:length(ut)
    st(t,:)=mean(s(ut(t)==temp,:),1);
  end
  sf=[];
  for k=1:size(s,2)
    pp=pchip(ut,st(:,k)');
    sf(:,k)=ppval(pp,alltemp);
    sfg(:,k)=gradient(sf(:,k));
  end
  subplot(size(v.all.scaled,2),2,dye*2-1);
  plot(ut,st,'r');
  hold on;
  plot(alltemp,sf,'g');
  subplot(size(v.all.scaled,2),2,dye*2);
  plot(alltemp,-sfg,'-');
  xlabel('Temperature');
  ylabel('delta Fluorescence');
  title(sprintf('%s - %s',v.filename,v.WFFP(dye).dye));
  if length(wellnames(wells))<20
    legend(wellnames{wells},'Location','NorthWest');
  end
end

