% Plot ct's from miner on 
function minerplot(miner,opd)
setfig('minerplot');clf;
wells=wellnames(opd);
nx=ceil(sqrt(length(wells)));
ny=ceil(length(wells)/nx);
for i=1:length(wells)
  subplot(nx,ny,i);
  opdcheck(opd,{wells{i}},'firststage',true);
  legend off;
  c=axis;
  hold on;
  plot(miner.CT(i)*[1,1],c(3:4),':');
  title(sprintf('%s: Ct=%.1f',wells{i},miner.CT(i)));
end
