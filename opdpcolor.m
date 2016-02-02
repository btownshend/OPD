% Draw a pcolor plot of the OPD Ct's
function opdpcolor(opd)
%clf;
cttmp=opd.ctgrid;
cttmp(9,:)=nan;
cttmp(:,13)=nan;
pcolor(cttmp);
axis ij;
title('Ct');
set(gca,'XTick',1.5:12.5);
set(gca,'XTickLabel',{'1','2','3','4','5','6','7','8','9','10','11','12'});
set(gca,'YTick',1.5:8.5);
set(gca,'YTickLabel',{'A','B','C','D','E','F','G','H'});
for i=1:8
  for j=1:12
    text(j+0.15,i+0.5,sprintf('%.1f',cttmp(i,j)));
  end
end
colorbar

