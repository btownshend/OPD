% Display well factor information
function wf=opdwellfactors(opd)
w=wellnames(opd);
wf=nan(8,12);
for j='A':'H'
  for i=1:12
    ind=find(strcmp(w,sprintf('%c%d',j,i)));
    if ~isempty(ind)
      wf(j-'A'+1,i)=opd.WFFP.val(ind);
    end
  end
end

c=wf;
c(9,13)=nan;
pcolor(c);
axis ij;
title('Well Factors');
set(gca,'XTick',1.5:12.5);
set(gca,'XTickLabel',{'1','2','3','4','5','6','7','8','9','10','11','12'});
set(gca,'YTick',1.5:8.5);
set(gca,'YTickLabel',{'A','B','C','D','E','F','G','H'});
for i=1:8
  for j=1:12
    text(j+0.15,i+0.5,sprintf('%.2f',c(i,j)));
  end
end

colorbar