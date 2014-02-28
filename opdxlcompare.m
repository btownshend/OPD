% Assume xl contains all points data read from excel file
% opdcheck has been run on same dataset to get all.*
% Compare xl(:,2) (first well) with all.scaled(1,:,:)
a=xl(:,2);
b=all.scaled(:,1,1);
match=zeros(length(a),length(b),'uint8');
maxdiff=0.005;
vbest=nan*zeros(length(a),1);
for i=1:length(a)
  match(i,:)=abs(a(i)-b)<maxdiff;
  nm(i)=sum(match(i,:));
  if nm(i)==1
    vbest(i)=find(match(i,:));
  end
end
figure(gcf);
clf;
subplot(211);
plot(vbest);
xlabel('Excel data points');
ylabel('OPD data points');
subplot(212);
xlpos=(1:length(vbest))';
plot(xlpos,xlpos-vbest);
ylabel('Excel Position - OPD Position');
xlabel('Excel position');




