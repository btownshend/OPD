% Compare miner fit to FU data
function minerfit(opd,miner,j)
a=miner.Logistic_a(j);
b=miner.Logistic_b(j);
X0=miner.Logistic_X0(j);
Y0=miner.Logistic_Y0(j);
ct=miner.CT(j);
cycle=opd.avg.cycle+0.63;
mdl=Y0+a./(1+(cycle/X0).^b);
wnames=wellnames(opd);
ti=sprintf('Minerfit %s',wnames{j})
setfig(ti);clf;
plot(cycle,opd.avg.scaled(:,:,j));
hold on;
plot(cycle,mdl);
c=axis;
c(3:4)=[0,max(opd.avg.scaled(:))];
axis(c);
plot(ct*[1,1],c(3:4),':');
plot(c(1:2),miner.DynamicThreshold(j)*[1,1],':');
fprintf('a=%.2f, b=%.2f, X0=%.2f, Y0=%.2f, Ct=%.2f\n', a,b,X0,Y0,ct);
xlabel('Cycle');
ylabel('FU');
title(ti);
