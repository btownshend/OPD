% Compare miner fit to FU data
function minerfit(opd,miner)
for j=1:length(miner.SampleNames)
  a=miner.Logistic_a(j);
  b=miner.Logistic_b(j);
  X0=miner.Logistic_X0(j);
  Y0=miner.Logistic_Y0(j);
  ct=miner.CT(j);
  cycle=opd.avg.cycle;
  mdl=Y0+a./(1+(cycle/X0).^b);
  plot(cycle,opd.avg.scaled(:,:,j));
  hold on;
  plot(cycle,mdl);
  xlabel('Cycle');
  ylabel('FU');
  fprintf('a=%.2f, b=%.2f, X0=%.2f, Y0=%.2f, Ct=%.2f\n', a,b,X0,Y0,ct);
  break;
end
