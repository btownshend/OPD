% opdprotocol - show protocol used in OPD file
function opdprotocol(opd)
steps=opd.PIFB.step;
for i=1:length(steps)
  s=steps(i);
  if s.duration==32767
    dur='HOLD';
  else
    dur=sprintf('%2d:%02d',floor(single(s.duration)/60.0),mod(s.duration,60));
  end
  if s.gradient>0.1
    if s.repeat==32768
      temp=sprintf('%.1f-%.1f',s.temperature+[0,s.gradient]);
    else
      temp=sprintf('%.1f:%.1f:%.1f',s.temperature,s.gradient,s.temperature+s.gradient*(single(s.repeat)));
    end
  else
    temp=sprintf('%.1f',s.temperature);
  end
  if s.step==1
    fprintf('Cycle %2d.1, duration=%s, temp=%s, rpt=%d\n', s.cycle,dur,temp,s.repeat);
  else
    fprintf('        .%d, duration=%s, temp=%s\n', s.step,dur,temp);
  end
end
