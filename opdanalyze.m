% opdanalyze - pull out data into regular structures:
%  v.all - all data points
%    structure with fields:
%      fluor(i,j,k,l) - raw fluorescence for timepoint i, filter j, reading k (1=in well, 2=outside well), sample l
%      netfluor(i,j,l) - rescaled by exposure time and pure dye caliburation, net fluor (in well-out well)
%      scaled(i,j,l) - further scaled by well factors
%      cycle(i) - fractional cycle number
%      time(i) - elapsed time from run start in minutes
%      clipped(i,j,l) - number of clipped pixels/sample
%      valid(i) - valid flag (unsure what this really is)
%  v.avg - average values per cycle
%      scaled(c,j,l) - fully scaled average value at cycle c, filter j, sample l
%      cycle(c) - integer cycle number

function v=opdanalyze(v)
% Get measurements from first rep only
%sel=[v.data.rep]==1;
if isempty(v.data)
  disp('No data to analyze.');
  return;
end

d=v.data;
all=struct();
% Need to figure out how to index this when there are multiple ones...
for i=1:length(v.RMEL)
  for j=1:length(v.RMEL(i).RMEF)
    % scale(i,j) - fluorophore i, using filter j
    pure(i,j)=v.RMEL(i).RMEF(j).response;   
  end
end
v.pure=pure;
scale=4*2500*inv(pure);  % Empirical; scales to normalize RMEF to 2500 (or sometimes 10000)
v.scale=scale;
currep=-1;
all.filters=unique([d.filter1]);
all.wscaled=[];
i=0;
lastfilter=-1;
for j=1:length(d)
  filter=find(d(j).filter1==all.filters);
  if filter~=lastfilter || length(all.filters)==1
    i=i+1;
  end
  lastfilter=filter;
  all.fluor(i,filter,:,:)=d(j).v;
  all.clipped(i,filter,:)=d(j).clipped;
  all.exposure(i,filter)=d(j).exposure;
  all.netfluor(i,filter,:)=(d(j).v(1,:)-d(j).v(2,:))*640/d(j).exposure;
  for k=1:size(all.netfluor,3)
    % Scale by well factor
    all.wscaled(i,filter,k)=all.netfluor(i,filter,k)/v.WFFP(filter).val(k);
  end
  all.valid(i,filter)=d(j).validFlag;
  if d(j).rep ~= currep
    cyclestart=d(j).time;
    currep=d(j).rep;
  end
  all.cycle(i)=double(d(j).rep-1)+(60*(d(j).time-cyclestart)+0.5)/double(d(j).dwellTime);
  all.time(i)=d(j).time;
  all.stage(i)=d(j).stage;
  all.temp(i)=d(j).temperature1;
end
all.netfluor(all.netfluor==0)=nan;
if 1
% Interpolate missing values of wscaled
for j=1:size(all.wscaled,2)
  missing=all.wscaled(:,j,1)==0;
  fmissing=find(missing);
  for k=1:length(fmissing)
    % Find closest non-missing sample
    [mt,pos]=min(abs(all.time(~missing)-all.time(fmissing(k))));
    funmiss=find(~missing);  
    repl=funmiss(pos);
    all.wscaled(fmissing(k),j,:)=all.wscaled(repl,j,:);
  end
end
end
for i=1:size(all.wscaled,3)
  all.scaled(:,:,i)=all.wscaled(:,:,i)*v.scale;
end
v.all=all;
v.stageavg={};
for stage=unique(all.stage)
  % Per cycle stats
  avg=struct();
  stagesel=all.stage==stage;
  if any(~stagesel)
    fprintf('Using only %d/%d samples that are for stage %d\n', sum(stagesel),length(stagesel),stage);
  end

  avg.cycle=unique(floor(all.cycle(stagesel)));
  for i=1:length(avg.cycle)
    cycle=avg.cycle(i);
    for f=1:size(all.scaled,2)
      sel=floor(all.cycle)==cycle & all.stage==stage;
      avg.scaled(i,f,:)=mean(all.scaled(sel,f,:),1);
      avg.temp(i)=mean(all.temp(sel));
    end
  end
  v.stageavg{end+1}=avg;
  if stage==all.stage(1) 
    % First stage, copy into all, avg for backward compat
    v.avg=avg;
  end
end
