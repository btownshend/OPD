  % locate file
function file=opdlocate()
  dname='/Volumes/smolke-lab$/Brent';
  if ~exist(dname,'dir')
    error('Smolke server not mounted\n');
  end
  d=dir([dname,'/Data*']);
  if length(d)<1
    error('No data files found in %s\n',dname);
  end
  [~,ord]=sort(datenum({d.date}));
  d=d(ord);
  file=[dname,'/',d(end).name];
  fprintf('Tracking %s\n', file);
