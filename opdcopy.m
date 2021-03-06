% Copy latest OPD file to local directory
function opdcopy(force)
fn=opdlocate();
slash=find(fn=='/',1,'last');
basename=fn(slash+1:end);
if exist(basename,'file')
  if nargin<1 || ~force
    error('%s already exists in current dir.',basename);
  end
end
cmd=sprintf('cp "%s" .',fn);
[s,r]=system(cmd);
if r~=0
  error('Failed execution of: %s\n%s\n', cmd, s);
end
  