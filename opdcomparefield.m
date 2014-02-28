% opdcomparefield - compare each given structure's value for given
% field
function opdcomparefield(s,f)
maxmatches=5;
if nargin>1
  fexpand=strrep(f,'.',''',''');
  cmd=sprintf('getfield(s{i},''%s'')',fexpand);
  for i=1:length(s)
    if ~isempty(s{i}) && isfield(s{i},fexpand)
      x{i}=eval(cmd);
    end
  end
else
  x=s;
  f='';
end
if iscell(x{1})
  disp('Cell array - comparing first elements only');
  for i=1:length(x)
    x{i}=x{i}{1};
  end
end
if isstruct(x{1})
  if length(x{1})>1
    disp('Struct array - comparing first element only');
  end
  % Recursive structure
  fprintf('\n<%s\n',f);
  fields=fieldnames(x{1});
  for i=1:length(fields)
    opdcomparefield(x,char(fields(i)));
  end
  fprintf('/%s>\n',f);
  return;
end
done=zeros(size(x));
for i=1:length(x)
  if done(i) || isempty(x{i})
    continue;
  end
  done(i)=1;
  matches=i;
  for j=i+1:length(x)
    if done(j)
      continue;
    end
    if isequal(x{i},x{j})
      matches(end+1)=j;
      done(j)=1;
    end
  end
  if length(matches)==length(x)
    % Suppress printing if all match
    break;
  end
  maxlen=20;
  if isinteger(x{i})
      fmt='%d ';
  elseif isfloat(x{i})
    fmt='%f ';
  else
    fmt='%s ';
    maxlen=200;
  end
  fprintf('%s: ',f);
  if (length(x{i})>maxlen)
    fprintf(fmt,x{i}(1:maxlen));
    fprintf('... ');
  else
    fprintf(fmt,x{i});
  end
  fprintf('(%dx) [',length(matches));
  fprintf('%d ',matches(1:min(length(matches),maxmatches)));
  if length(matches)>maxmatches
    fprintf('...');
  end
  fprintf(']\n');
end
