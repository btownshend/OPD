% Print out info appropriate for pasting into Miner
function result=minerrun(v)
data='Cycle\t';
for i=1:length(v.WIRT)
  pos=v.WIRT(i).platepos;
  wellnames{i}=sprintf('%c%d',floor(double(pos)/12)+'A',mod(pos,12)+1);
  if i<length(v.WIRT)
    data=[data,sprintf('%s\t',wellnames{i})];
  else
    data=[data,sprintf('%s\n',wellnames{i})];
  end      
end

sc=v.avg.scaled;
for i=1:size(sc,1)
  data=[data,sprintf('%f\t',v.avg.cycle(i))];
  for j=1:size(sc,3)-1
    data=[data,sprintf('%f\t',sc(i,1,j))];
  end
  data=[data,sprintf('%f\n', sc(i,1,end))];
end
url1='http://ewindup.info/miner/data_submit.htm';
[r1,s1]=urlread2(url1);

% Execute miner
%url='http://www.htmlcodetutorial.com/cgi-bin/mycgi.pl';
url='http://ewindup.info/miner/miner_submit.asp';
args={'MaxPvalue','Default 0.05',
      'MinEfficiency','Default 0.0 (0%)',
      'MaxEfficiency','Default 2.0 (200%)',
      'UserEmail','bst@stanford.edu',
      'PlatformSelect','BioRad_iCycler',
      'Platform','BioRad_iCycler',
      'ExampleDataByColumn','',
      'ExampleDataByRow','',
      'inputDirection','0',
      'outputDirection','0',
      'RawDataSubmit','Submit raw data and Run!',
      'RawData',data
     };
args=args';
args=args(:);
fprintf('Sending request to %s...',url)
tic;
refheader=http_createHeader('Referer',url1);
[querystring,header]=http_paramsToString(args,1);
[result,status]=urlread2(url,'POST',querystring,[header,refheader]);
elapsed=toc;
fprintf('got status %d (%s) after %.1f seconds\n',status.status.value,status.status.msg,elapsed);
keyboard


