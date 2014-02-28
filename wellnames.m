function w=wellnames(v)
if ~isfield(v,'WIRT')
  w={};
  return;
end
for i=1:length(v.WIRT)
  pos=v.WIRT(i).platepos;
  w{i}=sprintf('%c%d',floor(double(pos)/12)+'A',mod(pos,12)+1);
end
