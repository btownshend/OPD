% Follow OPD output in real-time
function opdtrack(file)
if nargin<1
  file=opdlocate();
end

while true
  z=opdread(file);
  opdcheck(z);
  figure(gcf);
  pause(30);
  fprintf('Updating...\n');
end