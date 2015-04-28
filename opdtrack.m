% Follow OPD output in real-time
function opdtrack(file)
if nargin<1
  file=opdlocate();
end

setfig('opdtrack'); clf;
while true
  z=opdread(file);
  setfig('opdtrack');
  opdcheck(z);
  pause(30);
  fprintf('Updating...\n');
end