% Get a grid of wellnames
function w=wellgrid()
w={};
for j=1:8
  for i=1:12
    w{j,i}=sprintf('%c%d',j+'A'-1,i);
  end
end
end

