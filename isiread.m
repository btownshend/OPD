% Read an iCycler ISI file (exposure image)
function isi=isiread(filename)
fd=fopen(filename,'r');
x=fread(fd,'*uint16');
hdrlength=8;
trailerlength=1342+1050;
hdr=x(1:hdrlength);
trailer=x(end-trailerlength+1:end);
x=double(x(hdrlength+1:end-trailerlength));
r1=342;
r2=floor(length(x)/r1);
fprintf('len(x)=%d = %d x %d + %d\n', length(x),r1,r2,length(x)-r1*r2);
xim=x(1:r1*r2);
im=reshape(xim,r1,r2);
im=im/max(im(:));
setfig('isiread');clf;
imshow(im);
setfig('isi2');clf;
plot(x);
isi=struct('hdr',hdr,'trailer',trailer,'x',x,'im',im);

