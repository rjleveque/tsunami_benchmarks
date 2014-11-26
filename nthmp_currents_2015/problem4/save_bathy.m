% Read in the Matlab bathy file, select the onshore portion only, 
% and then save as ascii files for use in Python.

clf
clear all

fid(1)=fopen('../all_data/4_Seaside_OSU_Model_Lab/bathy/seaside_matrix.bin', 'r');
dum=fread(fid(1),1,'int32');
m=fread(fid(1),1,'int32');
n=fread(fid(1),1,'int32');
x=fread(fid(1),m,'float32');
y=fread(fid(1),n,'float32');
z=fread(fid(1),[m,n],'float32');
fclose all;

dx=x(2)-x(1);
xmax=x(m);
xmin=x(1);

I = find(x>=32.5);   % onshore portion
x_onshore = x(I);
z_onshore = z(I,:);

save('z_onshore.txt','z_onshore','-ascii')
save('x_onshore.txt','x_onshore','-ascii')
save('y_onshore.txt','y','-ascii')
