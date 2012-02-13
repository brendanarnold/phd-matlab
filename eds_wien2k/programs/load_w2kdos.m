function [E DOS]=load_W2kDOS(filename);

fid=fopen(filename,'r');
data=textscan(fid,'%f %f %*f %*f %*f %*f','headerLines',3);
E=data{1};
DOS=data{2};
