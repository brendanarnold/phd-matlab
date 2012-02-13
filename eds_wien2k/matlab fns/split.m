function [up,dn]=split(data,col);
%function [up,dn]=split(data,col);

%splits single table of data into up,dn sweeps by finding max value of data in column col

maxind=find(data(:,col)==max(data(:,col)),1);
up=data(1:maxind,:);
dn=data(maxind+1:end,:);
