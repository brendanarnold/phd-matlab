function [norms dirs]=normdir(m);
%function normdir(m);
%
%m is an n*3 matrix containing plane normals
%returns magnitudes and dirns separately

for n=1:size(m,1)
    norms(n)=norm(m(n,:));
    dirs(n,:)=m(n,:)/norms(n);
end