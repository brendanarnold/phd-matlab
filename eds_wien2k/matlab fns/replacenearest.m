function output=replacenearest(source,values,tol);
%function output=replacenearest(source,values,tol);
%
%replaces numbers in source with nearest values from values, writes a nan if none of values is
%within tol of source element

for n=1:numel(source)
    output(n)=values(abs(values-source(n))==min(abs(values-source(n))));
    if min(abs(values-source(n)))>tol
        output(n)=nan;
    end
end
output=reshape(output,size(source));