function value=propval(array,name);
%function value=propval(array,name);
%
%returns property value corresponding to name in cell array of property name/value pairs

value=[];
for n=1:2:numel(array)
    if isstr(array{n}) & strcmp(name,array{n})
        value=array{n+1};
        break;
    end
end

