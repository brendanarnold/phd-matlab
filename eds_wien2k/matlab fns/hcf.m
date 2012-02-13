function ccf=hcf(d);

notatend=true; ccf=1;
while notatend
    cfs=factor(d(1));
    for n=1:numel(d);
        cfs=intersect(cfs,factor(d(n)));
    end
    if ~isempty(cfs)
        d=d./prod(cfs);
        ccf=ccf*prod(cfs);
    else
        notatend=false;
    end
end


    