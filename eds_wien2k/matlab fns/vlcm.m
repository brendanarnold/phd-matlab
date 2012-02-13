function vlcm=vlcm(n)
%function vlcm=vlcm(n)
%
%returns lowest common multiple of all numbers contained in n

plcm=n(1);
for i=1:length(n)-1
    vlcm=lcm(plcm,n(i+1));
    plcm=vlcm;
end