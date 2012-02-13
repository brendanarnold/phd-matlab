function vlcm=LCM(n)

plcm=n(1)
for i=1:length(n)-1
    vlcm=lcm(plcm,n(i+1));
end