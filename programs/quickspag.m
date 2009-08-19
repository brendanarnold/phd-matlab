function quickspag(bandsdata,bandnums,special_kpts,coord_sys,num_kpts,labels)
%function quickspag(bandsdata,bandnums,special_kpts,coord_sys,num_kpts,labels)
%
%makes spaghetti plot of bands bandnums on path defined by special_kpts. 
%coordsys='C' (cartesian) or 'R' rec latt
%total no. of pts given by num_kpts, labels is cell array of special kpt labels
%quick way of making spag plot but better to run lapw1 with a proper .klist for the path

dirnstr={'up','dn'};
Epath{1}=[];
if isempty(bandnums)
    bandnums=1:length(bandsdata.Wien2k_bandnums);
end

%k route in rec_latt_vec basis:
kroute=special_kpts;
%convert from rec_latt_vec coords to Cartesian:
if strcmp(coord_sys(1),'R')
    kroute=kroute*bandsdata.rec_latt_vecs;
end

%make the kpath between the given symm points
[kpath givenpts]=make_kpath(kroute,num_kpts);

%wrap kpath points into cuboid and interpolate to get E
wrappedkpath=wrapcoords(kpath(:,1),kpath(:,2),kpath(:,3),bandsdata);
for bandnum=1:length(bandnums)
    Epath{bandnum}=interp3(permute(bandsdata.cartX,[2 1 3]),permute(bandsdata.cartY,[2 1 3]),permute(bandsdata.cartZ,[2 1 3]),permute(bandsdata.cartE{bandnums(bandnum)},[2 1 3]),...
        wrappedkpath.X,wrappedkpath.Y,wrappedkpath.Z);
    Epath{bandnum}=13.61*(Epath{bandnum}-bandsdata.FermiLevel(bandsdata.spindirns(bandnums(bandnum))));
end

figure;
hold on;
for bandnum=1:length(bandnums)
    colchars='rgbkymcw';
    colchar=colchars(mod(bandsdata.Wien2k_bandnums(bandnums(bandnum)),7)+1);
    switch bandsdata.spindirns(bandnums(bandnum))
        case 1
            style=['-' colchar];
        case 2
            style=[':' colchar];
    end           
    hplot(bandnum)=plot(Epath{bandnum},style);
    legendstr{bandnum}=['Band ' num2str(bandsdata.Wien2k_bandnums(bandnums(bandnum))) ', ' dirnstr{bandsdata.spindirns(bandnums(bandnum))}];   
end

legend(hplot,legendstr);

ylim([-2 1]);
currylim=ylim;
for n=1:length(givenpts)
    line([givenpts(n); givenpts(n)],[currylim(1); currylim(2)],'LineWidth',2,'Color','k');
    text(givenpts(n),currylim(1)-(currylim(2)-currylim(1))/10,labels(n));
end
line(xlim,[0 0],'LineWidth',1,'Color','k');

end

function [kpath givenpts]=make_kpath(kpointlist,numpoints);
kpath=[]; givenpts=[];
for segment=1:size(kpointlist,1)-1
    seglength(segment)=sqrt((kpointlist(segment+1,1)-kpointlist(segment,1))^2 +...
        (kpointlist(segment+1,2)-kpointlist(segment,2))^2+...
        (kpointlist(segment+1,3)-kpointlist(segment,3))^2);
end

totlength=sum(seglength);
for segment=1:size(kpointlist,1)-1
    numsegpts=round(seglength(segment)/totlength*numpoints);
    currnumpathpts=size(kpath,1);
    givenpts(end+1)=currnumpathpts+1;
    kpath(currnumpathpts+1:currnumpathpts+numsegpts,1)=linspace(kpointlist(segment,1),kpointlist(segment+1,1),numsegpts);
    kpath(currnumpathpts+1:currnumpathpts+numsegpts,2)=linspace(kpointlist(segment,2),kpointlist(segment+1,2),numsegpts);
    kpath(currnumpathpts+1:currnumpathpts+numsegpts,3)=linspace(kpointlist(segment,3),kpointlist(segment+1,3),numsegpts);
end
givenpts(end+1)=size(kpath,1);
end
