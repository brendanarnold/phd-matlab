function masses=calcmasses(orbitdata,bandsdata,deltaE);
%function calcmasses(orbitdata,bandsdata,deltaE);
%
%calculates cyclotron masses of orbits in orbitdata
%requires columns: bandnum, Bdirn, kcentre(x,y,z), kpara, frequency

Ryd2Joules=13.6*1.6e-19;

if size(orbitdata,2)==16 |size(orbitdata,2)==17
    col_kpara=14;
    cols_kcentre=11:13;
    cols_Bdirn=2:4;
    col_bandnum=5;
    col_freq=10;
end

%max permissible change in kcentre of orbit as E_F -> E_F+deltaE
%(~data spacing)
kdiff_threshold=sqrt(bandsdata.dx^2+bandsdata.dy^2+bandsdata.dz^2)/2;

masses=zeros(size(orbitdata,1),1);
%recalculate F at E_F+deltaE to get mass of each extremal orbit
for orbitnum=1:size(orbitdata,1)
    kpara=orbitdata(orbitnum,col_kpara);
    freq=orbitdata(orbitnum,col_freq);
    slice=sliceFS(orbitdata(orbitnum,cols_Bdirn),kpara,bandsdata,orbitdata(orbitnum,col_bandnum));
    coplanar_orbits=orbitareas(slice,bandsdata,bandsdata.FermiLevel+deltaE);
    %find corresponding orbit 
    deltacentre=sqrt((orbitdata(orbitnum,cols_kcentre(1))-coplanar_orbits{1}(:,3)).^2+...
        (orbitdata(orbitnum,cols_kcentre(2))-coplanar_orbits{1}(:,4)).^2+...
        (orbitdata(orbitnum,cols_kcentre(3))-coplanar_orbits{1}(:,5)).^2);
    deltafreq=abs(freq-coplanar_orbits{1}(:,2));
    matching_kcentre=deltacentre<kdiff_threshold;
    if length(find(matching_kcentre))>1
        %another orbit may share kcentre so distinguish by value of F
        ind1=find(matching_kcentre);
        ind2=find(deltafreq(matching_kcentre)==min(deltafreq(matching_kcentre)));
        matching_index=ind1(ind2);
    elseif length(find(matching_kcentre))==1
        matching_index=find(matching_kcentre);
    else
        matching_index=[];
    end
    if length(matching_index)==1
        masses(orbitnum)=1.6e-19*((6.626e-34)/2/pi)/(9.1e-31)*(coplanar_orbits{1}(matching_index,2)-freq)/(deltaE*Ryd2Joules);
    else
        %indicate didn't find unambiguous matching orbit at E_F+deltaE
        masses(orbitnum)=999;
    end    
end
