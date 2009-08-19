function results=orbitbranches(normal,bandsdata,bandnums,pathfilename);
%makes slices at intervals along normal and finds orbit areas as function
%of kpara
%tries to establish which orbits belong to which continuous branches
%result=branches is a 1D cell array, each cell containing 2D array describing orbits at different
%kpara belonging to the same branch
%branches column format is:
%1-5: same as orbits
%6: kpara
%7: slicenum

minfreq=100; %tesla
num_slices=100;

if isempty(bandnums)
    bandnums=1:length(bandsdata.Wien2k_bandnums);
end
%generate kpara series within first BZ
cart2reclatt=inv(bandsdata.rec_latt_vecs);
normalreclatt=normal*cart2reclatt;
kparamin=-sqrt(dot(normal,normal))*(0.75/max(abs(normalreclatt)));
kparamax=sqrt(dot(normal,normal))*(0.75/max(abs(normalreclatt)));
dkpara=(kparamax-kparamin)/(num_slices-1);
%max kcentre difference between orbits on consecutive slices for identifying as same
%branch:
kdiff_threshold=dkpara*5;
fdiff_threshold=0.3;

kparas=kparamin:dkpara:kparamax;

%get orbit areas slice by slice
allorbits{1}=[];
for bandnum=1:length(bandnums)
for slicenum=1:size(kparas,2)
    slice=sliceFS(normal,kparas(slicenum),bandsdata,bandnums(bandnum));
    orbits=orbitareas(slice,bandsdata,bandsdata.FermiLevel);    
    %add to the list over all slices
    orbits{1}=[orbits{1} kparas(slicenum)*ones(size(orbits{1},1),1) slicenum*ones(size(orbits{1},1),1)];
    allorbits{bandnum}=[allorbits{bandnum}; orbits{1}];
end

%discard orbits with f outside limits
allorbits{bandnum}=allorbits{bandnum}(allorbits{bandnum}(:,2)>=minfreq,:);

%identify orbits as fn of kpara and calculate derivative

%step through slices, assigning orbits to branch on basis of smallest
%change in orbit centre position and frequency
numbranches=0;
prevsliceorbits=[]; branchassignments=[]; branches{bandnum}=[];
for slicenum=1:num_slices
    if slicenum>1
        prevsliceorbits=sliceorbits;
        prevbranchassignments=branchassignments;
    end
    sliceorbits=allorbits{bandnum}(allorbits{bandnum}(:,7)==slicenum,:);
    %generate nxm array of centre-centre distances where n is numorbits of
    %current slice and m is numorbits of prev slice
    if slicenum==16
        a=1;
    end
    if size(prevsliceorbits,1)>0 & size(sliceorbits,1)>0
        %if this is not the first slice with orbits, must do branch identifying
        deltacentre=[]; deltafreq=[];
        for orbitnum=1:size(sliceorbits,1)
            deltacentre(orbitnum,:)=sqrt((sliceorbits(orbitnum,3)-prevsliceorbits(:,3)).^2+...
                (sliceorbits(orbitnum,4)-prevsliceorbits(:,4)).^2+...
                (sliceorbits(orbitnum,5)-prevsliceorbits(:,5)).^2);
            deltafreq(orbitnum,:)=(sliceorbits(orbitnum,2)-prevsliceorbits(:,2))/sliceorbits(orbitnum,2);
        end
        %find best matches from matrix of kcentre differences
        [minima matchingnum]=min(abs(deltacentre.*deltafreq),[],2);
        assigned=1:length(matchingnum);
        %remove best matches that are not good enough
        %i.e. those for which either kcentre or freq difference is too large
        goodmatch=[]; goodmatch(1:orbitnum)=1;
        for orbitnum=1:size(sliceorbits,1)
            if ((deltacentre(orbitnum,matchingnum(orbitnum)))>kdiff_threshold) ...
                | ((abs(deltafreq(orbitnum,matchingnum(orbitnum))))>fdiff_threshold)
                goodmatch(orbitnum)=0;
            end
        end
        matchingnum=matchingnum(goodmatch==1);
        assigned=assigned(goodmatch==1);
        minima=minima(goodmatch==1);            

        %assigned contains indexes into sliceorbits for successfully
        %matched orbits; matchingnum contains corresponding indexes into
        %prevsliceorbits for matched orbits
        %branchassignments has length size(sliceorbits,1) and contains the
        %branch numbers that the orbits are assigned to
        
        if length(assigned)>size(prevsliceorbits,1) | length(matchingnum)~=length(unique(matchingnum))
            %more good matches than orbits on prev slice OR >1 good match
            %to particular orbit
            ['Ambiguity identifying branches.']
            for matchnum=1:length(assigned)
                %match is considered ambiguous if another orbit in current slice is matched
                %to the same orbit in previous slice
                ambiguous=length(find(matchingnum==matchingnum(matchnum)))>1;
                if ambiguous
                    %get the indexes into assigned and matchingnum of the
                    %ambiguously assigned current slice orbits
                    possmatchingnums=find(matchingnum(matchnum));
                    
                    %discard the worse of the good matches (better than
                    %random)
                end
            end
        else
            %first make branch assignments for successfully matched orbits
            for matchedorbitnum=1:length(assigned)
                branchassignments(assigned(matchedorbitnum))=prevbranchassignments(matchingnum(matchedorbitnum));
                branches{bandnum}{branchassignments(assigned(matchedorbitnum))}(end+1,:)=sliceorbits(assigned(matchedorbitnum),:);
            end

            %assign orbits that do not have match in prev slice to new
            %branches
            unassigned=setdiff(1:size(sliceorbits,1),assigned);
            if length(unassigned>0)
                %means that new branch(es) have appeared
                for newbranchorbit_num=1:length(unassigned)
                    newbranch=1+length(branches{bandnum});
                    branches{bandnum}{newbranch}(1,:)=sliceorbits(unassigned(newbranchorbit_num),:);
                    branchassignments(unassigned(newbranchorbit_num))=newbranch;
                end
            end
            if size(sliceorbits,1)<size(prevsliceorbits,1)
                %means existing branches have disappeared        
            end
        end
    elseif size(sliceorbits,1)>0 
        %'initialise' branches with orbits on a nonempty slice that
        %immediately follows an empty slice. i.e. either first nonempty
        %slice or nonempty slice after an empty slice
        for orbitnum=1:size(sliceorbits,1)
            branchassignments(orbitnum)=length(branches{bandnum})+1;
            branches{bandnum}{branchassignments(orbitnum)}(1,:)=sliceorbits(orbitnum,:);
        end
    end    
end
end

if nargin>3
    outputbranches(branches,pathfilename);
end

results=branches;
