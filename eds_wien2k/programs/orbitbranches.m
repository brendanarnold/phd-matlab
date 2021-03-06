function results=orbitbranches(normal,bandsdata,bandnums,pathfilename,params,kparas);
%makes slices at intervals along normal and finds orbit areas as function
%of kpara

%arguments:
%   normal:         direction of B field (Cartesian coords in reciprocal space)
%   bandsdata:      structure output by convert_W2kE containing bandstructure data
%   bandnums:       band numbers to do calc for (empty means all)
%   pathfilename:   file to save orbit branch info to (empty means don't save)
%   params:         structure containing fields that control orbit search (see calcrotplot.m)
%   kparas:         list of values of k_para_B to use (empty means function will calc kparas for params.num_slices in the 1st BZ interval)

%tries to establish which orbits belong to which continuous branches
%result=branches is a 1D cell array, each cell containing 2D array describing orbits at different
%kpara belonging to the same branch
%branches column format is:
%1-5: same as orbits
%6: kpara
%7: slicenum
%8: dk (point spacing within slice plane)
%18/08/06 Note line 35, transpose operator removed as wrong for hex.

normal=normal/norm(normal);

if isempty(bandnums)
    bandnums=1:length(bandsdata.Wien2k_bandnums);
end
%generate kpara series within first BZ
cart2reclatt=inv(bandsdata.rec_latt_vecs);
normalreclatt=normal*cart2reclatt;

if ~exist('kparas') || isempty(kparas)
    %kpara values have not been provided in input arguments
    %kpara range is found by calc intersection of normal with all planes that bound BZ
    bvecs=bandsdata.BZfacenormals*bandsdata.rec_latt_vecs;
    for planenum=1:size(bvecs,1)
        if dot(normal,bvecs(planenum,:))==0
            kparas(planenum)=NaN;
        else
            kparas(planenum)=norm(normal)*dot(bvecs(planenum,:),bvecs(planenum,:))/dot(normal,bvecs(planenum,:));
        end
    end
    kparamin=max(kparas(kparas<0));
    kparamax=min(kparas(kparas>0));
    dkpara=(kparamax-kparamin)/(params.num_slices-1);

    %max kcentre difference between orbits on consecutive slices for identifying as same
    %branch:
    kdiff_threshold=dkpara*5;

    %range of kpara needs to extend a bit beyond BZ so orbits close to or on
    %boundary can be tested for extremal
    kparas=kparamin-dkpara*4:dkpara:kparamax+dkpara*4;
else
    %kpara values have been provided in input arguments
    kdiff_threshold=min(diff(kparas))*5;
end
   
fdiff_threshold=0.3;
    
%getinplanedirn finds dirn in plane perp to normal and passing through
%gamma for which searchvol polygon is narrowest
[inplanedirns minwidth]=getinplanedirn(normal,bandsdata);
dk=minwidth/(params.num_pts-1);
%get orbit areas slice by slice
allorbits{1}=[];
for bandnum=1:length(bandnums)
    for slicenum=1:size(kparas,2)
        if slicenum==103
            debug=true;
        end
        slice=sliceFS(normal,inplanedirns,kparas(slicenum),dk,bandsdata,bandnums(bandnum));
        orbits=orbitareas(slice,bandsdata,bandsdata.FermiLevel(bandsdata.spindirns(bandnums(bandnum))),false,false);    
        %add to the list over all slices
        orbits{1}=[orbits{1} kparas(slicenum)*ones(size(orbits{1},1),1) slicenum*ones(size(orbits{1},1),1)];
        allorbits{bandnum}=[allorbits{bandnum}; orbits{1}];
    end

    %discard orbits with f outside limits
    allorbits{bandnum}=allorbits{bandnum}(allorbits{bandnum}(:,2)>=params.minfreq,:);


    if ~isempty(allorbits{bandnum})
        %identify orbits as fn of kpara 
        % FIXME: Seems to be a problem with the branch identification, just
        % allocate 1 to all branches
%         branchnumbers=grouporbits(allorbits{bandnum},kparas,6,2,fdiff_threshold,3,kdiff_threshold,1);
        branchnumbers = ones(size(allorbits{bandnum}, 1), 1);
        for branchnum=1:max(branchnumbers)
            branches{bandnum}{branchnum}(:,:)=allorbits{bandnum}(branchnumbers==branchnum,:);
            branches{bandnum}{branchnum}=[branches{bandnum}{branchnum} repmat(dk,[size(branches{bandnum}{branchnum},1) 1])];
        end
    else
        branches{1}=[];
    end
end

if ~isempty(pathfilename)
    outputbranches(branches,pathfilename);
end

results=branches;

