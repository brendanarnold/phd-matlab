function results=calcrotplot(bandsdata,bandnums,FermiLevelOffsets,initBdirn,rotaxis,angleincORlist,numsteps,masscalc_deltaE,pathfilename,params,deltaM_bandsdata,deltaM);
%function results=calcrotplot(bandsdata,bandnums,FermiLevelOffsets,initBdirn,rotaxis,angleincORlist,numsteps,masscalc_deltaE,pathfilename,params,deltaM_bandsdata,deltaM);
%
%calls extremalorbits.m at a series of angles to make rotation plot
%arguments:
%bandsdata:         structure as output by convert_W2kE
%bandnums:          bands to include in rot plot calc (empty means all)
%FermiLevelOffsets: offset in Ryd to apply to Wien2k calc Fermi level (empty means zero). 2x1 matrix required for spin pol case
%initBdirn:         B dirn corresponding to zero angle in rot plot (Cartesian coord in reciprocal space)
%                       e.g. for hex this means [1 0 0], [0.5 sqrt(3)/2 0], [0 0 1] are parallel to reciprocal lattice vectors a, b and c respectively.
%rotaxis:           axis that B is rotated about in rot plot (Cartesian coord in reciprocal space)
%angleincORlist:    increment angle (if single num) or 1xn matrix of angles
%numsteps:          number of angles to calc at
%masscalc_deltaE:   E difference in Ryd used to calc band masses
%pathfilename:      output path and filename (empty means don't save!)
%params:            struct with fields that control orbit searches
%                       minfreq - smallest orbit that will be counted (def 100)
%                       num_pts - number of pts across kspace slice at kpara where slice is biggest (def 200)
%                       num_slices - number of slices at different kpara (def 100)
%deltaM_bandsdata:  bandstructure at close value of M (only used in spin pol calc for finding df/dM)
%deltaM:            difference in M between deltaM_bandsdata and bandsdata

%structure of output:
%col 1: angle
%2-4: Bdirn(x,y,z)
%5: bandnum
%6: Wien2k_bandnum
%7: spindirn
%8-17: same as extremalorbits

colnames={'angle','B_x','B_y','B_z','bandnum','W2k_bandnum','spindirn','branch','area','freq',...
    'kx','ky','kz','kpara','slicenum','dk','curvature','m_b','dFdM'};
if ~isempty(FermiLevelOffsets)
    bandsdata.FermiLevel=bandsdata.FermiLevel+FermiLevelOffsets;
end
if isempty(bandnums)
    bandnums=1:length(bandsdata.Wien2k_bandnums);
end
if isempty(params)
    load('defsliceparams.mat');
    disp(['Using default parameters: largest slice=' num2str(params.num_pts) 'x' num2str(params.num_pts) ...
        ', no. of slices per B dirn=' num2str(params.num_slices) ', min_F=' num2str(params.minfreq)]);
end

% Hacked slightly so that if cell array is passed, can have just one angle
% specified
if iscell(angleincORlist)
    angles = [angleincORlist{1}];
elseif length(angleincORlist)>1
    angles=angleincORlist;
    numsteps=length(angles);
else
    angles=(0:numsteps-1)*angleincORlist;
end

for bandnum=1:length(bandnums)
    Bdirn=initBdirn;
    disp(['Searching band ' num2str(bandnum) '/' num2str(length(bandnums)) ' (' num2str(bandsdata.Wien2k_bandnums(bandnums(bandnum))) ').']);
    rotplotdata{bandnum}=[];
    for stepnum=1:numsteps
        disp(['    Angle ' num2str(angles(stepnum))]);
        [Bdirn(1) Bdirn(2) Bdirn(3)]=rotate3axis(initBdirn(1),initBdirn(2),initBdirn(3),rotaxis,angles(stepnum)/360*2*pi);    
        eos=extremalorbits(Bdirn,bandsdata,bandnums(bandnum),params);
        numorbits=size(eos{1},1);
        if numorbits>0
            rotplotdata{bandnum}=[rotplotdata{bandnum}; [ones(numorbits,1)*[angles(stepnum) Bdirn bandnums(bandnum) bandsdata.Wien2k_bandnums(bandnums(bandnum)) bandsdata.spindirns(bandnums(bandnum))] eos{1}(:,:)]];
        end
        disp(['    Dirn ' num2str(stepnum) '/' num2str(numsteps) '. Num orbits:' num2str(numorbits)]);        
        pause(0.1); drawnow;
    end

    %calculate masses
    disp('Calculating masses (and df/dMs if deltaM_bandsdata supplied) for all extremal orbits.'); 
    if ~isempty(masscalc_deltaE) & abs(masscalc_deltaE)>0
        [masses dfdMs band_vertices{bandnum,1}]=calcmasses(rotplotdata{bandnum},bandsdata,masscalc_deltaE,deltaM_bandsdata,deltaM);
    else
        masses=zeros(size(rotplotdata{bandnum},1),1);
    end

    rotplotdata{bandnum}=[rotplotdata{bandnum} masses dfdMs];
    
    %assign orbits to branches according to continuity with angle
    %(replace existing 'branch' column)
    kdiff_threshold=sqrt(bandsdata.dx^2+bandsdata.dy^2+bandsdata.dz^2);
    fdiff_threshold=10000;
    rotplotdata{bandnum}(:,8)=grouporbits_test(rotplotdata{bandnum},angles,1,10,fdiff_threshold,11,kdiff_threshold,2);

    if size(rotplotdata{bandnum},1)==0
        rotplotdata{bandnum}=zeros(0,19);
    end
    if ~isempty(pathfilename)
        outputrotplot(rotplotdata{bandnum},[pathfilename '.orbits.txt'],params,bandnum~=1);
        %define rplotdata containing data for all bands so far
        rplotdata=cat(1,rotplotdata{:});
        vertices=cat(1,band_vertices{:});
        [pathstr,name,ext,versn] = fileparts(pathfilename);
        if ~isempty(pathstr)
            pathstr=['\' pathstr];
        end
        save([pathstr name '.orbits.mat'],'rplotdata','vertices','params','colnames');
    end
end

results=rotplotdata;