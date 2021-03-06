function results=extremalorbits(normal,bandsdata,bandnums,params);
%
%function results=extremalorbits(normal,bandsdata,bandnums);
%calls orbitbranches.m to get kpara dependence of orbits
%establishes which of these have extremal points
%structure of output:
%col 1: branch (currently meaningless outside orbitbranches)
%cols 2-9: same as orbitbranches
%col 10: curvature factor d^2F/dkpara^2

%18/08/06 Note line 53, transpose operator removed as wrong for hex.

maxdist=max([bandsdata.dx bandsdata.dy bandsdata.dz]);
%maxdist outside bandsdata.centresvol to include extremal orbits

if isempty(bandnums)
    bandnums=1:length(bandsdata.Wien2k_bandnums);
end

optsobj = fitoptions('Method','NonlinearLeastSquares');
ftypeobj=fittype('f0+a*(k-k0).^2+b*(k-k0).^4','options',optsobj,'ind','k','coefficients',{'f0','k0','a','b'});

branches=orbitbranches_vol(normal,bandsdata,bandnums,[],params);
extremalorbits{1}=[];
for bandnum=1:length(bandnums)
for branch=1:length(branches{bandnum})
    %find indexes of extremal f points within branch
    numpoints=size(branches{bandnum}{branch},1);
    extremalindexes=[]; quadcoeff=[];
    curbranch=branches{bandnum}{branch};
    if numpoints>4
        for slice_ind=3:numpoints-2
            %slice_ind because it is index within branch i.e. slice_ind=1
            %is where branch starts not necessarily slicenum=1
            fdiffs=curbranch(slice_ind+2:-1:slice_ind-1,2)-curbranch(slice_ind+1:-1:slice_ind-2,2);
            if (fdiffs(2)*fdiffs(3))<=0 & (fdiffs(3)-fdiffs(2))*(fdiffs(2)-fdiffs(1))>=0 & ...
                    (fdiffs(3)-fdiffs(2))*(fdiffs(4)-fdiffs(3))>=0
                extremalindexes(end+1)=slice_ind;
                %fit cubic to 5 points around extremum
                ks=curbranch(slice_ind-2:slice_ind+2,6);
                poly=polyfit(ks,curbranch(slice_ind-2:slice_ind+2,2),3);
                k_ext=(-poly(2)+(-1:2:1)*sqrt(poly(2)^2-3*poly(1)*poly(3)))/(3*poly(1));
                %cubic might have >1 extremum so take the one nearest
                %centre
                k_ext=k_ext(abs(k_ext-mean(ks))==min(abs(k_ext-mean(ks))));
                if length(k_ext)~=1
                    disp(['Extremum fit gone wrong. slice_ind=' num2str(slice_ind) ' Branch=' num2str(branch)]);
                end
                f_ext=polyval(poly,k_ext);
                curv_fact=6*poly(1)*k_ext+2*poly(2);
                k_centre=[interp1(curbranch(:,6),curbranch(:,3),k_ext) interp1(curbranch(:,6),curbranch(:,4),k_ext) interp1(curbranch(:,6),curbranch(:,5),k_ext)];
                %only add to list of extremal orbits if centre is inside BZ
                %for n=1:size(bandsdata.centresvol,1)
                %    centre_kps(n)=norm(bandsdata.centresvol(n,:));
               %     centre_normals(n)=bandsdata.centresvol(n,:)/centre_kps(n);
               % end
               % outsideBZ=isoutside(centre_kps,centre_normals,bandsdata.ptinvol,k_centre,maxdist);
               % if ~outsideBZ
                    extremalorbits{bandnum}(end+1,:)=[branch 0 f_ext k_centre k_ext curbranch(slice_ind,7) curbranch(slice_ind,8) curv_fact];                    
               % end
            end
        end
    end
       
        %note output data includes branch assignment made by
        %orbitbranches.m which is only meaningful wrt different kpara for
        %same direction of normal
        
end
end

results=extremalorbits;
