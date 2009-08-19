function branches=grouporbits(orbitdata,param_values,param_col,freq_col,fdiff_threshold,centre_col,kdiff_threshold,criterion);
%function branches=grouporbits(orbitdata,param_values,param_col,freq_col,fdiff_threshold,centre_col,kdiff_threshold,criterion);
%
%identifies orbit branches that are continuous wrt changes in param_values
%criterion: 1 means we are identifying branches wrt variations in kpara
%2 means wrt different field directions (for rot plot)
%3 means wrt different energies

numbranches=0; prevsetrows=[];
branches=zeros(size(orbitdata,1),1);

%cycle through distinct parameter values
for param_counter=1:length(param_values)
    if param_counter>1
        prevsetrows=currsetrows;
    end
    if param_counter==67
        debug=true;
    end
    currsetrows=find(orbitdata(:,param_col)==param_values(param_counter));
    matchingnum=[]; goodmatch=[];
    if length(prevsetrows)>0 & length(currsetrows)>0
        %if this is not the first slice with orbits, must do branch
        %identifying
        deltacentre=[]; deltafreq=[];
        %generate nxm array of centre-centre distances, freq differences
        %where n is numorbits in current parameter value set, m is
        %numorbits of prev set
        for orbitnuminset=1:length(currsetrows)
            deltacentre(orbitnuminset,:)=sqrt(...
                (orbitdata(currsetrows(orbitnuminset),centre_col)-orbitdata(prevsetrows,centre_col)).^2+...
                (orbitdata(currsetrows(orbitnuminset),centre_col+1)-orbitdata(prevsetrows,centre_col+1)).^2+...
                (orbitdata(currsetrows(orbitnuminset),centre_col+2)-orbitdata(prevsetrows,centre_col+2)).^2);
                
            deltafreq(orbitnuminset,:)=(orbitdata(currsetrows(orbitnuminset),freq_col)-orbitdata(prevsetrows,freq_col))/orbitdata(currsetrows(orbitnuminset),freq_col);
        end
        switch criterion
            case {1,2,3}
                %find best matches from matrix of kcentre differences                
                normedfs=repmat(max(abs(deltafreq),[],1),[size(deltacentre,1) 1]);
                proxmatrix=deltacentre.*(1-abs(deltafreq)./normedfs.*kdiff_threshold./deltacentre);
                [minima matchingnum]=min(proxmatrix,[],1);
                goodmatch=true(numel(prevsetrows),1);
                %numel(matchingnum) equals numel(prevsetrows);
                %matchingnum(index1) contains index of closest orbit in currsetrows to orbit index1
                %in prevsetrows

                %remove multiple matches to same orbit
                counter2=0;
                while numel(matchingnum(goodmatch))~=numel(unique(matchingnum(goodmatch)))                   
                    for counter=1:numel(matchingnum)
                        if numel(matchingnum==matchingnum(counter))>1
                            multimatchrows=find(matchingnum==matchingnum(counter));
                            [minminima keepmatch]=min(minima(multimatchrows));
                            removematch=setdiff(multimatchrows,multimatchrows(keepmatch));
                            goodmatch(removematch)=false;
                        end
                    end
                    counter2=counter2+1;
                    if counter2>100
                        debug=true;
                    end
                end                  
                %remove best matches that are not good enough
                %i.e. those for which either kcentre or freq difference is too large
                goodmatches=matchingnum(goodmatch);
                prevsetindex=find(goodmatch);
                for counter=1:numel(goodmatch)
                    if goodmatch(counter) & ((deltacentre(matchingnum(counter),counter))>kdiff_threshold) %...
                        %| ((abs(deltafreq(orbitnuminset,matchingnum(orbitnuminset))))>fdiff_threshold)
                        goodmatch(counter)=false;
                    end
                end
        end
    %make branch assignments for successfully matched orbits
    branches(currsetrows(matchingnum(goodmatch)))=branches(prevsetrows(goodmatch));  
    end    
                            
    %assign orbits that do not have match in prev set to new branches
    unassigned=setdiff(1:numel(currsetrows),matchingnum(goodmatch));
    if length(unassigned>0)
        %means that new branch(es) have appeared
        branches(currsetrows(unassigned))=max(branches)+(1:length(unassigned));
    end
end           

        