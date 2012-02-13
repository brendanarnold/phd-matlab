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
    if param_counter==31
        debug=true;
    end
    currsetrows=find(orbitdata(:,param_col)==param_values(param_counter));
    assigned=[]; matchingnum=[];
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
                [minima matchingnum]=min(abs(deltacentre.*deltafreq),[],1);
                assigned=matchingnum;
                %remove multiple matches to same orbit
                while length(matchingnum)~=length(unique(matchingnum))
                    for counter=1:length(matchingnum)
                        multimatchrows=find(matchingnum==matchingnum(counter));
                        if length(multimatchrows)>1
                            [minminima keepmatch]=min(minima(multimatchrows));
                            removematch=setdiff(multimatchrows,keepmatch);
                            minima(removematch)=[];
                            matchingnum(removematch)=[];
                            assigned(removematch)=[];
                            break;
                        end
                    end
                end                  
                %remove best matches that are not good enough
                %i.e. those for which either kcentre or freq difference is too large
                goodmatch=[]; goodmatch(1:length(assigned))=1;
                for orbitnuminset=1:length(assigned)
                    if ((deltacentre(assigned(orbitnuminset),orbitnuminset))>kdiff_threshold) %...
                        %| ((abs(deltafreq(orbitnuminset,matchingnum(orbitnuminset))))>fdiff_threshold)
                        goodmatch(orbitnuminset)=0;
                    end
                end
                matchingnum=matchingnum(goodmatch==1);
                assigned=assigned(goodmatch==1);
                minima=minima(goodmatch==1);            
        end
   
        %assigned contains indexes into currsetrows for successfully
        %matched orbits; matchingnum contains corresponding indexes into
        %prevsetrows for matched orbits
        if length(assigned)>length(prevsetrows) | length(matchingnum)~=length(unique(matchingnum))
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
        end
    end    
        
    %make branch assignments for successfully matched orbits
    branches(currsetrows(assigned))=branches(prevsetrows(matchingnum));
                    
    %assign orbits that do not have match in prev set to new branches
    unassigned=setdiff(1:length(currsetrows),assigned);
    if length(unassigned>0)
        %means that new branch(es) have appeared
        branches(currsetrows(unassigned))=max(branches)+(1:length(unassigned));
    end
end           

        