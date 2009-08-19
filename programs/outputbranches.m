function outputbranches(branches,pathfilename);
%write kpara dependent orbits results array to file organised into branches

ofile=fopen(pathfilename,'w');
fprintf(ofile,'bandnum,branch,area [Angstrom^-2],freq [T],k_x,k_y,k_z,k_parallel,slicenum\r\n');
for bandnum=1:numel(branches);
for branch=1:numel(branches{bandnum})
for orbit=1:size(branches{bandnum}{branch},1)
    fprintf(ofile,'%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\r\n',bandnum,branch,branches{bandnum}{branch}(orbit,:));
end
end
end
fclose(ofile);