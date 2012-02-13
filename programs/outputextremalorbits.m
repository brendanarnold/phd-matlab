function outputextremalorbits(extremalorbits,pathfilename);
%write extremal orbit results to file

ofile=fopen(pathfilename,'w');
fprintf(ofile,'branch,area [Angstrom^-2],freq [T],k_x,k_y,k_z,k_parallel,slicenum,dF/dkpara^2\r\n');
for orbit=1:size(extremalorbits,1)
    fprintf(ofile,'%e,%e,%e,%e,%e,%e,%e,%e,%e\r\n',extremalorbits(orbit,:));
end
fclose(ofile);