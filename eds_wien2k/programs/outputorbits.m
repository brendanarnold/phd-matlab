function outputorbits(orbits,pathfilename);
%write orbits results array to file

ofile=fopen(pathfilename,'w');
fprintf(ofile,'area [Angstrom^-2],freq [T],k_x,k_y,k_z,k_parallel,slicenum\r\n');
for row=1:size(orbits,1)
    fprintf(ofile,'%e,%e,%e,%e,%e,%e,%e\r\n',orbits(row,:));
end
fclose(ofile);