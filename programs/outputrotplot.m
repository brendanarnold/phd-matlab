function outputrotplot(rotplotdata,pathfilename,params,append);
%write extremal orbit results to file

if exist(pathfilename,'file') & append==true
    ofile=fopen(pathfilename,'a');
else
    ofile=fopen(pathfilename,'w');
    fprintf(ofile,'min F=%s; num_slices=%s; num_pts=%s\r\n',num2str(params.minfreq),num2str(params.num_slices),num2str(params.num_pts));
    fprintf(ofile,'angle,B_x,B_y,B_z,band number,Wien2k band number,spin dirn,branch,area [Angstrom^-2],freq [T],k_x,k_y,k_z,k_parallel,slicenum,inplane dk (Angstrom^-1),dF/dkpara^2,mass (m_e),df/dM [T/(mu_B/f.u.)]\r\n');
end

for extremalorbit=1:size(rotplotdata,1)
    fprintf(ofile,'%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\r\n',rotplotdata(extremalorbit,:));
end
fclose(ofile);