function startmfit;

mfit;
inputpath='c:\work\ZrZn2\YellowMag\fieldmod\Au_calib\rawdata\';
set(findobj('Tag','mf_DataDir'),'string',inputpath);

filenames={'ey280505e'};

initialpars=[...
    0.0147 ...      %amp
    -3.7395e-5 ...   %amp_
    -5.2522e-4 ...   %amp__
    1.5444e3 ...    %freq
    -16.4687 ...    %phase
    -0.0043 ...     %backgd
    3.8278e-4 ...   %backgd_
    2.222e-7 ...  %backgd__
    -0.0085 ...     %fieldshift
    0.0686 ...     %harm2frac
    0.06 ...     %harm2phaseoffset
    0];     %harm2phaseoffset

Bstart=8.7; Bend=9.3;

mf_batch('FitFuncFile LKsingleF2');
argstr=['load ' filenames{1} '.dat'];
mf_batch(argstr);


for parnum=1:length(initialpars)
    mf_batch(['par ' num2str(parnum) ' ' num2str(initialpars(parnum))]);
end

end
function opt_T;
    mf_batch('fix all');
    mf_batch('free 5,14');
    mf_batch('fit');
end
function opt_backgd;
    mf_batch('fix all');    
    mf_batch('free 9,10');
    mf_batch('fit');
end
function opt_phases;
    mf_batch('fix all');
    mf_batch('free 4,8');
    mf_batch('fit');
end

function saveresultstofile(results);
outfile=fopen('c:\work\ZrZn2\yellowmag\fieldmod\phaseres1_8p7to9p3_downsweeps_101006.txt','w');
fprintf(outfile,'filename,T(K),phase1,err(phase1),phase2,err(phase2),chisqd,delta(phase1-phase2),err(delta phase)\r\n');
results.deltadeltaphase=(results.fitvals(:,4)-results.fitvals(:,8));
results.deltadeltaphase=results.deltadeltaphase(:)-results.deltadeltaphase(1);
results.errddphase=sqrt(results.fiterrs(:,4).^2+results.fiterrs(:,8).^2);
for n=1:length(results.stats.r2)
    
    
    fprintf(outfile,'%s,%f,%f,%f,%f,%f,%f,%f,%f\r\n',results.filenames{n},abs(results.fitvals(n,14)),results.fitvals(n,4),results.fiterrs(n,4),...
        results.fitvals(n,8),results.fiterrs(n,8),results.stats.chisqd(n,1),results.deltadeltaphase(n),...
        results.errddphase(n));    
end
fclose(outfile);
end