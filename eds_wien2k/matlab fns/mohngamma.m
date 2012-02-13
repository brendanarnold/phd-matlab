function result=MohnGamma(t,stonerishness,outfilename);
%calc mag sp heat/T of weak ferromag Eq(10) of Mohn PRB 1989
%t is red temperature=T/T_curie
%stonerishness is parameter tc=Tc/Tc(s) in Mohn paper where Tc is exptl Curie T,
%Tc(s) is Curie T predicted by pure Stoner theory
%thus stonerishness=1 corresponds to pure Stoner (no spin flucts)
%chi=0.018 (mu_b per f.u. per T) (from Pfleiderer data)

M0=0.17*9.274e-24*(6.022e23); %u_B per fu => Am^2 per mole
mu0_M0=0.17*(0.0231/0.1); %in Tesla
chi0=0.018*(0.0231/0.1); %=> dimless
Tc=27.5;
mag=-mu0_M0*M0/(2*chi0*Tc); %=gamma_m * T_c

result=zeros(size(t));
low_t=t(t<1); hi_t=t(t>=1);
result(t<1)=mag.*low_t.*(stonerishness.^2-((1-stonerishness.^2).^2)/5-3.*(low_t.^2).*...
    (stonerishness.^4)-(6/5).*(stonerishness^2).*low_t.*(1-stonerishness.^2));
result(t>=1)=mag.*hi_t.*stonerishness^2.*(9/5*(1-stonerishness^2).*hi_t+3/10*(1-stonerishness^2)^2/stonerishness^2);

if nargin==3
    ofile=fopen(outfilename,'w');
    fprintf(ofile,'Mohn-Hilscher specific heat calc for ''stonerishness'' t_c=%s\r\n',num2str(stonerishness));
    fprintf(ofile,'t=T/T_c,C_m\r\n');
    for row=1:length(result)
        fprintf(ofile,'%f,%f\r\n',t(row),result(row));
    end
    fclose(ofile);
end
