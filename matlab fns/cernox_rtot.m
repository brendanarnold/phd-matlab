function T=cernox_RtoT(R,serialstr)
%Lakeshore calibration no. 238211 for CX-1030-SD serial X01098
%Lakeshore calibration no. 236915 for CX-1010-SD serial X01838

serialstrs={'X01098' 'X01838'};
calibsections=[3 1];
%Chebychev coeffs (first row is lowest range):
coeffs{1}=[1.035413,-1.189990,.629961,-.303532,.143654,-.063592,.023552,-.012732,.001089,-.002841;
    9.311905,-9.218812,3.201323,-.837415,.152070,-.008370,-.007076,.002724,0,0;
    59.356416,-50.003219,9.670778,-.94533,.036504,.01699,-.000421,-.005237,.003094,0];
coeffs{2}=[1.204000,-1.360708,0.567231,-0.194692,0.058114,-0.015259,0.002456,0.000009,0,0];
%logR limits:
logR_fitextents{1}=[2.8256911,5.8292394;2.2580526,2.9480072;1.8016259,2.3133477];
logR_fitextents{2}=[2.73033130819,3.29154645053];
%values of R over which fits are useful (corresponding to
%0.3K,3K,20.1K,100K
R_limits{1}=[245010,760.20,194.16,72.109];
R_limits{2}=[1956.80,537.4416];

therm_num=find(strcmp(serialstrs,serialstr));
if ~isempty(therm_num)
    if (R<R_limits{therm_num}(calibsections(therm_num)+1) | R>R_limits{therm_num}(1))
        range=0; %outside calib range
    else
        range=find(R>R_limits{therm_num},1,'first')-1;
    end
    
    T=0;
    if (range~=0)
        logR=log10(R);
        X=((logR-logR_fitextents{therm_num}(range,1))-(logR_fitextents{therm_num}(range,2)-logR))/(logR_fitextents{therm_num}(range,2)-logR_fitextents{therm_num}(range,1));
    
        for term=0:1:9
            T=T+coeffs{therm_num}(range,term+1)*cos(term*acos(X));
        end
    end
else
    'Serial number not recognised'
end