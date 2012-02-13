function [dirnresult, minwidthresult]=getinplanedirn(normal,bandsdata)
%note this will need mod so kpara does not need to be zero

%find direction in slice-plane at kpara=0 parallel to which width of
%polygon is minimum. this will be chosen as direction corresponding to row ->
%row+1 in the 2D slice matrix
%Make sure normal is normalised
%bandsdata.searchvol should be in Cartesian units
phimin=fminbnd(@(phi) getminwidth(phi,normal,bandsdata),0,pi);
minwidth=getminwidth(phimin,normal,bandsdata);
dirnresult{1}=planedirn(normal,phimin);
dirnresult{2}=cross(normal,dirnresult{1});
dirnresult{2}=dirnresult{2}/norm(dirnresult{2});

minwidthresult=minwidth;        

function result=getminwidth(phi,normal,bandsdata);
%finds intersection of faces of searchvol with plane perp to normal (and kpara=0 at present)
q=planedirn(normal,phi);
for facenum=1:size(bandsdata.searchvol,1)
    bndrynormal=bandsdata.searchvol(facenum,:);
    costheta=dot(bndrynormal/norm(bndrynormal),normal);       
    r0(facenum,:)=(costheta*norm(bndrynormal)*normal-bndrynormal);
    denom=(costheta^2-1);   
    if abs(denom)>1e-10
        r0(facenum,:)=r0(facenum,:)/denom;
    else
        r0(facenum,:)=sign(r0(facenum,:))*sign(denom)*Inf;
    end
    lamda(facenum)=dot(r0(facenum,:),bndrynormal);
    denom=dot(q,bndrynormal);
    if abs(denom)>1e-10
        lamda(facenum)=lamda(facenum)/denom;
    else
        lamda(facenum)=sign(lamda(facenum))*sign(denom)*Inf;
    end
end
%result=max(lamda)-min(lamda); %need to think here how to allow for [0 0 0] point not in volume
result=min(lamda(lamda>0))-max(lamda(lamda<0)); %(nov 06)

function result=planedirn(normal,phi);
%returns a dirn in plane perp to normal whose azimuthal angle is
%parameterized by phi
result=cross(normal,[cos(phi) sin(phi) 0]);
if isequal(result,[0 0 0])
    result=cross(normal,[0 cos(phi) sin(phi)]);
end
result=result/norm(result);
