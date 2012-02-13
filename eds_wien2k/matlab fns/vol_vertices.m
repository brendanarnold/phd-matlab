function [faces,vertices]=vol_vertices(kparas,planenormals,incl_point);
%function [faces,vertices]=vol_vertices(kparas,planenormals,incl_point);
%
%generates vertices of surface enclosed by set of planenormals

vertices=[]; faces=[];
for plnum=1:numel(kparas)
    otherplanes=(1:numel(kparas))~=plnum;
    plist=boundingpoly(kparas(plnum),planenormals(plnum,:),kparas(otherplanes),planenormals(otherplanes,:),incl_point);
    for ptnum=1:size(plist,1)
        vertices(end+1,:)=plist(ptnum,:);
        if ptnum>size(faces,2) & size(faces,1)>0
            %bug in matlab: should pad faces with <max vert with nans but
            %this causes extra edges to appear
            faces=[faces faces(:,size(faces,2))];
        end
        faces(plnum,ptnum)=size(vertices,1);
    end
    padcols=(size(plist,1)+1):size(faces,2);
    if ~isempty(padcols)
        faces(end,padcols)=faces(end,size(plist,1));
    end
end
