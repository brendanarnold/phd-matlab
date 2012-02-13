function fid=safe_fopen(filename,openmode);
%performs fopen until a valid fid is returned
%avoids crash e.g. due to user copying partial file

while(1) %loop in case open unsuccessful 
    fid=fopen(filename,openmode);
    if (fid>0) 
        break;
    end
end

