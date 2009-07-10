function [map,grayscale]=create_colormap(vals,cols,flag)

% function [map,grayscale]=create_colormap(vals,cols,flag)
% create RGB interpolated colourmap matrix (np,3) from a set of basic colours
% default colormap is 
%    vals=[0 0.1 0.2 0.4 0.5 1];
%    cols=str2mat('white','blue','darkgreen','yellow','red','black'); 
% edit subfunction 
% function val=rbg(colour_name) 
% at the end of this function to include RGB conversion for colours not 
% at corners of the RGB cube 
% see http://www.soapplab.auckland.ac.nz/info/formats/rgbcube.htm
% more RGB colours from http://web.njit.edu/~kevin/rgb.txt.html
% Radu Coldea 03-Sep-2004

if ~exist('vals','var')|isempty(vals)|~exist('cols','var')|~ischar(cols),
    % default colormap is TVtueb original 
    % tvtueb diffraction programme at HMI Berlin
    % http://www.hmi.de/bensc/misc/flat-cone/tvtueb/index_en.html
    vals=[0 0.1 0.2 0.4 0.5 1];
    cols=str2mat('white','blue','darkgreen','yellow','red','black');
end 
if ~exist('flag','var')|isempty(flag)|~islogical(flag),
    flag=(1<0); % do not show grayscale and RGB cube representation
end    
% check that the numbers of colours matches the numbers of values 
n=length(vals);
if size(cols,1)~=n,
    disp(sprintf('Number of interpolating colours (%d) and values (%d) are not compatible.',...
        size(cols,1),n));
    map=[];
    return;
end    
% identify RGB equivalent of the basic colours
rgb_cols=zeros(n,3);
for i=1:n,
    rgb_value=rbg(cols(i,:));
    if ~isempty(rgb_value),
        rgb_cols(i,:)=rgb_value;
    else
        disp(sprintf('could not find RGB equivalent for colour %s',...
            cols(i,:)));
        map=[];
        return;
    end    
end
%rgb_cols,
% construct interpolating colourmap with np scale values, increase this
% number to have finer levels
np=200;
t=linspace(0,1,np);
map=zeros(np,3);
index=interp1(vals,1:length(vals),t);
for i=1:np,
    i1=floor(index(i));
    i2=ceil(index(i));
    map(i,:)=rgb_cols(i1,:)+rem(index(i),1)*(rgb_cols(i2,:)-rgb_cols(i1,:));
end

% calculate the equivalent grayscale 
% determine the effective luminance of a pixel as per  
% Y=0.3RED+0.59GREEN+0.11Blue
% http://www.bobpowell.net/grayscale.htm
y=0.3*map(:,1)+0.59*map(:,2)+0.11*map(:,3);
grayscale=[y y y];
if flag, % visualize colour map
    colormap2rgbcube(map);
end

function val=rbg(colour_name)

% returns RGB 3-value vector equivalent of a colour name
% basic colours at the corners of RGB cube defined as per 
% http://www.soapplab.auckland.ac.nz/info/formats/rgbcube.htm
% rgbcube.jpg
% colour converter at http://www.draac.com/colorconvert.html

switch lower(deblank(colour_name)),
    case {'black'},  val=[0 0 0];
    case {'red'},    val=[1 0 0];
    case {'magenta'},val=[1 0 1];
    case {'blue'},   val=[0 0 1];
    case {'green'},  val=[0 1 0];
    case {'green4'}, val=[0 139/255 0];
    case {'greenyellow'}, val=[173 255 47]/255;
    case {'yellow'}, val=[1 1 0];    
    case {'white'},  val=[1 1 1];
    case {'cyan'},   val=[0 1 1];
    case {'orange'}, val=[1 165/255 0];
    case {'indigo'}, val=[75 0 130]/255;        
    case {'violet'}, val=[238 130 238]/255;
    case {'darkgreen'}, val=[0 0.5 0];  
    case {'rgb1'}, val=[0   0    1];  
    case {'rgb2'}, val=[0   1    1];
    case {'rgb3'}, val=[1   1    1];        
    otherwise, val=[];        
end        