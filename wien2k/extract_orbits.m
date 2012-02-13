function [orbits2d orbits3d] = extract_orbits(kx, ky, kz, energies, fe, angle, axis, slice_shift, interp_grid_pts)

% Rotates mesh around centre point by <angle> about <rot_axis>. Then shifts
% mesh by <slice_shift>
%
% Automatically expands mesh to fill Brillouin zone

% Max size will ever need the grid is the cartE cubic diagonal.
max_x = max(kx(:));
min_x = min(kx(:));
max_y = max(ky(:));
min_y = min(ky(:));
max_z = max(kz(:));
min_z = min(kz(:));
cubic_diagonal = sqrt((max_x - min_x)^2 + (max_y - min_y)^2 + (max_z - min_z)^2);

% Will need centre of the zone later
bz_centre = [(max_x + min_x)/2 (max_y + min_y)/2 (max_z + min_z)/2];

% Make huge grids
[ix, iy] = meshgrid(linspace(-cubic_diagonal/2, cubic_diagonal/2, interp_grid_pts), ...
    linspace(-cubic_diagonal/2, cubic_diagonal/2, interp_grid_pts));
iz = zeros(size(ix));

% Perform shift, then rotate
iz = iz + slice_shift;
[rot_x, rot_y, rot_z] = rotate_about_axis(ix(:), iy(:), iz(:), axis, angle);

% Shift to bz co-ordinates
rot_x = rot_x + bz_centre(1);
rot_y = rot_y + bz_centre(2);
rot_z = rot_z + bz_centre(3);

% Bypass next step if all the values lie outside of Brillouin zone
in_zone = (rot_x > min_x) & (rot_x < max_x) & (rot_y > min_y) & (rot_y < max_y) & (rot_z > min_z) & (rot_z < max_z);
if ~any(in_zone(:))
    orbits2d = {};
    orbits3d = {};
    return
end

% Interpolate onto energies
ienergies = interp3(kx, ky, kz, energies, rot_x, rot_y, rot_z, 'linear', NaN);

% Put onto a grid
ienergies = reshape(ienergies, [interp_grid_pts interp_grid_pts]);

% Take contour and extract orbits in 2d plane
C = contourc(ix(1,:), iy(:,1), ienergies, [fe fe]);
orbits2d = {};
orbits3d = {};
while ~isempty(C)
    num_pts = C(2,1);
    orbits2d{end+1} = C(:,2:num_pts+1)';
    % Put the points into co-ordinates of the BZ
    contour_z = iz(1,1) * ones(size(orbits2d{end}, 1), 1);
    [pts_x, pts_y, pts_z] = rotate_about_axis(orbits2d{end}(:,1), orbits2d{end}(:,2), contour_z, axis, angle);
    orbits3d{end+1} = [pts_x pts_y pts_z];
    % Remove this set of points from the contour result
    C = C(:,num_pts+2:end);
end


end



