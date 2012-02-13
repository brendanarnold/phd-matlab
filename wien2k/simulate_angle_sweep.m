function [out, full_out] = simulate_angle_sweep(kx, ky, kz, energies, ef, axis, step, num_steps, interp_grid, slice_dk)
% Calculates a simulated angle sweep through a volume
% Returns 'out' a set of data suitable for saving as ascii and 'full_out'
% which also includes orbit vertices in 2d and 3d

% Factor 10e20 due to inverse Ang^2 times hbar/e
H_BAR_BY_E_METRE = 1.05457148e5/1.60217646;

% Max size will ever need the grid is the cartE cubic diagonal.
max_x = max(kx(:));
min_x = min(kx(:));
max_y = max(ky(:));
min_y = min(ky(:));
max_z = max(kz(:));
min_z = min(kz(:));

% Note the Brillouin zone corners used alter in determining the extent of
% the zone w.r.t. the slice.
bz_corners = [min_x min_y min_z;
    min_x min_y max_z;
    min_x max_y min_z;
    min_x max_y max_z;
    max_x min_y min_z;
    max_x min_y max_y;
    max_x max_y min_z;
    max_x max_y max_z];
bz_corners(:,1) = bz_corners(:,1) - (min_x + max_x)/2;
bz_corners(:,2) = bz_corners(:,2) - (min_y + max_y)/2;
bz_corners(:,3) = bz_corners(:,3) - (min_z + max_z)/2;

% Sweep angles
out = [];
full_out = {};
for angle = 0:step:step*(num_steps-1)
    disp(sprintf('Calculating for angle %.1f', angle));
    
    % Find out extent of BZ in slice direction - do this so can proportion 
    % slices so as to cover whole zone.
    [ux, uy, uz] = rotate_about_axis(0, 0, 1, axis, angle);
    slice_u = [ux; uy; uz];
    bz_extent = max(bz_corners * slice_u);
    
    for slice_shift = -bz_extent:slice_dk:bz_extent
        [orbits2d orbits3d] = extract_orbits(kx, ky, kz, energies, ef, angle, axis, slice_shift, interp_grid);
        for i = 1:length(orbits2d)
            orbit2d = orbits2d{i};
            orbit3d = orbits3d{i};
            if is_closed_orbit(orbit2d)
                area = polyarea(orbit2d(:,1), orbit2d(:,2));
                f = area * H_BAR_BY_E_METRE / (2 * pi);
                % Calculate the centre of the polygon
                [x_centre, y_centre] = poly_centre(orbit2d(:,1), orbit2d(:,2));
                [x_rot_centre, y_rot_centre, z_rot_centre] = rotate_about_axis(x_centre, y_centre, slice_shift, axis, angle);
                % Collect data for ascii file
                out = [out; [angle slice_shift f area x_rot_centre y_rot_centre z_rot_centre]];
                % Collect data for matlab file
                full_out{end+1}.area = area;
                full_out{end}.freq = f;
                full_out{end}.slice_shift = slice_shift;
                full_out{end}.angle = angle;
                full_out{end}.orbit2d = orbit2d;
                full_out{end}.orbit3d = orbit3d;
                full_out{end}.centre = [x_rot_centre y_rot_centre z_rot_centre];
            end
        end
    end
end

end