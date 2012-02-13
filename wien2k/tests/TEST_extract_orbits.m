load('FreeElectron_cartE.mat');

fe = 1;
angle = 0;
axis = [1 0 0];
slice_shift = 0;
interp_grid_pts = 100;

disp('Normal test ...')

[orbits2d orbits3d] = extract_orbits(kx, ky, kz, cartE, fe, angle, axis, slice_shift, interp_grid_pts);
area = polyarea(orbits2d{1}(:,1), orbits2d{1}(:,2));
disp(sprintf('Area should be close to pi: %.10f', area));

disp('Rotation test...');
for angle = 0:45
    [orbits2d orbits3d] = extract_orbits(kx, ky, kz, cartE, fe, angle, axis, slice_shift, interp_grid_pts);
    if isempty(orbits2d)
        disp('Area is empty');
    else
        area = polyarea(orbits2d{1}(:,1), orbits2d{1}(:,2));
        disp(sprintf('Area should be close to pi: %.10f', area*cos(angle*pi/180)));
    end
end

angle = 0;
disp('Shift test ...');
for slice_shift = 0:0.1:0.5
    [orbits2d orbits3d] = extract_orbits(kx, ky, kz, cartE, fe, angle, axis, slice_shift, interp_grid_pts);
    if isempty(orbits2d)
        disp('Area is empty');
    else
        area = polyarea(orbits2d{1}(:,1), orbits2d{1}(:,2));
        disp(sprintf('Area should be close to pi: %.10f', area));
    end
end


    