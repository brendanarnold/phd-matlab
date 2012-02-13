NUM_PTS = 92;
MAX_X = 3;
MIN_X = 0;
MAX_Y = 3;
MIN_Y = 0;
MAX_Z = 1;
MIN_Z = 0;
FILENAME = 'FreeElectron_cartE.mat';


half_x = (MAX_X - MIN_X) / 2.0;
half_y = (MAX_Y - MIN_Y) / 2.0;

[kx, ky, kz] = meshgrid(linspace(MIN_X, MAX_X, NUM_PTS), ...
    linspace(MIN_Y, MAX_Y, NUM_PTS), ...
    linspace(MIN_Z, MAX_Z, NUM_PTS));

% Use free electron model

[ktheta, krho] = cart2pol(kx - half_x, ky - half_y);
cartE = krho .^ 2;

% Fermi level at 1, means krho = 1 which means area of pi
areas = [];
for z=1:NUM_PTS
    C = contour(squeeze(kx(1,:,1)), squeeze(ky(:,1,1)), squeeze(cartE(:,:,1)), [1 1]);
    areas = [areas polyarea(C(1,2:end), C(2,2:end))];
end
disp(sprintf('Area is: %.5f with stdev: %.7f', mean(areas), std(areas)));
disp(sprintf('Difference from pi is: %.10f', pi-mean(areas)));
% Even with 'perfect' example at 92pts, area is wrong in 5th s.f. or 4 d.p.

save(FILENAME, 'kx', 'ky', 'kz', 'cartE');



