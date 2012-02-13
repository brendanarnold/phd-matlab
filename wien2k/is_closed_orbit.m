function [result] = is_closed_orbit(orbit)
% Accepts Nx2 array of x,y values. Guesses at whether the orbit is closed
% or not.

% Arbitrary factor to determine allowed difference in spacing from mean
DIFF_FACTOR = 2;

% Takes second to last point, last point, first and second point
if size(orbit, 1) < 4
    result = true;
    return
end
dxy = diff(orbit([end 1:end], :));
spacings = hypot(dxy(:,1), dxy(:,2));
avg_spacings = mean(spacings(2:end));
if spacings(1) > (avg_spacings * DIFF_FACTOR)
    result = false;
else
    result = true;

end