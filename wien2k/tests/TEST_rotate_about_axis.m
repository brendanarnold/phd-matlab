% Test
[x y z] = meshgrid(0:10, 0:10, 0);
[xf, yf, zf] = rotate_about_axis(x(:), y(:), z(:), [1 0 0], 35, 0);
scatter3(x(:), y(:), z(:), 'b');
xlabel('x');
ylabel('y');
zlabel('z');
hold on;
scatter3(xf, yf, zf, 'r');
[xf, yf, zf] = rotate_about_axis(xf, yf, zf, [1 0 0], -35, 0);
scatter3(xf, yf, zf, 'g');


figure;
[x y z] = meshgrid(0:10, 0:10, 0);
[xf, yf, zf] = rotate_about_axis(x(:), y(:), z(:), [1 1 0], 35, 0);
scatter3(x(:), y(:), z(:), 'b');
xlabel('x');
ylabel('y');
zlabel('z');
hold on;
scatter3(xf, yf, zf, 'r');
[xf, yf, zf] = rotate_about_axis(xf, yf, zf, [1 1 0], -35, 0);
scatter3(xf, yf, zf, 'g');


figure;
[x y z] = meshgrid(0:10, 0:10, 0);
[xf, yf, zf] = rotate_about_axis(x(:), y(:), z(:), [1 1 1], 35);
scatter3(x(:), y(:), z(:), 'b');
xlabel('x');
ylabel('y');
zlabel('z');
hold on;
scatter3(xf, yf, zf, 'r');
[xf, yf, zf] = rotate_about_axis(xf, yf, zf, [1 1 1], -35);
scatter3(xf, yf, zf, 'g');