Creating a 3D plot in MATLAB is straightforward and can be done using functions like plot3, mesh, surf, or scatter3, depending on the type of 3D visualization you need. Here's an example for each:

1. Line Plot in 3D (plot3)
Copy the code
% Define data
x = 0:0.1:10;
y = sin(x);
z = cos(x);

% Create 3D line plot
plot3(x, y, z, 'LineWidth', 2);
grid on;
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('3D Line Plot');

2. 3D Surface Plot (surf)
Copy the code
% Define grid
[X, Y] = meshgrid(-5:0.5:5, -5:0.5:5);
Z = X.^2 + Y.^2;

% Create surface plot
surf(X, Y, Z);
shading interp; % Smooth shading
colormap jet;   % Color scheme
colorbar;       % Add color bar
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('3D Surface Plot');

3. 3D Scatter Plot (scatter3)
Copy the code
% Generate random data
x = rand(1, 100) * 10;
y = rand(1, 100) * 10;
z = rand(1, 100) * 10;

% Create scatter plot
scatter3(x, y, z, 50, z, 'filled'); % Color based on z-values
grid on;
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('3D Scatter Plot');


These examples should help you get started with 3D plotting in MATLAB. You can customize the plots further by adjusting properties like colors, markers, and labels.