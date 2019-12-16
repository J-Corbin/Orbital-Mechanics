function plot_planet(radius, color, location)
[xx, yy, zz] = sphere(100);
surface(radius*xx + location(1), radius*yy + location(2),...
    radius*zz + location(3), 'FaceColor', color, 'EdgeColor', 'none');
end