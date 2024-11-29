function integral = integrate(curve)
% This function integrates one curve to calculate the performance of the solver.

integral = 0;
% Since the stair function is right continuous, when we calculate the performance, 
% we need to multiply the height of the left point by the length of the stair.
for i = 2 : size(curve, 2)
    integral = integral + (curve(1, i) - curve(1, i-1))...
        *curve(2, i-1);
end

end

