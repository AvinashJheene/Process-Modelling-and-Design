function F = minSolvent(x, coeff1, coeff2, x2, y2)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
F(1) = (x(2) - y2) - ( (x(1) -x2) * coeff1 / ((1 + coeff2*x(1))^2));
F(2) = (coeff1*x(1)/(1 + coeff2*x(1))) - x(2);
end