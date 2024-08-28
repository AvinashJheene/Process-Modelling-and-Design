function d=solve(x_fit,x,y)
d=(x_fit(1).*x)./(1 + (x_fit(2).*x)+(x_fit(3).*x.^2))-y;
end