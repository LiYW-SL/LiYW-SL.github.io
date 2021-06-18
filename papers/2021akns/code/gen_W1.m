function [y] = gen_W1(E,e)
type = iscolumn(e);
if type
    e = e';
end
N = (length(E)-1)/2;
trivial_E = (-N:N)*2*pi;
trivial_E(N+1) = 1;
poly_fun1 = 1 - e./E;
poly_fun2 = 1 - e./trivial_E';
poly_fun2(N+1,:) = 1/2*e;
poly_fun = poly_fun1./poly_fun2;
y = prod(poly_fun,1).*sin(1/2*e);
if type
    y = y';
end



end

