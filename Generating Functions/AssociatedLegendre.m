function P = AssociatedLegendre(L,M)

syms x 
PL = Legendre(L);
rhs = PL;

for i = 1:M
    rhs = diff(rhs);
end

P = (-1)^M*(1-x^2)^(M/2)*rhs;

end
