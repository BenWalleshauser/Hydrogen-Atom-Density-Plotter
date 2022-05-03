function PL = Legendre(L)

syms x PL

rhs = (x^2-1)^L;
for i = 1:L
    rhs = diff(rhs);
end

PL = (1/(2^L*factorial(L)))*rhs;
PL = simplify(PL);

end
