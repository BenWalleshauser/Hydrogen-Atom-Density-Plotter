function L = AssociatedLaguerre(q,p,N,a)

syms L x r 

rhs = exp(-x)*x^(p+q);

for i = 1:q
    rhs = diff(rhs);
end

L = x^(-p)*exp(x)*rhs/(factorial(q));
L = simplify(L);

L = subs(L,x,2*r/(N*a));
end
