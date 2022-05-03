function Y = SphericalHarmonics(L,M)

syms x theta
%ignoring phi bc multiplying by conjugate
PLM = AssociatedLegendre(L,M);


Y = sqrt((2*L+1)*factorial(L-M)/(4*pi*factorial(L+M)))*subs(PLM,x,cos(theta));

end