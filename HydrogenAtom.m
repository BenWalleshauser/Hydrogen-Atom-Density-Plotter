%B.W.
%Density plots for hydrogen atom
%% Picking Eigenstate
clear
clc

N = 4;
L = 2;
M = 1;

%% Finding wavefunction
syms a r theta

%Wavefunction, not including phi term for convenience
psi = sqrt((2/(N*a))^3*factorial(N-L-1)/(2*N*factorial(N+L)))*exp(-r/(N*a))*(2*r/(N*a))^L*AssociatedLaguerre(N-L-1,2*L+1,N,a)*SphericalHarmonics(L,M);
psi_squared = psi^2;

%% Evaluating
%Bohr radius
a_val = 0.529e-10;

%Discretize r and theta
del_theta = 0.1;
del_r = 10e-12;

%Maximum value of r
max_r = a_val*18;

r_val = 0:del_r:max_r;
theta_val = 0:del_theta:pi;

%Matrix of probability densities:
A = zeros(length(theta_val),length(r_val));


%Evaluating probability density
it = 1;
f = waitbar(0,'Please wait...');
for i = 1:length(r_val)
    for j = 1:length(theta_val)
        
        dP = psi_squared;
        dP_val = subs(dP,[a r theta],[a_val r_val(i) theta_val(j)]);

        A(j,i) = dP_val;
        it = it+1;
        waitbar(it/(length(r_val)*length(theta_val)),f,['Percent Done: ',num2str(100*it/(length(r_val)*length(theta_val)))]);
    end
end
close(f)
%% Plotting

figure(1)
for i = 1:length(r_val)-1
    hold on
     R=[r_val(i);r_val(i+1);r_val(i+1);r_val(i)]*ones(1,length(theta_val));
     
     theta_val2 = [theta_val theta_val(end)+del_theta];
     theta_val_adj=[theta_val2(1:end-1);theta_val2(1:end-1);theta_val2(2:end);theta_val2(2:end)];
     y = R.*cos(theta_val_adj);   
     z = R.*sin(theta_val_adj); 
    intensity=ones(4,1)*A(:,i)';
    patch(z,y,intensity,'EdgeColor','none'), axis equal
    patch(-z,y,intensity,'EdgeColor','none'), axis equal    %utilizing cylindrical symmetry
end

xticks([-max_r:a_val:max_r]);
Str(1) = "hold";
for i = -max_r/a_val:max_r/a_val
    Str(i+max_r/a_val+1) = [num2str(i),'a'];
end
xticklabels(Str)

yticks([-max_r:a_val:max_r]);
Str(1) = "hold";
for i = -max_r/a_val:max_r/a_val
    Str(i+max_r/a_val+1) = [num2str(i),'a'];
end
yticklabels(Str);
grid on

xlim([-max_r max_r])
ylim([-max_r max_r])
xlabel('y')
ylabel('z')
title(['Cross Section of Density Plot for \Psi_',num2str(N),'_',num2str(L),'_',num2str(M)])
hold off


%% Generating Functions

function L = AssociatedLaguerre(q,p,N,a)

syms L x r 

rhs = exp(-x)*x^(p+q);

for i = 1:q
    rhs = diff(rhs);
end

L = x^(-p)*exp(x)*rhs/(factorial(q));
L = simplify(L);

%xval = 2r/na
L = subs(L,x,2*r/(N*a));
end


function Y = SphericalHarmonics(L,M)

syms x theta
%ignoring phi bc multiplying by conjugate
PLM = AssociatedLegendre(L,M);


Y = sqrt((2*L+1)*factorial(L-M)/(4*pi*factorial(L+M)))*subs(PLM,x,cos(theta));

end


function P = AssociatedLegendre(L,M)
%xval = cos(theta)
syms x 
PL = Legendre(L);
rhs = PL;

for i = 1:M
    rhs = diff(rhs);
end

P = (-1)^M*(1-x^2)^(M/2)*rhs;

end

function PL = Legendre(L)

syms x PL

rhs = (x^2-1)^L;
for i = 1:L
    rhs = diff(rhs);
end

PL = (1/(2^L*factorial(L)))*rhs;
PL = simplify(PL);

end




