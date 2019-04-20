%--------------------------------------------------------------------------
% Evaluates the constitutive tensor (in Voigt notation) for material type 9.
%--------------------------------------------------------------------------
function c   = ctens9(kinematics,properties,dim)
mu1          = properties(1);
mu2          = properties(2);
k            = properties(3);
F            = kinematics.F;
F_bar        = J^(-1/3)*F;
J            = kinematics.J;
b            = kinematics.b;
b_bar        = F_bar*(F_bar.')
Ib           = kinematics.Ib;
C            = F'*F;
C_bar        = J^(-2/3)*C;
C_bar_inv    = C_bar^-1;
T            = kinematics.n; 
I1           = Ib;
I1_bar       = J^(-2/3)*I1
I2           = 0.5.*(Ib.^2 - trace(b.^2));
I2_bar       = J^(-2/3)*I2;
I3           = J.^2;
I3_bar       = J^(-2/3)*I3
I            = eye(size(J))

c            = zeros(dim,dim,dim,dim);
for l=1:dim
    for k=1:dim
        for j=1:dim
            for i=1:dim
                Cvol(i,j,k,l) = J^(-1/3)*(2*J-1)*k*C_bar_inv(i,j)*C_bar_inv(k,l)-2*J^(-1/3)*(J-1)*k*C_bar_inv(i,l)*C_bar_inv(j,k);
                Ciso(i,j,k,l) = 2*J^(-4/3)*mu2*(I*(I-1)-(2/3)*J^(-4/3))*(mu1+2*mu2*I1_bar)*(C_bar_inv(i,j)*I+I*C_bar_inv(k,l))...
                    + (4/3)*J^(-4/3)*mu2*(C_bar_inv(i,j) * C_bar(k,l) + C_bar(i,j)*C_bar_inv(k,l))...
                    + (2/9)*J^(-4/3)*(mu1*I1_bar+4*mu2*I2_bar)*(C_bar_inv(i,j) * C_bar_inv(k,l))...
                    + (2/3)*J^(-4/3)*(mu1*I1_bar+2*mu2*I2_bar)*(1/2)*((C_bar_inv(i,j) * C_bar_inv(j,l))*(C_bar_inv(i,l) * C_bar_inv(j,k)));
                 
                c(i,j,k,l) = Cvol(i,j,k,l) + Ciso(i,j,k,l);
            end
        end
    end    
end
end