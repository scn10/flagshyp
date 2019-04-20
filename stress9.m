%--------------------------------------------------------------------------
% Evaluates the Cauchy stress tensor for material type 9.
%--------------------------------------------------------------------------
function Cauchy = stress9(kinematics,properties,dim)
mu1          = properties(2);
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
Svol         = J^(1/3)*(J-1)*k*C_bar_inv;
Siso         = J^(-2/3)*((-1/3)*-(mu1*I1_bar+2*mu2*I2_bar)*C_bar_inv +(mu1+mu2*I1_bar)*I-mu2*C_bar)

PK2          = Svol + Siso;
Cauchy       = J^-1 * PK2;

end

