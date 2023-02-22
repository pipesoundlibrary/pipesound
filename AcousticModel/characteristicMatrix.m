% CHARACTERISTIC MODE MATRIX REQUIRED BY KZSOLVER.

% COPYRIGHT: 
%
%   MIT License
%   © 2023 Dario Chiantello <pipesoundlibrary@gmail.com>
%

function A = characteristicMatrix(n, omega, kz, ss)

R1      = ss.R1_m;
R2      = ss.R2_m;
ro_l    = ss.ro_l_kg_m3;
c_phi   = ss.c_phi_m_s;
ro_s    = ss.ro_s_kg_m3;
c_gamma = ss.c_gamma_m_s;
c_psi   = ss.c_psi_m_s;


k_phi    = omega/c_phi;
k_gamma  = omega/c_gamma;
k_psi    = omega/c_psi;

q_phi   = sqrt(k_phi^2   - kz^2);
q_gamma = sqrt(k_gamma^2 - kz^2);
q_psi   = sqrt(k_psi^2   - kz^2);


 if n == 0
 
    %n=0
    d_31 = besselj(1,q_gamma*R1)*(-2*i*kz*q_gamma);
    d_32 = bessely(1,q_gamma*R1)*(-2*i*kz*q_gamma);
    d_33 = besselj(1,q_psi*R1)  *(kz^2 - q_psi^2);
    d_34 = bessely(1,q_psi*R1)  *(kz^2 - q_psi^2);

    d_41 = besselj(1,q_gamma*R2)*(-2*i*kz*q_gamma);
    d_42 = bessely(1,q_gamma*R2)*(-2*i*kz*q_gamma);
    d_43 = besselj(1,q_psi*R2)  *(kz^2 - q_psi^2);
    d_44 = bessely(1,q_psi*R2)  *(kz^2 - q_psi^2);

    W = ro_l*omega^2/(2*ro_s*c_psi^2) * besselj(0,q_phi*R1)/(q_phi*besselj(1,q_phi*R1));
    
    d_51 = besselj(0,q_gamma*R1)*(-kz^2 + k_psi^2/2)        + besselj(1,q_gamma*R1)*(-q_gamma*(1+R1*W)/R1);
    d_52 = bessely(0,q_gamma*R1)*(-kz^2 + k_psi^2/2)        + bessely(1,q_gamma*R1)*(-q_gamma*(1+R1*W)/R1);
    d_53 = besselj(0,q_psi*R1)  *(i*kz*q_psi)               + besselj(1,q_psi*R1)  *(-i*kz*(1+R1*W)/R1);
    d_54 = bessely(0,q_psi*R1)  *(i*kz*q_psi)               + bessely(1,q_psi*R1)  *(-i*kz*(1+R1*W)/R1);

    d_61 = besselj(0,q_gamma*R2)*(-kz^2 + k_psi^2/2)        + besselj(1,q_gamma*R2)*(-q_gamma/R2);
    d_62 = bessely(0,q_gamma*R2)*(-kz^2 + k_psi^2/2)        + bessely(1,q_gamma*R2)*(-q_gamma/R2);
    d_63 = besselj(0,q_psi*R2)  *(i*kz*q_psi)               + besselj(1,q_psi*R2)  *(-i*kz/R2);
    d_64 = bessely(0,q_psi*R2)  *(i*kz*q_psi)               + bessely(1,q_psi*R2)  *(-i*kz/R2);


     A = [d_31 d_32 d_33 d_34
          d_41 d_42 d_43 d_44
          d_51 d_52 d_53 d_54
          d_61 d_62 d_63 d_64];
  
else %n!=0

    d_11 = besselj(n,q_gamma*R1)*(2*i*n*(n-1)/R1^2)               + besselj(n+1,q_gamma*R1)*(-2*i*q_gamma*n/R1);
    d_12 = bessely(n,q_gamma*R1)*(2*i*n*(n-1)/R1^2)               + bessely(n+1,q_gamma*R1)*(-2*i*q_gamma*n/R1);
    d_13 = besselj(n,q_psi*R1)  *(-kz*q_psi)                      + besselj(n+1,q_psi*R1)  *(2*kz*(n+1)/R1);
    d_14 = bessely(n,q_psi*R1)  *(-kz*q_psi)                      + bessely(n+1,q_psi*R1)  *(2*kz*(n+1)/R1);
    d_15 = besselj(n,q_psi*R1)  *(2*i*n*(1-n)/R1^2 + i*q_psi^2)   + besselj(n+1,q_psi*R1)  *(-2*i*q_psi/R1);
    d_16 = bessely(n,q_psi*R1)  *(2*i*n*(1-n)/R1^2 + i*q_psi^2)   + bessely(n+1,q_psi*R1)  *(-2*i*q_psi/R1);

    d_21 = besselj(n,q_gamma*R2)*(2*i*n*(n-1)/R2^2)               + besselj(n+1,q_gamma*R2)*(-2*i*q_gamma*n/R2);
    d_22 = bessely(n,q_gamma*R2)*(2*i*n*(n-1)/R2^2)               + bessely(n+1,q_gamma*R2)*(-2*i*q_gamma*n/R2);
    d_23 = besselj(n,q_psi*R2)  *(-kz*q_psi)                      + besselj(n+1,q_psi*R2)  *(2*kz*(n+1)/R2);
    d_24 = bessely(n,q_psi*R2)  *(-kz*q_psi)                      + bessely(n+1,q_psi*R2)  *(2*kz*(n+1)/R2);
    d_25 = besselj(n,q_psi*R2)  *(2*i*n*(1-n)/R2^2 + i*q_psi^2)   + besselj(n+1,q_psi*R2)  *(-2*i*q_psi/R2);
    d_26 = bessely(n,q_psi*R2)  *(2*i*n*(1-n)/R2^2 + i*q_psi^2)   + bessely(n+1,q_psi*R2)  *(-2*i*q_psi/R2); 
    
    d_31 = besselj(n,q_gamma*R1)*(2*i*kz*n/R1)                    + besselj(n+1,q_gamma*R1)*(-2*i*kz*q_gamma);
    d_32 = bessely(n,q_gamma*R1)*(2*i*kz*n/R1)                    + bessely(n+1,q_gamma*R1)*(-2*i*kz*q_gamma);
    d_33 = besselj(n,q_psi*R1)  *(q_psi*n/R1)                     + besselj(n+1,q_psi*R1)  *(kz^2 - q_psi^2);
    d_34 = bessely(n,q_psi*R1)  *(q_psi*n/R1)                     + bessely(n+1,q_psi*R1)  *(kz^2 - q_psi^2);
    d_35 = besselj(n,q_psi*R1)  *(-i*kz*n/R1);
    d_36 = bessely(n,q_psi*R1)  *(-i*kz*n/R1);

    d_41 = besselj(n,q_gamma*R2)*(2*i*kz*n/R2)                   + besselj(n+1,q_gamma*R2)*(-2*i*kz*q_gamma);
    d_42 = bessely(n,q_gamma*R2)*(2*i*kz*n/R2)                   + bessely(n+1,q_gamma*R2)*(-2*i*kz*q_gamma);
    d_43 = besselj(n,q_psi*R2)  *(q_psi*n/R2)                    + besselj(n+1,q_psi*R2)  *(kz^2 - q_psi^2);
    d_44 = bessely(n,q_psi*R2)  *(q_psi*n/R2)                    + bessely(n+1,q_psi*R2)  *(kz^2 - q_psi^2);
    d_45 = besselj(n,q_psi*R2)  *(-i*kz*n/R2);
    d_46 = bessely(n,q_psi*R2)  *(-i*kz*n/R2);

    W = ro_l*omega^2/(2*ro_s*c_psi^2) *  besselj(n,q_phi*R1)/(-n/R1*besselj(n,q_phi*R1) + q_phi*besselj(n+1,q_phi*R1));

    d_51 = besselj(n,q_gamma*R1)*(-kz^2 + k_psi^2/2 + n*(1-n+R1*W)/R1^2)         + besselj(n+1,q_gamma*R1)*(-q_gamma*(1+R1*W)/R1);
    d_52 = bessely(n,q_gamma*R1)*(-kz^2 + k_psi^2/2 + n*(1-n+R1*W)/R1^2)         + bessely(n+1,q_gamma*R1)*(-q_gamma*(1+R1*W)/R1);
    d_53 = besselj(n,q_psi*R1)  *(i*kz*q_psi)                                    + besselj(n+1,q_psi*R1)  *(-i*kz*(n+1+R1*W)/R1);
    d_54 = bessely(n,q_psi*R1)  *(i*kz*q_psi)                                    + bessely(n+1,q_psi*R1)  *(-i*kz*(n+1+R1*W)/R1);
    d_55 = besselj(n,q_psi*R1)  *(n*(-1+n-R1*W)/R1^2)                            + besselj(n+1,q_psi*R1)  *(-n*q_psi/R1);
    d_56 = bessely(n,q_psi*R1)  *(n*(-1+n-R1*W)/R1^2)                            + bessely(n+1,q_psi*R1)  *(-n*q_psi/R1);

    d_61 = besselj(n,q_gamma*R2)*(-kz^2 + k_psi^2/2 + n*(1-n)/R2^2)              + besselj(n+1,q_gamma*R2)*(-q_gamma/R2);
    d_62 = bessely(n,q_gamma*R2)*(-kz^2 + k_psi^2/2 + n*(1-n)/R2^2)              + bessely(n+1,q_gamma*R2)*(-q_gamma/R2);
    d_63 = besselj(n,q_psi*R2)  *(i*kz*q_psi)                                    + besselj(n+1,q_psi*R2)  *(-i*kz*(n+1)/R2);
    d_64 = bessely(n,q_psi*R2)  *(i*kz*q_psi)                                    + bessely(n+1,q_psi*R2)  *(-i*kz*(n+1)/R2);
    d_65 = besselj(n,q_psi*R2)  * n*(-1+n)/R2^2                                  + besselj(n+1,q_psi*R2)  *(-n*q_psi/R2);
    d_66 = bessely(n,q_psi*R2)  * n*(-1+n)/R2^2                                  + bessely(n+1,q_psi*R2)  *(-n*q_psi/R2);


    A = [d_11 d_12 d_13 d_14 d_15 d_16
         d_21 d_22 d_23 d_24 d_25 d_26
         d_31 d_32 d_33 d_34 d_35 d_36
         d_41 d_42 d_43 d_44 d_45 d_46
         d_51 d_52 d_53 d_54 d_55 d_56
         d_61 d_62 d_63 d_64 d_65 d_66];  
    
 end
      