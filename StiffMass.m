function [m11, m12, m21, m22, n11, n12, n21, n22, k11, k12, k21, k22, K, M, N] = StiffMass (A1, A2, den1, den2, c1, s1, c2, s2,L1, L2, E1, E2, Mcons)  
  %Lumped mass Matrix
m11 = ((den1*A1*L1*c1^2)/2 + (den2*A2*L2*c2^2)/2)+ Mcons;
m12 = ((den1*A1*L1*c1*s1)/2 + (den2*A2*L2*c2*s2)/2);
m21 = ((den1*A1*L1*c1*s1/2) + (den2*A2*L2*c2*s2/2));
m22 = ((den1*A1*L1*s1^2)/2 + (den2*A2*L2*s2^2)/2) + Mcons;
M = [m11 m12;
     m21 m22];
 
 %Distributed mass Matrix
n11 = ((den1*A1*L1*c1^2)/3 + (den2*A2*L2*c2^2)/3) + Mcons;
n12 = ((den1*A1*L1*c1*s1)/3 + (den2*A2*L2*c2*s2)/3);
n21 = ((den1*A1*L1*c1*s1/3) + (den2*A2*L2*c2*s2/3));
n22 = ((den1*A1*L1*s1^2)/3 + (den2*A2*L2*s2^2)/3) + Mcons;
N = [n11 n12;
     n21 n22];

 %Stiffness Matrix
 k11 = ((E1*A1*c1^2)/L1 + (E2*A2*c2^2)/L2);
 k12 = ((E1*A1*c1*s1)/L1 + (E2*A2*c2*s2)/L2);
 k21 = ((E1*A1*c1*s1)/L1 + (E2*A2*c2*s2)/L2);
 k22 = ((E1*A1*s1^2)/L1 + (E2*A2*s2^2)/L2);
 K = [k11 k12;
     k21 k22];
end
 
