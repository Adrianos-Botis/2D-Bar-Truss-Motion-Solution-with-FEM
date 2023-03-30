%-----Main Code---------------------
%Input Data for the bars-all units SI
clear all
clc
den1 = 8000;                           %Density  kg/m^3
den2 = 8000;
A1 = 0.0019625;                        %Area m^2
A2 = 0.0019625;
L1 = 1;                                %Length m
L2 = 2; 
E1 = 10^9;                             %Elastic Modulus Pa
E2 = 10^9;
Angle1 =   20;                          %Angle of 1st Bar in Degrees
Angle2 =  120;                          %Angle of 1st Bar in Degrees
InitC = [0.01 0.01 0 0];                %Vector of initial conditions [u, du, v, dv]
Ntms = 2;                               %Ntms times the total truss mass
Mcons = Ntms*(den1*A1*L1 + den2*A2*L2); %Concentrated extra mass  

if (Angle1 == Angle2)
    c1 = 1;
    s1 = 0;
    c2 = 1;
    s2 = 0;   
c1n = cosd(Angle1);
s1n = sind(Angle1);
inc(1) = InitC(1)*c1n + InitC(3)*s1n;
inc(2) = InitC(2)*c1n + InitC(4)*s1n;
inc(3) = -InitC(1)*s1n + InitC(3)*c1n;
inc(4) = -InitC(2)*s1n + InitC(4)*c1n;
%Matrices
[m11, m12, m21, m22, n11, n12, n21, n22, k11, k12, k21, k22, K, M, N] = StiffMass(A1, A2, den1, den2, c1, s1, c2, s2,L1, L2, E1, E2, Mcons);

 %Eigenvalues Calculation
 EigLumped = eig(K,M)   %Lumped mass
 EigDistr = eig(K,N)    %Consistent mass
 
 %Solution of system
    tspan = [0 0.1];
    y0 = inc;
    opts = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);
    [time,yloc] = ode45(@(t,y) odefun(t,y,m11,m12,m21,m22,k11,k12,k21,k22),tspan, y0,opts);     %Lumped 
    [ts,zloc] = ode45(@(t,y) odefun(t,y,n11,n12,n21,n22,k11,k12,k21,k22),tspan, y0,opts);    %Consistent
   
    
    %Calculation of internal stresses
    [numRows,numColls] = size(yloc);
    for i = 1:numRows
    y(i,1) = c1n*yloc(i,1) - s1n*yloc(i,3);
    y(i,2) = c1n*yloc(i,2) - s1n*yloc(i,4);
    y(i,3) = s1n*yloc(i,1) + c1n*yloc(i,3);
    y(i,4) = s1n*yloc(i,2) + c1n*yloc(i,4);
    DL1(i) = y(i,1); 
    DL2(i) = y(i,3);
    Stress1(i) = E1*DL1(i)/L1;
    Stress2(i) = E2*DL2(i)/L2;
    end
    
    [numz, collz] = size(zloc);
    for i = 1:numz
    z(i,1) = c1n*zloc(i,1) - s1n*zloc(i,3);
    z(i,2) = c1n*zloc(i,2) - s1n*zloc(i,4);
    z(i,3) = s1n*zloc(i,1) + c1n*zloc(i,3);
    z(i,4) = s1n*zloc(i,2) + c1n*zloc(i,4);
    end
    
elseif  (abs(Angle1 - Angle2) > 0 ) && (abs(Angle1 - Angle2) ~= 180 ) && (abs(Angle1 - Angle2) ~= 90)
    c1 = cosd(Angle1);
    s1 = sind(Angle1);
    c2 = cosd(Angle2);
    s2 = sind(Angle2);
    
 %Matrices
[m11, m12, m21, m22, n11, n12, n21, n22, k11, k12, k21, k22, K, M, N] = StiffMass(A1, A2, den1, den2, c1, s1, c2, s2,L1, L2, E1, E2, Mcons);

 %Eigenvalues Calculation
 EigLumped = eig(K,M)   %Lumped mass
 EigDistr = eig(K,N)    %Consistent mass
 
 %Solution of system
    tspan = [0 0.5];
    y0 = InitC;
    opts = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);
    [time,y] = ode45(@(t,y) odefun1(t,y,M,K),tspan, y0,opts);     %Lumped 
    [ts,z] = ode45(@(t,y) odefun1(t,y,N,K),tspan, y0,opts);       %Consistent
   
    
    %Calculation of internal stresses
    [numRows,numColls] = size(y);
    for i = 1:numRows
    DL1(i) = -c1*y(i,1) - s1*y(i,2); 
    DL2(i) = -c2*y(i,1) - s2*y(i,2);
    Stress1(i) = E1*DL1(i)/L1;
    Stress2(i) = E2*DL2(i)/L2;
    end
    
elseif (abs(Angle1 - Angle2) == 180 )
    c1 = 1;
    s1 = 0;
    c2 = -1;
    s2 = 0;   
c1n = cosd(Angle1);
s1n = sind(Angle1);
inc(1) = InitC(1)*c1n + InitC(3)*s1n;
inc(2) = InitC(2)*c1n + InitC(4)*s1n;
inc(3) = -InitC(1)*s1n + InitC(3)*c1n;
inc(4) = -InitC(2)*s1n + InitC(4)*c1n;
 
%Matrices
[m11, m12, m21, m22, n11, n12, n21, n22, k11, k12, k21, k22, K, M, N] = StiffMass(A1, A2, den1, den2, c1, s1, c2, s2,L1, L2, E1, E2, Mcons);

 %Eigenvalues Calculation
 EigLumped = eig(K,M)   %Lumped mass
 EigDistr = eig(K,N)    %Consistent mass
 
 %Solution of system
    tspan = [0 0.1];
    y0 = inc;
    opts = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);
    [time,yloc] = ode45(@(t,y) odefun(t,y,m11,m12,m21,m22,k11,k12,k21,k22),tspan, y0,opts);     %Lumped 
    [ts,zloc] = ode45(@(t,y) odefun(t,y,n11,n12,n21,n22,k11,k12,k21,k22),tspan, y0,opts);    %Consistent
   
    
    %Calculation of internal stresses
    [numRows,numColls] = size(yloc);
    for i = 1:numRows
    y(i,1) = c1n*yloc(i,1) - s1n*yloc(i,3);
    y(i,2) = c1n*yloc(i,2) - s1n*yloc(i,4);
    y(i,3) = s1n*yloc(i,1) + c1n*yloc(i,3);
    y(i,4) = s1n*yloc(i,2) + c1n*yloc(i,4);
    DL1(i) = y(i,1); 
    DL2(i) = y(i,3);
    Stress1(i) = E1*DL1(i)/L1;
    Stress2(i) = E2*DL2(i)/L2;
    end
    
    [numz, collz] = size(zloc);
    for i = 1:numz
    z(i,1) = c1n*zloc(i,1) - s1n*zloc(i,3);
    z(i,2) = c1n*zloc(i,2) - s1n*zloc(i,4);
    z(i,3) = s1n*zloc(i,1) + c1n*zloc(i,3);
    z(i,4) = s1n*zloc(i,2) + c1n*zloc(i,4);
    end
    
elseif (abs(Angle1 - Angle2) == 90)
    c1 = cosd(Angle1);
    s1 = sind(Angle1);
    c2 = cosd(Angle2);
    s2 = sind(Angle2);
    
%Matrices
[m11, m12, m21, m22, n11, n12, n21, n22, k11, k12, k21, k22, K, M, N] = StiffMass(A1, A2, den1, den2, c1, s1, c2, s2,L1, L2, E1, E2, Mcons);
%Eigenvalues Calculation
 EigLumped = eig(K,M)   %Lumped mass
 EigDistr = eig(K,N)    %Consistent mass
 
 %Solution of system
    tspan = [0 0.1];
    y0 = InitC;
    opts = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);
    [time,y] = ode45(@(t,y) odefun(t,y,m11,m12,m21,m22,k11,k12,k21,k22),tspan, y0,opts);     %Lumped 
    [ts,z] = ode45(@(t,y) odefun(t,y,n11,n12,n21,n22,k11,k12,k21,k22),tspan, y0,opts);    %Consistent
   
    %Calculation of internal stresses
    [numRows,numColls] = size(y);
    for i = 1:numRows
    DL1(i) = -c1*y(i,1) - s1*y(i,3); 
    DL2(i) = -c2*y(i,1) - s2*y(i,3);
    Stress1(i) = E1*DL1(i)/L1;
    Stress2(i) = E2*DL2(i)/L2;
    end 
   
end    

    if (abs(Angle1 - Angle2) > 0 ) && (abs(Angle1 - Angle2) ~= 180 ) && (abs(Angle1 - Angle2) ~= 90)
        j = 2;
        k = 3;
    else
        j = 3;
        k = 2;
    end  
 
   %Results Plotting
      figure(1)
      plot(time, y(:,1), ts, z(:,1))
      figure(2)
      plot(time,y(:,j), ts, z(:,j))
      figure(3)
      plot(time,Stress1,'blue',time, Stress2,'red')
      figure(4)
      plot(time,y(:,k),ts,z(:,k))
      figure(5)
      plot(time,y(:,4),ts,z(:,4))
     figure(6)
     plot(time, y(:,1), 'red', time, y(:,3),'blue')