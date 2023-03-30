%----Solver for specific angle Between the Bars----------
function dydt = odefun1(t,y,M,K)
% y = [u, v, du, dv]'
% dydt = [du, dv, d2u, d2v]'
A(1:2,1:2) = zeros(2);
A(1:2,3:4) = eye(2); 
A(3:4,1:2) = -linsolve(M,K);
A(3:4,3:4) = zeros(2);
dydt = A*y;
end