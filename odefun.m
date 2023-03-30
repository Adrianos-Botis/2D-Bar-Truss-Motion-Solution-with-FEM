%---Diferrential equation solver-----------------------------
function dydt = odefun(t,y,m11, m12, m21, m22, k11, k12, k21, k22)
 dydt = zeros(4,1);
  if (m11 == 0)
    dydt(1) = 0;
    dydt(2) = 0;
    dydt(3) = y(4);
    dydt(4) = (-k22/m22)*y(3); 
  elseif (m22 == 0)
    dydt(1) = y(2);
    dydt(2) = (-k11/m11)*y(1);
  elseif (m12 == 0) || (m21 == 0)
   dydt(1) = y(2);
   dydt(2) = (-k11/m11)*y(1);
   dydt(3) = y(4);
   dydt(4) = (-k22/m22)*y(3);
  end
      
end