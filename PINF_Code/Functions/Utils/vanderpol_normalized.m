function dy = vanderpol_normalized(t, y, c, xMax)
  % y = [y1; y2; u3; u4]  where u3 = y3/xMax(3), u4 = y4/xMax(4)
  % xMax = [xMax1; xMax2; xMax3; xMax4]

  % recover true variables for the 3–4 block
  Y3 = y(3)*xMax(3);
  Y4 = y(4)*xMax(4);

  dy = zeros(4,1);
  % dy1/dt = 1 (unchanged)
  dy(1) = 1;

  % dy2/dt = c2 + 2 c3 y1 + 3 c4 y1^2
  dy(2) = c(2) + 2*c(3)*y(1) + 3*c(4)*y(1)^2;

  % d(u3)/dt = d(y3/xMax3)/dt = (dY3/dt)/xMax3 = Y4 / xMax3
  dy(3) = Y4 / xMax(3);

  % d(u4)/dt = ( dY4/dt ) / xMax4
  %        = [ y2 (1−Y3^2) Y4 – Y3 ] / xMax4
  dy(4) = ( y(2)*(1 - Y3^2)*Y4 - Y3 ) / xMax(4);
end
