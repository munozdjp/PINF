function dy = saddlenodeTimeReescale(t,y,c,xMaxstate,reescaletime,maxBeta)
% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz

dy = [reescaletime;
  %-c(1)- 2*c(2)* (-y(1)+c(4)) - 3*c(3)*(-y(1)+c(4))^2;
 reescaletime*(c(1)+2*c(2).*(t)+3.*c(3).*t.^2)/maxBeta;
+reescaletime*(((y(2)*maxBeta)-(y(3)*xMaxstate(3)).^2))/xMaxstate(3)
];
