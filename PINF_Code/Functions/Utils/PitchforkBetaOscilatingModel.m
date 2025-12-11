function dy = PitchforkBetaOscilatingModel(t,y,c)
%here the coeficcient of the polynomial should be zero.
%Nevertheless we can get a 
% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz

dy = [1;
  %-c(1)- 2*c(2)* (-y(1)+c(4)) - 3*c(3)*(-y(1)+c(4))^2;
 %c(1)+2*c(2)*(y(1))+3*c(3)*y(1).^2;
 c(6)*c(4)*cos(y(1)*c(6));
+y(3).*(y(2)-y(3).^2)
];
