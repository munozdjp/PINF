function dy = hopfOscilatingModel(t,y,c,omega,A,traslation)
%function dy = hopfPolyOrder3(t,y,c,omega,A,traslation)
% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz
%y(3)=y(3)-100;
%y=y-traslation;
%y(3)=y(3);%-100;
y(3)=y(3)-traslation;
dy = [1;
-c(1)-2*c(2)*(-y(1)+c(4))-3*c(3)*(-y(1)+c(4))^2+c(7)*c(5)*cos(y(1)*c(7));
y(2)*y(3) - omega*y(4) - A*y(3)*(y(3)^2+y(4)^2);
omega*y(3) + y(2)*y(4) - A*y(4)*(y(3)^2+y(4)^2);
];