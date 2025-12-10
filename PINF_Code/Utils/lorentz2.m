function dy = lorentz2(t,y,c,omega,rho,beta,traslation)
%function dy = hopfPolyOrder3(t,y,c,omega,A,traslation)
% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz
%y(3)=y(3)-100;
%y=y-traslation;
%y(3)=y(3);%-100;


%y(3)=y(3)-traslation;
dy = [1;
c(2)+2*c(3)*y(1)+3*c(4)*y(1)^2;
%0;
omega*(y(4)-y(3));
y(3)*(y(2)-y(5))-y(4);
y(3)*y(4)-beta*y(5)
];

