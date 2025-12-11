function dy = lorentz2(t,y,c,omega,beta,traslation,xMax,reescaletime,maxBeta)
%function dy = hopfPolyOrder3(t,y,c,omega,A,traslation)
% Copyright 2023, All Rights Reserved
% Code by Juan P. Diaz
% For Paper, "IHCV"
% by Juan Diaz, Juan Bernal, Ali B, Jesper T
%y(3)=y(3)-100;
%y=y-traslation;
%y(3)=y(3);%-100;
y(2) = y(2)*maxBeta;
y3 = y(3)*xMax(3);
y4 = y(4)*xMax(4);
y5 = y(5)*xMax(5);

%y(3)=y(3)-traslation;
dy = [reescaletime;
reescaletime*(c(2)+2*c(3)*y(1)+3*c(4)*y(1)^2)/maxBeta;
%0;
reescaletime*(omega*(y4-y3))/xMax(3);
reescaletime*(y3*(y(2)-y5)-y4)/xMax(4);
reescaletime*(y3*y4-beta*y5)/xMax(5)
];

