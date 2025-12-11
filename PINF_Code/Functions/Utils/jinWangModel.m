function dy = jinWangModel(t,y,c,a1,a2,k1,k2,s,n,reescaleJin,maxBeta )
%function dy = hopfPolyOrder3(t,y,c,omega,A,traslation)
% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz
%y(3)=y(3)-100;
%y=y-traslation;
%y(3)=y(3);%-100;
dy = [reescaleJin;
%c(1)+2*c(2)*(y(1))+3*c(3)*y(1).^2;
reescaleJin*(c(1)+2*c(2).*(t)+3.*c(3).*t.^2)/maxBeta;
reescaleJin*((a1*(y(3))^n)/(s^n+(y(3))^n) + ((y(2)*maxBeta)*s^n)/(s^n+(y(4))^n)-k1*(y(3)));
reescaleJin*((a2*(y(4))^n)/(s^n+(y(4))^n) + ((y(2)*maxBeta)*s^n)/(s^n+(y(3)^n))-k2*(y(4)));
];