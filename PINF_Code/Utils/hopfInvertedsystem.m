function dy = hopfInvertedsystem(t,y,c,omega,A,xMaxstate,reescaletime,maxBeta)
%y(3)=y(3)-traslation;
%Here the normalization is on doned previously. 
y(2) = y(2) * maxBeta;
y(3) = y(3) * xMaxstate(3);
y(4) = y(4) * xMaxstate(4);

dy = [reescaletime;

reescaletime * ((-1)*(+c(2)+2*c(3)*(t)+3*c(4)*(t)^2))/maxBeta;
reescaletime * ((-(y(2))*y(3) - omega*y(4) - A*y(3)*(y(3)^2+y(4)^2))/xMaxstate(3));
reescaletime * ((omega*y(3) + (-(y(2)))*y(4) - A*y(4)*(y(3)^2+y(4)^2))/xMaxstate(4));
];