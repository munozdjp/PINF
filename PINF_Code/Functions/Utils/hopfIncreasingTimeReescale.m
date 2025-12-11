function dy = hopfIncreasingTimeReescale(t,y,c,xMaxstate,reescaletime,maxBeta)

y(2) = y(2) * maxBeta;
y(3) = y(3) * xMaxstate(3);
y(4) = y(4) * xMaxstate(4);

dy = [1;
reescaletime * (-c(1)-2*c(2)*(-y(1)+c(4))-3*c(3)*(-y(1)+c(4))^2 )/maxBeta;;    
%reescaletime * ((+c(2)+2*c(3)*(t)+3*c(4)*(t)^2))/maxBeta;
reescaletime * ((y(2)*(y(3)) - (y(4)) - (y(3))*((y(3))^2+(y(4))^2))/xMaxstate(3));
reescaletime * (((y(3)) + y(2)*(y(4)) - (y(4))*((y(3))^2+(y(4))^2))/xMaxstate(4));
];