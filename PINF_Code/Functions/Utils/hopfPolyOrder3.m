function dy = hopfPolyOrder3(t,y,c,omega,A,traslation,xMax,reescaletime,maxBeta)

y(2) = y(2)*maxBeta;
y(3) = y(3)*xMax(3);
y(4) = y(4)*xMax(4);

dy = [reescaletime;
reescaletime*(-c(1)-2*c(2)*(-t+c(4))-3*c(3)*(-t+c(4))^2)/maxBeta;
reescaletime*(y(2)*(y(3)) - omega*(y(4)) - A*(y(3))*((y(3))^2+(y(4))^2))/xMax(3);
reescaletime*(omega*(y(3)) + y(2)*(y(4)) - A*(y(4))*((y(3))^2+(y(4))^2))/xMax(4);
];
% 
% dy = [reescaletime;
% reescaletime*(-c(1)-2*c(2)*(-y(1)+c(4))-3*c(3)*(-y(1)+c(4))^2)/maxBeta;
% reescaletime*(y(2)*(y(3)) - omega*(y(4)) - A*(y(3))*((y(3))^2+(y(4))^2))/xMax(3);
% reescaletime*(omega*(y(3)) + y(2)*(y(4)) - A*(y(4))*((y(3))^2+(y(4))^2))/xMax(4);
% ];