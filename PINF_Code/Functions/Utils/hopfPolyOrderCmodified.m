function dy = hopfPolyOrderCmodified(t,y,c,omega,A,traslation,xMax,reescaletime,maxBeta)

y(2) = y(2)*maxBeta;
y(3) = y(3)*xMax(3);
y(4) = y(4)*xMax(4);

dy = [reescaletime;

% Polyorder 11
reescaletime*(-c(2)-2*c(3)*(-t+c(5))-3*c(4)*(-t+c(5))^2 ...
-4*c(6)*(-t+c(5))^3-5*c(7)*(-t+c(5))^4 ...
-6*c(8)*(-t+c(5))^5-7*c(9)*(-t+c(5))^6 ...
-8*c(10)*(-t+c(5))^7-9*c(11)*(-t+c(5))^8 ...
-10*c(12)*(-t+c(5))^9-11*c(13)*(-t+c(5))^10)/maxBeta;


reescaletime*(y(2)*(y(3)) - omega*(y(4)) - A*(y(3))*((y(3))^2+(y(4))^2))/xMax(3);
reescaletime*(omega*(y(3)) + y(2)*(y(4)) - A*(y(4))*((y(3))^2+(y(4))^2))/xMax(4);
];