function dy = hopfPolyOrder3(t,y,c,omega,A,traslation,xMax)


dy = [1;
-c(1)-2*c(2)*(-y(1)+c(4))-3*c(3)*(-y(1)+c(4))^2;
(y(2)*(y(3)*xMax(3)) - omega*(y(4)*xMax(4)) - A*(y(3)*xMax(3))*((y(3)*xMax(3))^2+(y(4)*xMax(4))^2))/xMax(3);
(omega*(y(3)*xMax(3)) + y(2)*(y(4)*xMax(4)) - A*(y(4)*xMax(4))*((y(3)*xMax(3))^2+(y(4)*xMax(4))^2))/xMax(4);
];