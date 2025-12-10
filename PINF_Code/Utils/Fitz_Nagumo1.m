function dy = Fitz_Nagumo1(t,y,c,a,b)

dy = [1;
c(2)+2*c(3)*y(1)+3*c(4)*y(1)^2;
y(3)+y(2)*y(4)-1/(3*y(2)^2)*y(3)^3;
-1/y(2) *(y(3)/y(2)+b*y(4)-a);
];
