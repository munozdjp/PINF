function dy = vanderpol(t,y,c)

dy = [1;
c(2)+2*c(3)*y(1)+3*c(4)*y(1)^2;
y(4)
y(2)*(1-y(3)^2)*y(4)-y(3);
];
