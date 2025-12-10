function dy = Fitz_Nagumo2th(t,y,c,...
    a,lambda,epsilon)

dy = [1;
c(2)+2*c(3)*y(1)+3*c(4)*y(1)^2;
epsilon*(y(3)*(y(3)-lambda)*(1-y(3)))-y(4)+y(2);
y(3)-a*y(4);
];
