function dy = Fitz_Nagumo2th(t,y,c,a,lambda,epsilon,xMax,reescaletime,maxBeta)
y(2) = y(2)*maxBeta;
y3 = y(3)*xMax(3);
y4 = y(4)*xMax(4);

dy = [reescaletime;
reescaletime*(c(2)+2*c(3)*y(1)+3*c(4)*y(1)^2)/maxBeta;
reescaletime*(epsilon*(y3*(y3-lambda)*(1-y3))-y4+y(2))/xMax(3);
reescaletime*(y3-a*y4)/xMax(4);
];
