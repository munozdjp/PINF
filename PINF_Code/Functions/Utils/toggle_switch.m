function dy = toggle_switch(t,y,c,alpha2,beta,gamma,xMax,reescaleJin,maxBeta)

y(2)= y(2)*maxBeta;
y(3) = y(3)*xMax(3);
y(4) = y(4)*xMax(4);

dy = [reescaleJin;
reescaleJin* (c(2)+2*c(3).*t+3*c(4)*t.^2)/maxBeta;
reescaleJin* (y(2)/(1+y(4)^beta)-y(3))/xMax(3);
reescaleJin* (alpha2/(1+y(3)^gamma)-y(4))/xMax(4);
];


