function dy = toggle_switch(t,y,c,alpha2,beta,gamma,reescaleJin,maxBeta)

y(2)= y(2)*maxBeta;

dy = [reescaleJin;
reescaleJin* (c(2)+2*c(3).*t+3*c(4)*t.^2)/maxBeta;
reescaleJin* (y(2)/(1+y(4)^beta)-y(3));
reescaleJin* (alpha2/(1+y(3)^gamma)-y(4));
];


