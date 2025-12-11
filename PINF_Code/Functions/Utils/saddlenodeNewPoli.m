function dy = saddlenodeNewPoli(t,y,c,xMaxstate,reescaletime,maxBeta)

dy = [reescaletime;
  %-c(1)- 2*c(2)* (-y(1)+c(4)) - 3*c(3)*(-y(1)+c(4))^2;
%This is mine normalized time with the constant (maxBeta)
 reescaletime*(c(1)+2*c(2).*(t)+3.*c(3).*t.^2)/maxBeta;
+reescaletime*((y(2)*maxBeta-(y(3)*xMaxstate(3)).^2))/xMaxstate(3)
];


