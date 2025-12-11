function dy = pitchforkPolyorder3Left2Right(t,y,c,xMax,reescaletime,maxBeta)

y(2) = y(2)*maxBeta;
y(3) = y(3)*xMax(3);

dy = [reescaletime;
%-c(1)-2*c(2)*(-y(1)+c(4))-3*c(3)*(-y(1)+c(4))^2;
reescaletime*(c(1)+2*c(2).*(y(1))+3.*c(3).*y(1).^2)/maxBeta;
reescaletime*(y(3).*(y(2)-(y(3)).^2))/xMax(3);
];