function dy = lotkavolterrav1(t,y,c,npar,xscale,yscale,betatras,xtras,ytras,scaletimeLV,maxBetaLV)


y(1) = y(1);
y(2) = y(2)*maxBetaLV + betatras;
y(3) = y(3)+xtras;
y(4) = y(4)+ytras;
%reescaling of the variables. 


dy = [scaletimeLV;
%(-1)*(-c(2)-2*c(3)*(-t+c(5))-3*c(4)*(-t+c(5))^2);
scaletimeLV*(c(2)+2*c(3)*(t)+3*c(4)*t.^2)/maxBetaLV;
scaletimeLV*(xscale* (((y(3)/xscale)^2*(1-(y(3)/xscale))/(npar+(y(3)/xscale)))...
 -  (y(3)/xscale)*(y(4)/yscale)));
scaletimeLV*(yscale*(-(y(4)/yscale)*((y(2))-(y(3)/xscale))));
];