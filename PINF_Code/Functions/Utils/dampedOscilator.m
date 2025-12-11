function dy = dampedOscilator(t,y,b,m,k,reescaletime)

dy =[reescaletime;
    reescaletime * y(3);
    reescaletime * (-(b/m)*y(3) - (k/m)*y(2));
    ];