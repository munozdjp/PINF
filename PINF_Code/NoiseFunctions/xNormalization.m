function [x,yzero]=xNormalization(x,yzero)
    xMaximum = max(abs(x));
    yzeroMaximum = max(abs(yzero));
    xNormalized = x./xMaximum;
    x = xNormalized;
%     yzeroNormalized = yzero./yzeroMaximum;
%     yzero = yzeroNormalized;
        yzero = yzero;
end 