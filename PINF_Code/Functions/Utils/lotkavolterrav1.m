function dy = lotkavolterrav1(t, y, c, npar, xMax, scaletimeLV, maxBetaLV)
% lotkavolterrav1_norm  Normalized‐entry / raw‐exit Lotka–Volterra variant
%
%   y = [ τ; β_norm; x3_raw; x4_raw ]
%     – τ is just t carried along
%     – β_norm is β scaled so that β = β_norm*maxBetaLV
%     – x3_raw, x4_raw are the “true” state variables
%
%   c        : 1×5 polynomial coefficients [c1,c2,c3,c4,Initinterval]
%   npar     : nonlinear parameter
%   xMax     : [1,1,xscale,yscale]
%   scaletimeLV : time‐scaling (for τ̇)
%   maxBetaLV   : scale for β

  % unpack & renormalize inputs
  beta = y(2)*maxBetaLV;      % recover true β
  U3   = y(3)/xMax(3);        % normalized x3
  U4   = y(4)/xMax(4);        % normalized x4

  % precompute the raw‐derivative factors
  D3 = U3^2*(1 - U3)/(npar + U3) - U3*U4;
  D4 = -U4*(beta - U3);

  dy = zeros(4,1);
  dy(1) = scaletimeLV;                                         % τ̇ = scaletimeLV
  dy(2) = scaletimeLV*(c(2) + 2*c(3)*t + 3*c(4)*t^2)/maxBetaLV; % β̇
  
  % multiply back by the scale to get raw‐derivatives of x3,x4
  dy(3) = scaletimeLV * xMax(3) * D3;  
  dy(4) = scaletimeLV * xMax(4) * D4;
end
% function dy = lotkavolterrav1(t,y,c,npar,xscale,yscale,betatras,xtras,ytras,scaletimeLV,maxBetaLV)
% 
% 
% y(1) = y(1);
% y(2) = y(2)*maxBetaLV + betatras;
% y(3) = y(3)+xtras;
% y(4) = y(4)+ytras;
% %reescaling of the variables. 
% 
% 
% dy = [scaletimeLV;
% %(-1)*(-c(2)-2*c(3)*(-t+c(5))-3*c(4)*(-t+c(5))^2);
% scaletimeLV*(c(2)+2*c(3)*(t)+3*c(4)*t.^2)/maxBetaLV;
% scaletimeLV*(xscale* (((y(3)/xscale)^2*(1-(y(3)/xscale))/(npar+(y(3)/xscale)))...
%  -  (y(3)/xscale)*(y(4)/yscale)));
% scaletimeLV*(yscale*(-(y(4)/yscale)*((y(2))-(y(3)/xscale))));
% ];