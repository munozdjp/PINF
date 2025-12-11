function ACDC = ACDC(t,y,c,alpha_x,beta_x, z_x, n_zx, alpha_y, beta_y, x_y, n_xy, delta_y, x_z, n_xz, y_z, n_yz, delta_z,traslation)
%   d_S = (0.008*log(10^-4)/(-4))*x(4); 
%   d_S = (0.008*log(10^(-4))/(-4))*x(4);

%y=y+traslation;
y(3:5)=log(y(3:5));
    dti     = 1;
    d_S     = c(2)+2*c(3)*y(1)+3*c(4)*y(1)^2;
    ACDC_X  = (alpha_x+beta_x*y(2))/(1+y(2)+(y(5)/z_x)^n_zx)-y(3);
    ACDC_Y  = (alpha_y+beta_y*y(2))/(1+y(2)+(y(3)/x_y)^n_xy)-delta_y*y(4);
    ACDC_Z  = 1/(1+(y(3)/x_z)^n_xz+(y(4)/y_z)^n_yz)-delta_z*y(5);
    
    ACDC    = [dti; d_S; ACDC_X; ACDC_Y; ACDC_Z];
end
