function [sig_eq,to,p] = DPeqstress(sig)
 I2=[1;1;1;0;0;0]; % Kroneker delta in voigt form
  p=(sig(1)+sig(2)+sig(3))/3;
  to=sig-(p.*I2); % Initial deviatoric stress 
 sig_eq=norm(to).*sqrt(1/2);
end
