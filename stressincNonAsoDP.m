function [dsig,depsp,dqvec,ddsdde,dsighydro,dmeanepsp,dseq] = stressincNonAsoDP(sig,qvec,deps,props)

 %Material parameters
  K= props(1);
  G = props(2);
  sigma_zero=props(3);
  H=props(4);
  a=props(5);
  b=props(6);
  lamb=props(7);
  mu=G;
  s=qvec;
  
  iter_counter=0;
  nint = length(qvec);
  dqvec = zeros(nint,1);
  depsp = zeros(6,1);
  I2=[1;1;1;0;0;0]; % Kroneker delta in voigt form
  P3=eye(6)-(1/3).*I2*I2';
  if (a==0&&b==0)
      stif = 2*mu*eye(6)+lamb*(I2*I2');
  else
      stif = (2*G*eye(6))+((K-(2*G/3))*(I2*I2')); % Elastic stiffness tensor in Voigt form
  end
  
  dsig_t=stif*deps; % change in trial stress
  sig_t=sig+dsig_t; % trial-stress
  [sig_eq,to,p]=DPeqstress(sig_t);
  Na=((to./sig_eq).*1/2)+(a*I2);
  Qb=((to./sig_eq).*1/2)+(b*I2);
  h=G+(9*K*a*b)+(H*(a+(1/sqrt(3)))*sqrt(((1/3)+(2*b*b))));
  sigmay=sigma_zero+(H*s);     % initial yield stress
  dlambda=0;
 f_trial=abs(sig_eq)+(3*a*p)-(a+sqrt(1/3))*sigmay;
 if f_trial <=1e-1
     sig=sig_t;
     dsig=dsig_t;
     dsighydro=p;
     dseq=sig_eq;
     dmeanepsp=0;
 else
     k=1;
       if (f_trial>1e-2&&iter_counter<100)
          x0=1e-5;
          Na=((to./sig_eq).*1/2)+(a*I2);
          Qb=((to./sig_eq).*1/2)+(b*I2);
          sigmay=sigma_zero+(H*s); 
          ddl=newton_rp(sig_t,Qb,stif,sigmay,x0,a)
%           ddl=NR(sig_t,stif,Qb,a,sigmay)
%           f_trial=@(ddl) ddlsolver(sig_t,stif,Qb,a,ddl,sigmay)
%           ddl=fzero(f_trial,1e-5)
%           dlambda=dlambda+ddl;
          ddepsp=ddl.*Qb;
          sig_t=sig_t-stif*ddepsp;
          [sig_eq,to,p]=DPeqstress(sig_t);
          h=G+(9*K*a*b)+(H*(a+(1/sqrt(3)))*sqrt(((1/3)+(2*b*b))));
          f_trial=abs(sig_eq)+(3*a*p)-(a+sqrt(1/3))*sigmay;
          depsp=depsp+ddepsp;
          dsighydro=p;
          dseq=sig_eq;
          dmeanepsp=(sqrt(2/3))*sqrt(ddepsp'*ddepsp);
          s=s+ddl;
          dqvec=s;
          iter_counter=iter_counter+1;
      end
      dsig=sig_t-sig;
      sig=sig_t;
      stif=stif-((1/h)*(stif'*(Qb*Na')*stif));   
 end   
  ddsdde=stif;
end





    
    
    function [f] = yld(sig_t,Nij,stif,sigmay,x,a)
    sig = sig_t - x.*(stif*Nij);
    [sig_eq,to,p] = DPeqstress(sig);
    f=abs(sig_eq)+(3*a*p)-(a+sqrt(1/3))*sigmay;
    end
    
    function [xn] = newton_rp(sig_t,Nij,stif,sigmay,x0,a)
    tol = 1e-6;
    x1 = x0+0.0001;
    f0 = yld(sig_t,Nij,stif,sigmay,x0,a);
    f1 = yld(sig_t,Nij,stif,sigmay,x1,a);
    df = (f0-f1)/(x0-x1);
    xn = x0 - (f0/df);
    f = yld(sig_t,Nij,stif,sigmay,xn,a);
    i = 1;
    maxiter = 100;
    
    while (f > tol)
        x0 = xn;
        x1 = x0 + 0.001;
        f0 = yld(sig_t,Nij,stif,sigmay,x0,a);
        f1 = yld(sig_t,Nij,stif,sigmay,x1,a);
        df = (f0-f1)/(x0-x1);
        xn = x0 - (f0/df);
        f = yld(sig_t,Nij,stif,sigmay,xn,a);
        i = i + 1;
        if ( i > maxiter-1)
           disp('solution not converged in ddlambda');
            break;
        end
    end
    end
    
%     function [f_trial]=ddlsolver(sig_t,stif,Qb,a,ddl,sigmay)
%     I2=[1;1;1;0;0;0];
%     sig = sig_t - ddl.*(stif*Qb);
%     p=(sig(1)+sig(2)+sig(3))/3;
%     to=sig-(p.*I2); 
%     sig_eq=norm(to).*sqrt(1/2);
%     f_trial=abs(sig_eq)+(3*a*p)-sigmay;
%     end

% function [y]=NR(sig_t,stif,Qb,a,sigmay)
% syms ddl
% I2=[1;1;1;0;0;0];
% sig = sig_t - ddl.*(stif*Qb);
% p=(sig(1)+sig(2)+sig(3))/3;
% to=sig-(p.*I2); 
% sig_eq=norm(to).*sqrt(1/2);
% f=abs(sig_eq)+(3*a*p)-sigmay;
% g=diff(f); %The Derivative of the Function
% tolerance=1e-5;
% x=0; % intial approximation;
% for i=1:100
%      fun=vpa(subs(f,ddl,x)); %Calculating the value of function at x0
%      fun_der=vpa(subs(g,ddl,x)); %Calculating the value of function derivative at x0
%      y=x-(fun/fun_der); % The Formula
% err=abs(y-x);
% if err<tolerance %checking the amount of error at each iteration
% break
% end
% x=y;
% end
% end