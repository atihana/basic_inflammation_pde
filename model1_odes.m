function dydt=model1_odes(t,y,p)

n=y(1);
a=y(2);
m=y(3);
c=y(4);

dydt(1)=c - p.nu*n;
dydt(2)=p.nu*n-p.gamma_a*a-p.phi*m.*a;
dydt(3)=c-p.gamma_m*m;
dydt(4)=p.gamma_a*(a.^2./(p.beta_a^2+a.^2))-c;

dydt=dydt';

end
