%model 1 space discreitzed system equation (ODEs in time) 

function dv=mod1(t,v,p)

n=v(1:p.L);
a=v(p.L+1:2*p.L);
m=v(2*p.L+1:3*p.L);
c=v(3*p.L+1:4*p.L);


dv(1:p.L)=-p.nu*n+c+p.rn*(n(p.left)-2*n+n(p.right));
dv(p.L+1:2*p.L)=p.nu*n-p.gamma_a*a-p.phi*m.*a;
dv(2*p.L+1:3*p.L)=-p.gamma_m*m+c+p.rm*(m(p.left)-2*m+m(p.right));
dv(3*p.L+1:4*p.L)=-c+p.rc*(c(p.left)-2*c+c(p.right))+p.alpha*p.dam(t)+p.gamma_a*a.^2./(a.^2+p.beta_a^2);

dv=dv';

