function [n1,a1,m1,c1,x,t]=solve1D(phi)

%clear all
%close all
clc
%profile on
nX=100;
x=linspace(0,1,nX);              

dx = x(2)-x(1);         
L=length(x);         

tmax=10000;            

BC='Periodic';

% Set parameter values
s=setParameterValues_mod1(phi,x,tmax);
s.ni=s.nu;

%=========================================
% FOR ICs AS IN CHAPTER 3
% Find corresponding steady states
% y0=[3.9406,0.0770,38.2178,0.3822]';
% ICs=fsolve(@(y)model1_odes(0,y,s,0),y0);
% n0=ICs(1);
% a0=ICs(2);
% m0=ICs(3);
% c0=ICs(4);
% damageRadius=0.25;
% n=zeros(size(x));a=n;m=n;c=n;
% n(((x-0.5).^2<damageRadius^2))=n0;
% a(((x-0.5).^2<damageRadius^2))=a0;
% m(((x-0.5).^2<damageRadius^2))=m0;
% c(((x-0.5).^2<damageRadius^2))=c0;
% s.dam=@(t)(0.*x');

%=========================================
% FOR ZERO ICs (WITH DAMAGE FUNCTION)
n=0.*x;
a=0.*x;
m=0.*x;
c=0.*x;
%=========================================

%boundary conditions
%diffusion
s.dmat=sparse(getDiffMatrix(L,dx,BC));

s.left=[L, 1:L-1];
s.right=[2:L, 1];

%initial conditions for ode15s
v=[n,a,m,c]';            
tspan=[0 tmax];
s.L=L;

s.rn=s.Dn/(dx*dx);
s.rm=s.Dm/(dx*dx);
s.rc=s.Dc/(dx*dx);

%options=odeset('Vectorized','off','JPattern',createJacobianMatrixChemo(L),'RelTol',1e-3,'AbsTol',1e-4);
options=odeset('Vectorized','off','RelTol',1e-3,'AbsTol',1e-4);
[t,dv]=ode15s(@(t,v)mod1(t,v,s),tspan,v,options);

n1=dv(:,1:L);
a1=dv(:,L+1:2*L);
m1=dv(:,2*L+1:3*L);
c1=dv(:,3*L+1:4*L);

% Windows...
%save(strcat('E:\Anahita\DataFiles\ParamSet2\',fname),'n1','a1','m1','c1','g1','t','Lx');
% Linux...
%save(strcat('output/',fname),'n1','a1','m1','c1','g1','t','Lx');
%profile viewer
end
