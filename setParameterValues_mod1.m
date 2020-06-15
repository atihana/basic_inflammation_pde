function p=setParameterValues_mod1(phiIn,x,maxT)

% Damage function
mu=0.5;
sigma=0.1;
p.tmax=maxT;
H=@(t)(t<p.tmax/3);
p.dam=@(t)((H(t)/(sigma*sqrt(2*pi)))*exp(-(x'-mu).^2/(2*sigma^2)));

%pro-inflammatory mediators (c) parameters
p.alpha=0.05;                     %production rate [/]
p.nu=0.7;                      % apoptosis rate [/]

%apoptotic neutrophils (a) parameters
p.gamma_a=1;                    %necrosis rate [/]
p.phi=phiIn;%0.001;                       %removal rate [/]
p.beta_a=0.1;                    %saturation constant [/]  

%macrophages (m) parameters
p.gamma_m=0.01;               %leaving tissue rate [/]

% Spatial parameters
p.Dc=1e-4;
p.Dn=1e-5;
p.Dm=1e-6;

end