%This program is aimed at computing the thermal-hydraulic-mechanical
%behaviors of a saturated porous layer with a thickness of h. A circular
%patch at the surface of the layer is heated by increasing the temperature to T0.
%Written by Zhu Song from Nanjing University of Science and Technology
%Email: songzhu@njust.edu.cn
close all; 
clear;
clc;
syms s xi;
global a alpha alphad betae cd eta G v h kpr kptr ktpr ktr kpz kptz ktpz ktz M md T0; 

alpha=0.119; %From Poroelasticity (Cheng AHD (2016)), Page107, Table 3.2, Rock salt
G=1.24*10^10;
v=0.25; %Poisson's ratio
vu=0.274; %Undrained Poisson's ratio
eta=(1-v)/(1-2*v);
poro=0.001;
betas=1.2*10^-4;%thermal expansion coefficient of skeletal material
%es=4; %ratio set by Reference as betae/betas
betaf=3.0*10^-4;%thermal expansion coefficient of fluid
betad=betas;%(betae-poro*(betaf-betas))/alpha;
betae=1.45*10^-5;%betas*(1-poro)+betaf*poro;%3*es*betas; %
alphad=2.48*10^6;%2*G*(1+v)/3/(1-2*v)*betad;%alphad=K*betad Page 607 (11.56) Cheng AHD
%alphad=3.6*10^4;% From Poroelasticity (Cheng AHD (2016)), Page633, Table 11.3, Berea sandstone
kpr=1.0*10^-18;%kp:permeability, kp=k/miuf, k: intrinsic permeability (1.0*10^-21m^2), miuf: viscosity of fluid (1.0*10^-3N.s/m^2)
kptr=1.0*10^-12;%Thermo-osmosis coefficient,
ktpr=1.0*10^-12;

kpz=1.0*10^-18;%kp:permeability, kp=k/miuf, k: intrinsic permeability (1.0*10^-21m^2), miuf: viscosity of fluid (1.0*10^-3N.s/m^2)
kptz=1.0*10^-12;%Thermo-osmosis coefficient,
ktpz=1.0*10^-12;

cd=1.89*10^6;%betas*(1-poro)+betaf*poro;%the drained specific heat at constant strain

ktr=6.6;  
ktz=6.6;%kpr*(2*G*v/(1-2*v)+2*G)/9800*cd;%2.24; adopted from Reference: kt/cd=2*G(1-v)/(1-2v)/c/Gammaw, 

M=1.81*10^11;%2*G*(vu-v)/alpha^2/(1-2*vu)/(1-2*v);

Tr=293;%the reference temperature, Tr=20"C=293'K

T0=1;%Increment of temperature

md=cd/Tr; %md=cd/T0

a=1;%Radius of heated circular area
h=1;%thickness of the finite layer

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Case: Pore pressure at the centerline 
h0=0.2;%Calculated depth
r0=0;%distance from the centerline
t0=cputime;
delta_h=0.02;
delta_r=0.05;
h1=0;%z=0

TT=[1,2,3,4,5,6,7,8,9,...
    1*10^1,2*10^1,3*10^1,4*10^1,5*10^1,6*10^1,7*10^1,8*10^1,9*10^1,...
    1*10^2,2*10^2,3*10^2,4*10^2,5*10^2,6*10^2,7*10^2,8*10^2,9*10^2,...
    1*10^3,2*10^3,3*10^3,4*10^3,5*10^3,6*10^3,7*10^3,8*10^3,9*10^3,...
    1*10^4,2*10^4,3*10^4,4*10^4,5*10^4,6*10^4,7*10^4,8*10^4,9*10^4,...
    1*10^5,2*10^5,3*10^5,4*10^5,5*10^5,6*10^5,7*10^5,8*10^5,9*10^5,...
    1*10^6,2*10^6,3*10^6,4*10^6,5*10^6,6*10^6,7*10^6,8*10^6,9*10^6,...
    1*10^7,2*10^7,3*10^7,4*10^7,5*10^7,6*10^7,7*10^7,8*10^7,9*10^7,...
    1*10^8,2*10^8,3*10^8,4*10^8,5*10^8,6*10^8,7*10^8,8*10^8,9*10^8,...
    1*10^9,2*10^9,3*10^9,4*10^9,5*10^9,6*10^9,7*10^9,8*10^9,9*10^9];

%%%%%Evolution of pore pressure with time
pc=zeros(1,90);

for i=1:90
    pc(1,i)=quadgk(@(xi) arrayfun(@(xi) talbot_inversion(@(s) ...
        (2.*psi1(s,xi).*m1n(s,xi).*(exp((h0-h).*xi)+exp(-h0.*xi-h.*xi))./(1+exp(-2.*h.*xi))...
        +2.*m3n(s,xi).*(exp((h0-h).*sqrt(lamda1(s,xi)))+exp(-h0.*sqrt(lamda1(s,xi))-h.*sqrt(lamda1(s,xi))))./(1+exp(-2.*h.*sqrt(lamda1(s,xi))))...
        +2.*m4n(s,xi).*(exp((h0-h).*sqrt(lamda2(s,xi)))+exp(-h0.*sqrt(lamda2(s,xi))-h.*sqrt(lamda2(s,xi))))./(1+exp(-2.*h.*sqrt(lamda2(s,xi)))))...
         .*besselj(0,r0.*xi).*xi,TT(i),10),xi), 0,100);
end

pc1=pc(1,1:90);

figure
semilogx(TT,pc1,'b') 
title('Pore pressure')
pc=pc';


