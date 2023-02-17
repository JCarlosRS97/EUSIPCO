clear; close all; clc;

Nt = 8;
Nr = 4;
nonDIR = 1;
Nris = 20^2;
treshold = 0.95;
L = 100;

if nonDIR
    file_name = sprintf('nonDIR_Nris_%d',Nris);
else
    file_name = sprintf('DIR_Nris_%d',Nris);
end
load(file_name);

% PGM complexity
Cpgm_max = Cpgm(end);
iter_pgm = min(find(treshold*Cpgm_max<Cpgm))-1

Comp_pgm_iter = 2*Nris*Nt*Nr+2*Nt^2*Nr+(3/2)*Nt*Nr^2+Nr^3+...
    Nr*Nris+Nt*Nris+3*Nris+(3/2)*Nt^3
Comp_pgm = Comp_pgm_iter*iter_pgm

% AO complexity
D = min(Nt,Nr);
Cao_max = Cao(end);
iter_ao = min(find(treshold*Cao_max<Cao))-1
Ioi = floor((iter_ao-1)/(Nris+1))+1
Comp_ao = (L+1)*Nr*Nt*Nris+L*(D^3+(1/2)*Nt^2*D)+...
        Ioi*(Nt^3+Nt^2*Nris+2*Nr*Nt*Nris+...
        (2*Nr^2*Nt+2*Nr^3)*Nris+D^3+(1/2)*Nt^2*D)
    


