function Se=stopping(Z1,M1,Z2,M2,E,k_correction)
%% Calculate stopping power
%% project Mat-TRIM (By Yang Yang)
% Inputs:
% Z1--Charge of incident particle [C]
% M1--Atomic Mass of incident particle [u] or [g/mol]
% Z2--Charge of target particle [C]
% M2--Atomic Mass of target particle  [u] or [g/mol]
% E---Incident ion energy [J]


% Basic Physics constant
u      = 1.67*10^-27;        % Atomic mass unit [kg]
me     = 9.10938291*10^(-31);%Electron mass [kg]
hbar   = 6.626*10^-34/2/pi;
e      = 1.6*10^(-19);       %Electron charge [C]
eV     = 1.6*10^(-19);       %Convert keV to Joul unit [J]
keV    = 1.6*10^(-16);       %Convert keV to Joul unit [J]
MeV    = 1.6*10^(-13);       %Convert MeV to Joul unit [J]
Na     = 6.02*10^23;         %Avogardo's constant [#/mol]
a0     = 0.529*10^(-10);     %Bohr Radius
c      = 2.99792458*10^8;    %Speed of light [m/s]
Ai     = 10^(-10);           % length unit [m]
ke     = 8.987551*10^9;      % 1/4/pi/electric_constant [N*m^2/C^2]


% Incident Particle
%M1 = 4 ;  %Atomic Mass of incident particle [u] or [g/mol]
%Z1 = 2 ;  %Charge of incident particle [C]
m1 = M1*u; %Atom mass [kg]


% Target Property
%M2 = 12 ;  %Atomic Mass of target particle  [u] or [g/mol]
m2 = M2*u;  %Atom mass [kg]
%Z2 = 6 ;   %Charge of target particle [C]
rho= 3;     %Target Material Density [g/cm^3]
n = rho*Na/M2*10^6; %atomic density [#/m^3]

% Other Useful Constant
a = 0.8853*a0/(Z1^(1/2)+Z2^(1/2))^(2/3); %Screening Length
V0 = c/137;               %Bohr velocity [m/s]
V_up1 = V0*Z1^(2/3);      %The lower limit 1 for high energy
V_up2 = V0*Z2^(1/2);      %The lower limit 2 for high energy
E_up1 = (M1*u*V_up1^2)/2; %The lower limit 1 for high energy [J]
E_up2 = (M1*u*V_up2^2)/2; %The lower limit 2 for high energy [J]

% Low energy region
k=k_correction*1.212*Z1^(7/6)*Z2/(Z1^(2/3)+Z2^(2/3))^(3/2)/M1^(1/2)*eV^(1/2)*Ai^2;
p=1/2;
S_low=k*E^p;


% High energy region
v1 = (2*E/m1)^(1/2); %velocity of incident particle [m/s]
if Z2<13 
    I0 = 12+7*Z2^(-1);
    I0 = I0*eV;
else
    I0 = 9.76+58.5*Z2^(-1.19);
    I0 = I0*eV;
end
Eb = 2*me*(v1^2)/Z2/I0;
%S_high=ke^2*8*pi*Z1^2*e^4/I0/Eb*log(Eb);

% Middle energy
if Z1<3
    C=100*Z1/Z2;
else
    C=5;
end
S_b  = ke^2*8*pi*Z1^2*e^4/I0/Eb*log(Eb+1+C/Eb);
S_mid= 1/(1/S_low + 1/S_b);

if E<0.5*MeV
    Se = S_low;
else
    Se  = S_mid;
end

end



