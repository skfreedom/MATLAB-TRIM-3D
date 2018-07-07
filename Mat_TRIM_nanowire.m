format long
clear all

%% This is the main function to simulate head-on ion radiation in a nanowire.
%% The example is: 100 keV H in to Au nanowire (radius = 50 nm)


%% Parameters -------------------------------------------------------------
tic
%%------------------------
%  Basic Physics constant
%%------------------------
u    = 1.67*10^-27;        % Atomic mass unit [kg]
me   = 9.10938291*10^(-31);%Electron mass [kg]
hbar = 6.626*10^-34/2/pi;
e    = 1.6*10^(-19);       %Electron charge [C]
eV   = 1.6*10^(-19);       %Convert keV to Joul unit [J]
keV  = 1.6*10^(-16);       %Convert keV to Joul unit [J]
MeV  = 1.6*10^(-13);       %Convert MeV to Joul unit [J]
Na   = 6.02*10^23;         %Avogardo's constant [#/mol]
a0   = 0.529*10^(-10);     %Bohr Radius
c    = 2.99792458*10^8;    %Speed of light [m/s]
Ai   = 10^(-10);           %length unit [m]
ke   = 8.987551*10^9;      %1/4/pi/electric_constant [N*m^2/C^2]

%%------------------------
%  Incident Particle
%%------------------------
M1 = 1 ;    %Atomic Mass of incident particle [u] or [g/mol]
Z1 = 1 ;    %Charge of incident particle [C]
m1 = M1*u;  %Atom mass [kg]

%%------------------------
%  Target(Slab) Property 
%%------------------------
M2 = 197 ;           %Atomic Mass of target particle  [u] or [g/mol]
m2 = M2*u;          %Atom mass [kg]
Z2 = 79 ;           %Charge of target particle [C]
rho= 19.3;          %Target Material Density [g/cm^3]
n  = rho*Na/M2*10^6;%atomic density [#/m^3]
k_correction=1.38;  %Correction factor for stoppingion=1.59;  %Correction factor for stopping 
rmax   = 50*10^-9;  %Radisu of nanowire [m]
rmax_2  = (rmax)^2;%square of Radisu of nanowire [m]
E_threshold=32*eV; % Threshold displacement energy
z_max = 6000*Ai;    %length of the nanowire

%%------------------------
%  Other Useful Constant
%%------------------------
a     = 0.8853*a0/(Z1^(1/2)+Z2^(1/2))^(2/3); %Screening Length
V0    = c/137;                               %Bohr velocity [m/s]
V_up1 = V0*Z1^(2/3);            %The lower limit 1 for high energy
V_up2 = V0*Z2^(1/2);            %The lower limit 2 for high energy
E_up1 = (M1*u*V_up1^2)/2;       %The lower limit 1 for high energy [J]
E_up2 = (M1*u*V_up2^2)/2;       %The lower limit 2 for high energy [J]
U_ref = ke*Z1*Z2*e^2/a;         %reference potential energy for dimensionless E


%% Monte Carlo Simulation -------------------------------------------------

%%------------------------
%  Initialization
%%------------------------
N = 1*10^4;% # of histories

%%------------------------
%  Tally bins
%%------------------------
Emax    = 1*MeV;
Emin    = 10*eV;%particle will die if energy is below Emin
bin_mesh= 20;
bin_size= 6000*10^(-10)/bin_mesh;%The size of lethargy bins
bin     = zeros(bin_mesh,1);%Tally
pka_density = bin;
pka_energy  = bin;  
E_elec      = bin;
E_coll      = bin;
side_out = 0;%Tally the particle exit from left
left_out= 0;%Tally the particle exit from left end of nanowire
right_out= 0;%Tally the particle exit from right end of nanowire
range_tally = 0;
E_out_total = 0; %


   
for i=1:N
    %%------------------------
    %  A new particle is born
    %%------------------------
    life = 1; %life of particle

    %%------------------------
    %  Position initialization
    %%------------------------
    angle    = 2*pi*rand;
    r_initial= rmax*rand;
    x = r_initial*cos(angle);
    y = r_initial*sin(angle);
    z = 0;
    id_old = 1;
    id_new =1;
    r_2 = 0;%lateral radius
    tr_x=0;
    tr_z=0;

    %%------------------------
    %  Energy initialization
    %%------------------------
    E = 100*keV; %Incident ion energy [J]

    %%------------------------
    %  Sample direction
    %%------------------------
    u = 0;
    v = 0;
    w = 1;
    theta=acos(w);

    bug=0;
    while life==1
        bug=bug+1;
        
        id_old=id_new;
        %%------------------------
        % Calculate Reduced energy 
        %%------------------------
        Ec        = E/(1+M1/M2);
        E_reduced = Ec/U_ref; %Dimensionless reduced energy

        %%------------------------
        % Calculate Stopping 
        %%------------------------
        Se = stopping(Z1,M1,Z2,M2,E,k_correction);%Stopping [J*m^2]
        distance_5percent=0.05*E/Se/n;        

        %%------------------------
        % Sample Pathinterval and Impact parameter
        %%------------------------
        if E_reduced>10
            distance = 0.02*(1+(M1/M2))^2/4/pi/a^2/n*(E_reduced^2+0.052*E_reduced^1.32)/log(1+E_reduced);
            if distance>distance_5percent
                distance=distance_5percent;
            end
            P=(-log(rand)/pi/n/distance)^(1/2);
         elseif E_reduced>0.03
            distance = n^(-1/3);
            P=(rand/(pi*n^(2/3)))^(1/2);

        else
            P=(rand/(pi*n^(2/3)))^(1/2);
            distance = n^(-1/3)-P*tan(theta/2);
        end
        
        %%------------------------
        % Correct the distance if it can leak out from the side
        %%------------------------
        C_k = x*u+y*v;
        C_a = u^2+v^2;
        C_c = x^2+y^2-rmax^2;
        d_min = (-C_k+(C_k^2-C_a*C_c)^0.5)/C_a;
        if distance>d_min
           distance =   d_min;
           life = 0; %kill thee particle
           side_out      = side_out+1;
           
        end
           
        
        
        %%------------------------
        %  Move to new location
        %%------------------------  
        x      = x+u*distance;
        y      = y+v*distance;
        dz     = w*distance;
        z_old  = z; %save the old z
        z      = z+dz;
        
        r_2    = x^2+y^2;
        id_new = ceil(z/bin_size);
        tr_x(bug+1,1)=x;
        tr_z(bug+1,1)=z;
        
        E_loss = distance*n*Se;
        
        %% add energy loss tally (Important)
        if id_new == id_old
            E_elec(id_new,1)=E_elec(id_new,1)+E_loss;
        else
            if dz>=0 
                if z<=z_max
                    if dz>=2*bin_size
                        E_elec(id_new,1)=E_elec(id_new,1)+E_loss*(z-(id_new-1)*bin_size)/dz;
                        E_elec(id_old,1)=E_elec(id_old,1)+E_loss*(id_old*bin_size-z_old)/dz;
                        for tt=(id_old+1):(id_new-1)
                            E_elec(tt,1)=E_elec(tt,1)+E_loss*(bin_size)/dz;
                        end
                    else 
                        E_elec(id_new,1)=E_elec(id_new,1)+E_loss*(z-(id_new-1)*bin_size)/dz;
                        E_elec(id_old,1)=E_elec(id_old,1)+E_loss*(id_old*bin_size-z_old)/dz;
                    end
                else 
                        E_elec(id_old,1)=E_elec(id_old,1)+E_loss*(id_old*bin_size-z_old)/dz;
                        for tt=(id_old+1):(bin_mesh)
                            E_elec(tt,1)=E_elec(tt,1)+E_loss*(bin_size)/dz;
                        end
                        life = 0; %kill thee particle
                        right_out      = right_out+1;
                end
            else
                if z>=0
                    if dz<=-2*bin_size
                        E_elec(id_new,1)=E_elec(id_new,1)+E_loss*(id_new*bin_size-z)/dz;
                        E_elec(id_old,1)=E_elec(id_old,1)+E_loss*(z_old-(id_old-1)*bin_size)/dz;
                        for tt=(id_new+1):(id_old-1)
                            E_elec(tt,1)=E_elec(tt,1)+E_loss*(bin_size)/dz;
                        end
                    else
                        E_elec(id_new,1)=E_elec(id_new,1)+E_loss*(id_new*bin_size-z)/dz;
                        E_elec(id_old,1)=E_elec(id_old,1)+E_loss*(z_old-(id_old-1)*bin_size)/dz;
                    end
                else
                        E_elec(id_old,1)=E_elec(id_old,1)+E_loss*(z_old-(id_old-1)*bin_size)/dz;
                        for tt=1:(id_old-1)
                            E_elec(tt,1)=E_elec(tt,1)+E_loss*(bin_size)/dz;
                        end
                        life = 0; %kill thee particle
                        left_out      = left_out+1;
                end
            end
        end
        
        

                    
%         if dz>=0
%             
%             z>6000*10^(-10)
%             E_elec(bin_mesh,1)=E_elec(bin_mesh,1)+E_loss;
%         elseif abs(dz)>2*bin_size
%             
%             
%             E_elec(id_new,1)=E_elec(id_new,1)+E_loss;
%         end
        E      = E-E_loss;
        Ec        = E/(1+M1/M2);
        E_reduced = Ec/U_ref; %Dimensionless reduced energy

        % Leak from left
%         if z<0 || r_2>=rmax_2
%             life = 0; %kill thee particle
%             side_out      = side_out+1;
%             %E_left(left_out,1)= E; %Memorize the exit particle energy


        % collision
        if life==0
            E_out_total   = E_out_total+E;
        else
            %%------------------------
            %% sample direction
            %%------------------------
            [r0,theta] = scattering(E_reduced,U_ref,P,a);
            mu_c = cos(theta);

            phi   = rand*2*pi; 
            if mu_c==-1
                theta_L=pi/2;
            else

                theta_L=atan(sin(theta)/(cos(theta)+(M1/M2)));
            end
            mu_L  = cos(theta_L);
            S_mu_L= sqrt(1-mu_L^2);
            C_phi = cos(phi);
            S_phi = sin(phi);
            s_w   = sqrt(1-w^2); 

            if bug==1
                u     = S_mu_L*C_phi;
                v     = S_mu_L*S_phi;               
            elseif s_w~=0
                u_new = mu_L*u+S_mu_L*(u*w*C_phi-v*S_phi)/s_w;
                v_new = mu_L*v+S_mu_L*(v*w*C_phi+u*S_phi)/s_w;            
                u     = u_new;
                v     = v_new;
            end
            w = mu_L*w-S_mu_L*s_w*cos(phi);%new w
            
            
            %%------------------------
            %% sample energy
            %%------------------------
            T = 4*M1*M2/((M1+M2)^2)*E*sin(theta/2)^2;
            if z>6000*10^(-10)
                E_coll(bin_mesh,1)=E_coll(bin_mesh,1)+T;
            elseif z>0
                E_coll(id_new,1)=E_coll(id_new,1)+T;
            end
            if T>E_threshold
                    if z>6000*10^(-10)
                        pka_density(bin_mesh,1)=pka_density(bin_mesh,1)+1;
                        pka_energy(bin_mesh,1)=pka_energy(bin_mesh,1)+T-E_threshold;
                    elseif z>0
                        pka_density(id_new,1)=pka_density(id_new,1)+1;
                        pka_energy(id_new,1)=pka_energy(id_new,1)+T-E_threshold;
                    end
            end
            E = E-T;

            %%------------------------
            %% kill particle if E is too small
            %%------------------------           
            if E<Emin
                life = 0;
            end           

        end
        

    end
    
%plot(tr_z*10^10,tr_x*10^10,'b-');drawnow;
    

if r_2<rmax_2
    if z>6000*10^(-10)
        bin(bin_mesh,1)=bin(bin_mesh,1)+1;
    elseif z>0
        bin(id_new,1)=bin(id_new,1)+1;
    end
end

    range_tally=range_tally+z;

end

range = range_tally/N;

pka_ave_energy=pka_energy./pka_density;
%E_out_total/side_out/eV/1000;
toc













