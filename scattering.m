function [r0,theta]=scattering(E_reduced,U_ref,P,a)
% Calculte r0 and scattering angle in the CMS frame
% project Mat-TRIM (By Yang Yang)
%-------------------------------------------------------
% Inputs:
% Z1--Charge of incident particle [C]
% M1--Atomic Mass of incident particle [u] or [g/mol]
% Z2--Charge of target particle [C]
% M2--Atomic Mass of target particle  [u] or [g/mol]
% E---Incident ion energy [J]
% E_reduced---Reduced energy
% P---Impact parameter [m]


B=P/a;


Ec=E_reduced*U_ref;


if E_reduced>10
    theta=2*asin((1+(2*E_reduced*B)^2)^(-1/2));
    r0   =0;
else
    R = B;
    
    %% Newton iteration to solve r0
    err=1;
    while err>10^-5
        V=U_ref/R*(0.35*exp(-0.3*R)+0.55*exp(-1.2*R)+0.1*exp(-6*R));
        Vprime=-U_ref/R^2*(0.35*exp(-0.3*R)+0.55*exp(-1.2*R)+0.1*exp(-6*R))-U_ref/R*(0.35*0.3*exp(-0.3*R)+0.55*1.2*exp(-1.2*R)+0.1*6*exp(-6*R));
        f=1-V/Ec-(B/R)^2;
        fprime=1-Vprime/Ec+2*B^2/R^3;
        R_new=R-f/fprime;
        err=abs((R_new-R)/R);
        R=R_new;
    end

    V=U_ref/R*(0.35*exp(-0.3*R)+0.55*exp(-1.2*R)+0.1*exp(-6*R));
    Vprime=(-U_ref/R^2*(0.35*exp(-0.3*R)+0.55*exp(-1.2*R)+0.1*exp(-6*R))-U_ref/R*(0.35*0.3*exp(-0.3*R)+0.55*1.2*exp(-1.2*R)+0.1*6*exp(-6*R)))/a;

    r0=R*a;
    R0=R;
    
    %% Solve Rc
    if V==0
        theta=0; 
    else
        rho_L=2*(Ec-V)/(-Vprime);
        Rc=rho_L/a;

        %% Solve delta
        c1=0.6743;
        c2=0.009611;
        c3=0.005175;
        c4=10;
        c5=6.314;

        alpha=1+c1*E_reduced^(-1/2);
        beta =(c2+E_reduced^(1/2))/(c3+E_reduced^(1/2));
        gamma=(c4+E_reduced)/(c5+E_reduced);

        A = 2*alpha*E_reduced*B^beta;
        G = gamma*((1+A^2)^(1/2)-A)^(-1);

        delta=A*(R0-B)/(1+G);

        %% Solve theta

        theta=2*acos((B+Rc+delta)/(R0+Rc));
    end
end

end