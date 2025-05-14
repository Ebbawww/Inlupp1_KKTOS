%% Inlämningsuppgift 2
clear all;
clc;

%Konstanter
rho_UO2=10.5*10^(-3); % Bränslets densitet [kg/cm^3]
d_stav=1; %Bränslestavsdiameter[cm]
l=0.0001; %De prompta neutronernas medellivslängd [s]
c_p=0.4*10^(3); %Bränslets specifika värmekapacitet [J/kg*K]
m_UO2 = (d_stav/2)^2*pi*1*rho_UO2; %Bränslets massa. volym * 1 cm bränslestav. [kg]
V_mV_U=2.8; %Volymföhållandet mellan moderator och bränsle
P_bstav=150; %Effekt som utvecklas i 1 cm bränslestav [W/cm]
U_H20=14.3; %Värmeövergångstalet mellan kapsling och vätska [W/(cmK)] 
delta_h_steam = 1506*10^3; %[J/kg]

V_mod=2.8*m_UO2/rho_UO2; %Moderatorvolym [cm^3]
P_mod=150; %Effekt som leds bort av moderatorkylningen [W/cm]

dk_dT= -4*10^(-5); %Förändringen i k med hänseende på temperaturen [K^(-1)]
dk_dalpha=(-100)*10^(-5); %Förändringen i k med hänseende på voiden [% void]

% Initialvärden för t=0
k0 = 1.02;
alpha0 = 0.25;
P0 = 150;
T0 = 1000+273;    

%Tidssteg
t = linspace(0, 1, 100000); %tidsvektorn
x=1/100000; %[s] tidssteg mellan varje iteration

%Bränslestavar faller ur, reaktorn blir prompt kritisk. Loopen simulerar
%effekterna detta har på olika variabler i reaktorn.
for i=1:length(t)
    if i==1 
        %simulerar det första tidsstegen med initialvärdena
        P(i) = P0;
        T(i)=T0;
        k(i)=k0;
        alpha(i)=alpha0;
        P_kapsling(i)=P0;

    else
        %simulerar de övriga tidsstegen
        P(i) = P(i-1) * exp(((k(i-1)-1) / k(i-1)) *x/l); %effekt i reaktorn [W/cm]
        dT_dt(i) = (P(i)-P_kapsling(i-1)) / (m_UO2 * c_p); %hur bränsletemperatur ändras över tid
        T(i) = T(i-1) + (dT_dt(i)*x); %den nya bränsletemperaturen [K]

        rho_steam(i) = XSteam('rho_pT', 70, (T(i)-273))*10^(-6); %ångans densitet för olika T, tryck konstant [kg/cm^3]
        P_kapsling(i)=(T(i-1)-(700+273)-(20+273)) * U_H20 * (1-alpha(i-1)); %effekt från kapslingen. (Temperaturen vid kapsling - kylvattnets T) [W/cm]

        dalpha_dt(i)=(((P_kapsling(i)-P_mod)/delta_h_steam)/rho_steam(i))/V_mod; %förändring i void under tid
        alpha(i) = alpha(i-1) + (dalpha_dt(i)*x); %den nya voiden

        dk_dt(i) = (dk_dalpha * dalpha_dt(i)*100)+(dk_dT * dT_dt(i)); %förändring i k-värde över tid
        k(i) = k(i-1) + (dk_dt(i)*x); %det nya k-värdet.
    end
end

%plot för T, P, void och k.
subplot(4,1,1)
plot(t, T)
title('Temperatur')
xlabel('Tid [s]')
ylabel('Temperatur [K]')

subplot(4,1,2)
plot(t,P)
title('Effekt')
xlabel('Tid [s]')
ylabel('Effekt [W/cm]')

subplot(4,1,3)
plot(t,alpha)
title('Void')
xlabel('Tid [s]')
ylabel('Void [%]')

subplot(4,1,4)
plot(t,k)
title('k')
xlabel('Tid [s]')
ylabel('k')
