%Inlämningsuppgift 1

%%%%%%%%%% Beräkning av konverteringskvoten
%C = (Sigma_a238/Sigma_a235)+(Sigma_f235/Sigma_a235)*epsilon*my_u235*P_s*(1-p);

clc;clear;

n = 1;
while n < 3
    
    %%%%%%%%%% Konstanter (typ)
    t = 60*60*24*365; %[s] tid ett år
    epsilon = 1.06;
    P_leakage = 0.97;
    P_s = 0.96; %icke-läckagefaktor för snabba neutroner - taget från första bästa källa på nätet
    P = 2928*10^6; %Termisk reaktoreffekt
    m_UO2 = 118000; %[kg]vikt för uranoxidet
    u = 1.66*10^-24; %[kg] universala massenheten
    E_fission = 200*10^6; %[MeV] utlöst energi per fission
    ny_pu239 = 2.9; % genomsnittligt antal för hur många neutroner som bildas per fission för Pu239
    ny_u235 = 2.42; % genomsnittligt antal för hur många neutroner som bildas per fission för U235
    sigma_f235 = 582.2*10^-24; %[barn] mikroskopiskt tvärsnitt för fission U235
    sigma_c235 = 93.6*10^-24; % [barn] mikroskopiskt tvärsnitt för infångning U235
    sigma_c238 = 2.7*10^-24; % [barn] mikroskopiskt tvärsnitt för infångning U238
    sigma_a238 = 2.7*10^-24; 
    rho_U = 10.97; % [g/cm3] urandioxids densitet, ändras när uran235 konverterar till plutonium?
    T_u = 1500+273.15; % [K] Temperatur ?, Maxtemp för Forsmark 1 är 1800 C
    N0_u235 = (0.025*rho_U)/((235+16*2)*u); %[kärnor/cm3] kärntäthet U235 vid start
    N0tot_u235 = N0_u235*(m_UO2/(10.97*10^-3)); % antal kärnor U235 totalt i härden vid start
    %N0_u235 = (e*rho_U)/((235+16*2)*u);
    N0_u238 = ((1-0.025)*rho_U)/((238+16*2)*u); %[kärnor/cm3] kärntäthet U238 vid start, KSU s102
    N0tot_u238 = N0_u238*(m_UO2/(10.97*10^-3)); % antal kärnor U238 totalt i härden vid start
    N0_pu239 = 0; %antal Pu239 kärnor vid start
    sigma_a235 = sigma_f235+sigma_c235; %%Tvärsnitt för absorption för U235, konverteras till Pu239 -> Sigma_abransle ändras därav över tid (idé: använd samma timestep som resten, splitta så att den går från 0% till 30% Pu över cykeln. Tillsammans med C?
    sigma_f239 = 744.4*10^-24; %tagen från KSU
    sigma_c239 = 268.8*10^-24; %tagen från KSU
    sigma_a239 = sigma_f239+sigma_c239;
    rho_Xe = 5.9; %[g/cm3] densitet för Xe135
    FR = P/(200*1.6021*(10^-13));%*(m_U02/(rho_U*10^-3))); %[1/cm3/s] fissionsrat för bränslet (?)
    X = linspace(0,10,1000);
    
        %%%%%%%%%% Beräkning av resonanspassagefaktorn
        
    if n == 1
        e = 0.025;
        N_u = ((1-e)*rho_U)/((238+16*2)*u);%antal kärnor U238 [kärnor/(cm3 bränsle)]
    elseif n == 2
        e = (Nt_u235+Nt_pu239) / (Nt_u238+Nt_u235+Nt_pu239); % anrikning efter t tid (ett år)
        N_u = ((1-e)*rho_U)/((238+16*2)*u); %antal kärnor U238 [kärnor/(cm3 bränsle)]
    end
    R_u = 1/2; % [cm] stavradie, data tagen från Forsmark 1
    V_u = 1; %Volym bränsle, sätter Vu/Vm till 1/3, behöver diskuteras vad det egentligen är
    V_m = 3; %Volym moderator
    xi = 1.26; %log. medelenergiförlusten för moderatorn, sätter xi*Sigma_s = 1.28 enligt KSU. Ändra detta till senare? 
    Sigma_s = 1; %makroskopiskt spridningstvärsnitt i moderatormaterialet [cm-1]

    sigma_Iabs_300K = 3 + 39.6/(sqrt(R_u*rho_U));
    B_1 = 6.1*10^(-3)+(0.94*10^(-2)/(R_u*rho_U));
    sigma_IabsT = sigma_Iabs_300K*(1+B_1*(sqrt(T_u)-sqrt(300))); %[barn]
    
    p = exp(-(1-e)*((N_u*sigma_IabsT*(10^-24)*V_u)/(xi*Sigma_s*V_m)));
    
    p_all = zeros(1,length(X));
    for i = 1:1:length(X)
        V_m = X(i);
        p_all(i) = exp(-(1-e)*((N_u*sigma_IabsT*(10^-24)*V_u)/(xi*Sigma_s*V_m)));
    end
    
    if n == 1
        e = 0.025; % anrikning vid start
        Sigma_a235 = ((e*rho_U)/(235*1.66043*10^-24))*sigma_a235; %makroskopiska tvärsnittet för absorption för uran235
        Sigma_a238 = (((1-e)*rho_U)/(238*1.66043*10^-24))*sigma_a238; %makroskopiska tvärsnittet för absorption för uran238
        Sigma_f235 = (((e)*rho_U)/(235*1.66043*10^-24))*sigma_f235; %makroskopiska tvärsnittet för fission för uran235
        %C_0 = (sigma_c238*N0_u238)/(sigma_a235*N0_u235+sigma_a239*N0_pu239); %koverteringskvot vid start
        %disp(['C_0 1 = ' num2str(C_0)]);
        C_0 = (Sigma_a238/Sigma_a235)+(Sigma_f235/Sigma_a235)*epsilon*ny_u235*(1-p);
        disp(['C_0 = ' num2str(C_0)]);
        disp(['e = ' num2str(e)]);
    elseif n == 2
        e = (Nt_u235+Nt_pu239) / (Nt_u238+Nt_u235+Nt_pu239); % anrikning efter t tid (ett år)
        Sigma_a235 = ((e*(1-andel_PuvsU)*rho_U)/(235*1.66043*10^-24))*sigma_a235; %makroskopiska tvärsnittet för absorption för uran235
        Sigma_a238 = (((1-e)*rho_U)/(238*1.66043*10^-24))*sigma_a238; %makroskopiska tvärsnittet för absorption för uran238
        Sigma_f235 = (((e*(1-andel_PuvsU))*rho_U)/(235*1.66043*10^-24))*sigma_f235; %makroskopiska tvärsnittet för fission för uran235
        Sigma_a239 = ((e*andel_PuvsU*rho_U)/(239*1.66043*10^-24))*sigma_a239;
        Sigma_f239 = ((e*andel_PuvsU*rho_U)/(239*1.66043*10^-24))*sigma_f239;
        %C_t = (sigma_c238*Nt_u238)/(sigma_a235*Nt_u235+sigma_a239*Nt_pu239); %konverteringskvot efter t tid (ett år)
        C_t = (Sigma_a238/(Sigma_a235*(1-andel_PuvsU)+Sigma_a239*andel_PuvsU))+((Sigma_f235*(1-andel_PuvsU)+Sigma_f239*andel_PuvsU)/(Sigma_a235*(1-andel_PuvsU)+Sigma_a239*andel_PuvsU))*epsilon*ny_u235*(1-p);
        disp(['C_t = ' num2str(C_t)]);
        disp(['e = ' num2str(e)]);
    end
    
    
    %%%%%%%%%% Beräkning av termiska fissionsfaktorn 
    if n == 1
        C = C_0;
        deltaA_u235 = ((P*t)/(E_fission*1.602*10^-19))*(1-C)*(sigma_a235/sigma_f235);
        Nt_pu239 = deltaA_u235*C; %förbrukade fissila kärnor multiplicerat med konverteringskvoten
        Nt_u235 = N0tot_u235 - deltaA_u235; % antal kärnor U235 efter t tid (ett år)
        Nt_u238 = N0tot_u238 - Nt_pu239; % antal kärnor U238 efter t tid (ett år), ingen hänsyn tagen till att X antal kärnor U238 har förbrukats genom fission
    elseif n == 2
        C = C_t;
        deltaA_u235 = ((P*t)/(E_fission*1.602*10^-19))*(1-C)*((sigma_a235/sigma_f235)*(1-andel_PuvsU)+(sigma_a239/sigma_f239)*andel_PuvsU);
        Nt_pu239 = deltaA_u235*C; %förbrukade fissila kärnor multiplicerat med konverteringskvoten
        Nt_u235 = N0tot_u235 - deltaA_u235*(1-C); % antal kärnor U235 efter t tid (ett år)
        Nt_u238 = N0tot_u238 - Nt_pu239; % antal kärnor U238 efter t tid (ett år), ingen hänsyn tagen till att X antal kärnor U238 har förbrukats genom fission
    end 
    
    
    disp(Nt_pu239);
    disp(Nt_u238);
    disp(Nt_u235);
    
    if n == 1
        andel_PuvsU = 0; % andel Pu239 av de fissila kärnorna, noll i början
    end
   
    ny = (1-andel_PuvsU)*ny_u235 + (andel_PuvsU)*ny_pu239;

    eta = ny * (e*(sigma_f235*(1-andel_PuvsU)+sigma_f239*andel_PuvsU)/(e*(sigma_a235*(1-andel_PuvsU)+sigma_a239*andel_PuvsU)+(1-e)*sigma_c238));




    %%%%%%%%%% Beräkning av termiska utnyttjandefaktorn
%Sätter alla fi till 1 då vi antar homogent neutronflöde
    
    if n == 1
        Sigma_abransle = Sigma_a235+Sigma_a238; %Makroskopiska tvärsnittet för bränslet, ändras över tiden
        insattning = 0.711; %[%] insättning av styrstaven, oklart vad den ska börja på, vi kan börja med noll och se vilken k_eff vi får
    elseif n == 2
        Sigma_abransle = Sigma_a235+Sigma_a238+Sigma_a239; %Makroskopiska tvärsnittet för bränslet, ändras över tiden
        insattning = 0.383; %[%] insättning av styrstaven, oklart vad den ska börja på, vi kan börja med noll och se vilken k_eff vi får
    end
    fi_bransle = Sigma_abransle/FR; %neutronflöde kopplat till absorption för U235
    %disp(['fi_bransle = 'num2str(fi_bransle)]);
    V_bransle = 1; % typiskt volymförhållande

    fi_mod = 1; 
    Sigma_amod = 0.02; %taget från KSU (I think)
    V_mod = 3; % typiskt volymförhållande

    fi_styrstav = 1;
    Sigma_astyrstav = 0.15; %[cm-1]
    styrstav_bredd = 270*10^-3; %[m]
    styrstav_tjocklek = 7*10^-3; %[m] vet inte vad tjockleken är för forsmarks styrstavar, men värdet stämmer överens med Sigma_astyrstav
    styrstav_langd = 3700*10^-3; %[m], antar att det inte spelar någon roll när vi har en längd insättning av styrstaven som vi styr
    N_styrstav = 161; %antal styrstavar
    V_styrstav = styrstav_bredd*styrstav_tjocklek*styrstav_langd*insattning*N_styrstav; %rätt?

    fi_gift = 1;
    sigma_agift = 2.65*10^6; %[barn] mikroskopiskt tvärsnitt för Xe135, taget från KSU
    gamma_TeU = 0.0321618; % fissionsproduktfördelning för Tellurium från U235
    gamma_IU = 0.0292737;% fissionsproduktfördelning för jod från U235
    T_I2 = 6.57*60*60;%[s] halveringstid för jod
    lambda_I = log(2)/T_I2;%sönderfallskonstanten för I135
    
    gamma_TePu = 0.0219297; % fission yield Tellurium från Pu239
    gamma_IPu = 0.0428742; % fission yield Iodine från Pu239
    gamma_XeU = 0.00178122; % fission yield Xenon från U235
    gamma_XePu = 0.007522799; % fission yield Xenon från Pu239
    
    if n == 1
        gamma_I = gamma_TeU+gamma_IU; %total fission yield för Jod, när vi antar att det inte finns någon Pu239 i reaktorn
        gamma_Xe = gamma_XeU; %total fission yield för Xenon, när vi antar att det inte finns någon Pu239 i reaktorn
    elseif n == 2
        gamma_I = (gamma_TeU+gamma_IU)*(1-andel_PuvsU)+(gamma_TePu+gamma_IPu)*andel_PuvsU; % total fission yield för Jod, när vi vet att det finns Pu239 i reaktorn
        gamma_Xe = gamma_XeU*(1-andel_PuvsU) + gamma_XePu*andel_PuvsU; % total fission yield för Xenon, när vi vet att det finns Pu239 i reaktorn
    end
    
    N_I = (gamma_I*FR)/lambda_I; %antalet kärnor Jod i härden
    N_A = 6.022*10^23; %Avogadros konstant
    M_I = 134.91006; %[u] relativ atommassa
    m_I = (N_I/N_A)*M_I; %[g] massa jod i härden
    %disp(['m_I = ' num2str(m_I) ' g']);
    
    T_Xe2 = 9.14*60*60; % [s] halveringstid Xe135
    lambda_Xe = log(2)/T_Xe2; % [s] sönderfallskonstant för Xe135
    sigma_aXe = 2.65*10^6*10^-24; % [cm-2] mikroskopiskt absorptionstvärsnitt för Xe135
    %N_Xe = (lambda_I*N_I+gamma_Xe*FR)/(sigma_aXe+lambda_Xe); % antal kärnor Xe135 i härden
    N_Xe = ((gamma_Xe+gamma_I)*FR)/(gamma_Xe+sigma_aXe*fi_bransle); % antal kärnor Xe135 i härden - denna tagen från slide (är nog denna)
    M_Xe = 134.9072075; % [u]Atomic mass Xe135
    m_Xe = (N_Xe/N_A)*M_Xe; % massa Xe135 i härden'
    %disp(['m_Xe = ' num2str(m_Xe) ' g']);

    V_gift =  N_Xe/((N0tot_u235+N0tot_u238)*2*16);%N_Xe / (N_tot); N_tot är endast rätt vid start, men ett bra antagande annars. Ignorera O2 ?
    Sigma_agift = N_Xe*sigma_aXe; % Makroskopiska tvärsnittet för 

    fi_an = 1; 
    Sigma_aan = 1;
    V_an = 1;

    %f = (fi_bransle*Sigma_abransle*V_bransle)/(fi_bransle*Sigma_abransle*V_bransle+fi_mod*Sigma_amod*V_mod+fi_styrstav*Sigma_astyrstav*V_styrstav+fi_gift*Sigma_agift*V_gift+fi_an*Sigma_aan*V_an); %V_an * mm = 1.1 (?)
    
    disp(['V_styrstav = ' num2str(V_styrstav)])
    
    %%%%%%%%%% Beräkning av neutronflödet
    %fi = P/(sigma_f235*gamma*N_235); 

    %%%%%%%%%% p*f plot
    
    f_all = zeros(1,length(X));
    for i = 1:1:length(X)
        V_m = X(i);
        f_all(i) = (Sigma_abransle*V_bransle)/(Sigma_abransle*V_bransle+Sigma_amod*V_m+Sigma_astyrstav*V_styrstav+Sigma_agift*V_gift+0.15);
    end
    
    f = (Sigma_abransle*V_bransle)/(Sigma_abransle*V_bransle+Sigma_amod*V_mod+Sigma_astyrstav*V_styrstav+Sigma_agift*V_gift+0.15);
    k = eta * epsilon * p * f * P_leakage;
    
    disp(['P_leakage = ' num2str(P_leakage)]);
    disp(['epsilon = ' num2str(epsilon)]);
    disp(['f = ' num2str(f)]);
    disp(['p = ' num2str(p)]);
    disp(['eta = ' num2str(eta)]);
    disp(['k_eff = ' num2str(k)]);
    
    pf_all = f_all.*p_all;
    [Max,Index] = max(pf_all);
    
    figure(n)
    plot(X,f_all);
    axis([0 10 0 1.1])
    if n == 1
        title('Driftstart');
    elseif n == 2
        title('Slutet av drifttid');
    end
    ylabel('');
    xlabel('V_{mod}/V_{bransle}');
    hold on
    plot(X,p_all);
    plot(X,pf_all);
    plot(Index/100,Max,'o');
    if n == 1
        plot(3,0.539,'o');
    elseif n == 2
        plot(3,0.55,'o');
    end
    legend('Termisk Utnyttjandefaktor (f)','Resonanspassagefaktor (p)','f*p','Maximum p*f', 'Reaktorns V_{mod}/V_{bransle}')
    hold off
    
    
    disp([' ']);
    
    if n == 2
        andel_PuvsU = Nt_pu239 / Nt_u235; % andel Pu239 av de fissila kärnorna
    end
    n = n + 1;
end