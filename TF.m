function [TF,Phase] = TF(M)
%This function recalculates thrust force at each step in order to plot
%Input: State Matrix
%Output: Thrust Vector

global Cdis At gam Vem R Pat Vair0 P0 mair0

for i = 1:size(M)
    Vair = M(i,8); %current volume of air in rocket
    mAir = M(i,7); %air mass
    if ~(Vair<Vem)
        Pend = P0*(Vair0/Vem)^gam; %Pressure 
        P2 = Pend*(mAir/mair0)^gam; %Internal pressure of rocket (phase 2)
        Vair = Vem;
    end
    if(Vair<Vem)
        %Phase 1:Water Propulsion
        Phase(i) = 1;
        P = ((Vair0/Vair)^gam)*P0; %air pressure function
        %calculating forces
        TF(i) = 2*Cdis*At*(P-Pat); %Thrust Force vector of rocket (N)
    elseif(Vair==Vem && P2>Pat)
        %Phase 2: Air Propulsion
        Phase(i) = 2;
        p = mAir/Vem;
        Pcrit = P2*((2/(gam+1))^(gam/(gam-1))); %Critical pressure
        T = P2/(R*p); %Air temperature

        %determining if flow choked
        if(Pcrit>Pat)
            %Exit Mach == 1
            Te = T*(2/(gam+1));
            pe = Pcrit/(R*Te);
            Ve = sqrt(gam*R*Te);
            Pe = Pcrit;
        else
            %Exit Mach != 1
            Me = sqrt((((P2/Pat)^((gam-1)/gam))-1)/((gam-1)/2));
            Te = T/(1+(((gam-1)/2)*(Me^2)));
            Ve = Me*sqrt(gam*R*Te);
            pe = Pat/(R*Te);
            Pe = Pat;
        end

        %variable rates
        TF(i) = (Cdis*pe*At*Ve*Ve)+((Pat-Pe)*At); %Thrust Force vector of rocket (N)
    else
        %Phase 3: Ballistic Flight
        Phase(i) = 3;
        TF(i) = 0;
    end
    
    
end


end

