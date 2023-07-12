function [colifreq,Hexp] = function_frequency(d32old,u0,airdensity,airviscosity,partdensity,Area,partvolumen_real,FColl)

%   Korina Terrazas Velarde 19.07.2007
    Re=(airdensity*u0*d32old)/airviscosity;%(-)
    Ar=((d32old^3)*airdensity*(partdensity-airdensity)*9.81)/(airviscosity^2);%(-)
    Eexpanded=(((18*Re)+0.36*(Re^2))/Ar)^0.21;
    %umf=(0.0008*9.81*(partdensity-airdensity)*(d32^2))/airviscosity;%(m/s)
    if Eexpanded<=0.39 
        msg='Bed can not be fluidized';
        disp(msg)
        return
    end
    if Eexpanded>=1 
        msg='Bed is being elutrated';
        disp(msg)
        return
    end

    %num_part_real=(masaparticulas/1000)/primpartmass;
    Hfixed=(partvolumen_real/(1-0.39))*(1/Area);%m
    Hexp=(((Hfixed*(1-0.39))/(1-Eexpanded)))*1000;%mm
    
    gfixed=0.61;
    gexp=partvolumen_real/((Hexp/1000)*Area);
    colifreq=FColl*((gexp/gfixed)^2)*(1-(gexp/gfixed))*u0;
    
    
