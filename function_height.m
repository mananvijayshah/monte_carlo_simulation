function [ho,A,vol_solid,dropletvolume,diamdropletonparticle] = function_height( partdensity,Tequilibrium,particle_porosity,Tgas,dropdiameter,teta,u0,airdensity,airviscosity,Y_bulk)

P=100400;                   %System Pressure [Pa]
liqdensity=1000;         %water density [Kg/m3]1021.07
%Tequilibrium=36;           %ES NECESARIO CALCULAR ESTA TEMPERATURA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Pvapor=(10^(8.0713-((1730.63/(Tequilibrium+233.426)))))*(101325/760);%Pascales
diffcoefficient=(2.252/P)*((Tgas+273.15)/273.15)^1.81; %m2/s

dropletvolume=(4*pi*((((dropdiameter/2)/1000)).^3))/3; %m3
diamdropletonparticle=((((3*(dropletvolume)/pi))*((sind(teta)).^3)/(2-(3*cosd(teta))+((cosd(teta)^3)))).^(1/3))*1000; %mm
area_drop = pi*(diamdropletonparticle/1000).^2;   %m2
ho = dropletvolume/(area_drop*particle_porosity);  %m
caraclenght=dropdiameter/1000;

% mass transfer coefficient 
Redrop=((caraclenght)*u0*airdensity)/airviscosity;
Sc=airviscosity/(airdensity*diffcoefficient);
Sh1=0.664*(Redrop.^0.5)*(Sc^(1/3));
Sh2=(0.037*(Redrop.^0.8)*Sc)/(1+2.443*(Redrop.^-0.1)*(Sc^(2/3)-1));
Sh=(8/pi)+((Sh1.^2)+(Sh2.^2)).^0.5;
beta=Sh*diffcoefficient./(caraclenght); %Coeficiente de transferencia de masa [m/s]
Ysat=0.620*(Pvapor/(P-Pvapor));
A = airdensity*beta*area_drop*(Ysat-Y_bulk); %Kg/sec
vol_solid = (1-particle_porosity)*dropletvolume/particle_porosity;
end