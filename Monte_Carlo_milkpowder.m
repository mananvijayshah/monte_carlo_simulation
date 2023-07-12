%Soft particle CVMC model (July-2022)
%                                                            By: Manan Shah
% Droplets:  Monosized
% Particles: Monosized, soft and porous
% No Overspray
% No Breakage


clear all
close all
clc
tic

%To save the entire simulation
file=['CVMC_Real_SMP_final-' date];
save_simulation='Yes';
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA THAT MUST BE GIVEN

FColl=35;                   %Pre-factor of the collision frequency
drying=1;                   %drying=1 when drying of deposited droplets is activated
breakage=0;                  %breakage=1 when breakage is activated
matrix_size=50;             %Matrix for monitoring St, St*

%Process parameters
Tgas=50;                       % Inlet gas  Temperature [°C]
Tequilibrium=18;            %Adiabatic saturation temperature [°C]
Hnozzle=150;                %Height of the tip of the nozzle respect to the distributor [mm]
Y_in = 0.001;                 %Inlet moisture content [kg/kg]
M_gas=24;                     %Mass flow rate of the fluidization gas [kg/h]
Dt=0.152;                      %inner diameter of granulator
M_spray=0.6;                   %Spraying rate [kg/h]
Y_bulk=0.0054;                      %moisture content of bulk gas [kg/kg]
% M_drygas=M_gas/(1+Y_in);
% Y_bulk=Y_in+(M_spray/M_drygas);  %moisture content of bulk gas [kg/kg]

%Solid
primpartdiameter=2.0000e-04; %Primary particle diameter [m]
bed_mass=600;               %Bed mass (batch) [gr]
partdensity=1300;           %Solid density [Kg/m3]
e=0.5;                        %Restitution coefficient [-]
ha=5e-6;                    %Height of particle asperities [m]
particle_porosity=0.75;      %Intra-particle porosity, substrate porosity [-]

%Binder(water)
diam_drop=0.04;             %dd1 droplet diameter [mm]
teta=40;                    %Contact angle [Degrees]

 
%Program parameters
Np=1000;                    %Initial number of primary particles and entities
collnumber=2000;              %ITERATIONS
clases=Np;                  %Clases for granulometric analysis (One particle per class)
h1max=Np;                   %Upper limit for granulometric analysis (particle based)
sec_PSD=10;                 %Time counting for a posterior PSD analysis,
%"sampling" each sec_PSD [s, process time]
save_file=10;               %Time conunting for saving the entire simulation [s, process time]
uc_sigma=0.1;               %collision velocity standard deviation [m/s]

%Glass transition temp, viscosity and solid temp constants 
C=-17.4;
B=51.6;                      %[K]
liqden = 1000;               %Water density [Kg/m3]
Tgs = 101 + 273.15;          %Glass transition temperature of solid [K]
Tgw = -137 + 273.15;         %Glass transition temperature of water [K]
const = 6.5;                   %Gordon-Taylor constant [-]
mug = 1e+12;               %viscosity at the glass transition temperature [Pa.s]
Cair = 29.15;                %Specifice [K heat of air [KJ/Kmol/K]
Cvap =  33.606;              %Specific heat of vapours [KJ/Kmol/K]
Hvap = 43774;                %Specific heat of vapourization [KJ/Kmol]

%Fractal properties
% Df =(0.0163*Tgas)-(7.097*M_spray)+ 1.698;                   %fractal dimension of the agglomerate
% kf = 5.85-(1.42*Df);                                                   %Prefactor of the agglomerate
Df =2.35;
kf =1.25;

% droplet addition
dropletrate=((M_spray*(1000/3600))/bed_mass)*(partdensity/liqden)*(((primpartdiameter*1000)/diam_drop)^3);         %Droplet per particle per second (From file PARAMETERS) [-]



dropdiameter=diam_drop ;     %droplet diameter [mm]
process_time=zeros(collnumber,1);
x_interval=zeros(collnumber,1);
diameter=zeros(collnumber,1);
primpartmass=((4*pi*((primpartdiameter/2)^3))/3)*partdensity;
num_part_real=(bed_mass/1000)/primpartmass;
partvolumen_real=num_part_real*primpartmass/partdensity; %Total volume
airdensity=0.00000000904*(Tgas)^3+0.00000944905*(Tgas)^2-0.00425096197*Tgas+1.2878411577; %[Kg/m3]
airviscosity=0.0000000389*Tgas+0.00001779; %[Pa.s]
granulom_numb=0;



Airvolumetricflow=M_gas/(3600*airdensity); %volumetric Air flow through the distributor plate 0.0238 [m3/s]
Area=(pi/4)*Dt^2;
u0=0.35;
% u0=(Airvolumetricflow/Area); %Fluidization velocity (m/s)
[ho,A,vol_solid,dropletvolume,diamdropletonparticle] = function_height( partdensity,Tequilibrium,particle_porosity,Tgas,dropdiameter,teta,u0,airdensity,airviscosity,Y_bulk);
Acap=pi*diamdropletonparticle*diamdropletonparticle*1e-6;
Ts = ((Tgas + 273.15)*(Cair + Y_in * Cvap) + (Y_in - Y_bulk)*Hvap)/( Cair + Y_bulk * Cvap); %% solid temperature[K]
M_solid = vol_solid*partdensity;             %mass of solid material[Kg]
M_wo= dropletvolume*liqden;                  %initial mass os water in droplet [Kg]
positions=ceil((pi.*(primpartdiameter).^2)./Acap);  %

%Single Matrices
a=1:Np;a=a';
b=ones(Np,1); %Initial number of particles per entity = 1
reportA=cat(2,a,b,b*positions,b*primpartdiameter,b*primpartmass);
reportB=zeros(Np,positions);
reportBage=zeros(Np,positions);
reportBmoisture=zeros(Np,positions);        % moisture content of the puddle
reportBmassfraction=zeros(Np,positions);          %glass transition proprtey of puddle
reportBinitial_moisture=zeros(Np,positions);
reportBfactA=zeros(Np,positions);
frac=0.25;
% Parameters for monitoring, e.g., agglomerate and droplets properties,
% number and nature of the events,
% internal counters for time dependent mechanisms,
% internal counters for saving simulation results.

time=0;
timerate=0;
time_counter_PSD=0;
coalescence_num=zeros(collnumber,2);
breakage_num=zeros(collnumber,3);
doublings=0;%Number of doublings of the simulation box (CVMC method)
breakage_per_event=0;%Number of wet breakages
breakage_per_event_dry=0;%Number of dry breakages
time_conter_file=0;
time_screen=0;

counter=1;
counter_Stokes_after=1;
counter_Stokes_before=1;


Stokes_before=zeros(collnumber,matrix_size);
Stokes_crit_before=zeros(collnumber,matrix_size);
Stokes_after=zeros(collnumber,matrix_size);
Stokes_crit_after=zeros(collnumber,matrix_size);


spray_droplets=0;
sprayed_droplets_matrix=zeros(collnumber,1);
sprayed_droplets_matrix_kind=zeros(collnumber,1);

droplets_arrive=0;
droplets_arrive_matrix=zeros(collnumber,1);
droplets_arrive_matrix_kind=zeros(collnumber,1);

droplets_lost=0;
droplets_lost_matrix=zeros(collnumber,1);
droplets_lost_matrix_kind=zeros(collnumber,1);

droplets_elut=0;
droplets_elut_matrix=zeros(collnumber,1);
droplets_elut_matrix_kind=zeros(collnumber,1);

Hexp_matrix=zeros(collnumber,1);
aditiontime=(1*1)/(dropletrate*Np); %time to add a droplet

Doubling_Factor=zeros(collnumber,1);
agglom_diameters=zeros(collnumber,2);

dropl_event=0;
counter_spray_kind=1;
counter_arrive_kind=1;
counter_lost_kind=1;
counter_elut_kind=1;

%=====================================
%% %Begining of the simulation
for k=1:collnumber %Number of pair-wise particle collisions
    
    if k==1
        dianm_old=primpartdiameter;
    end
%     =============== Stepwise constant volume Monte Carlo Method ================
    x=size(reportA,1);
    if x<=(Np/2)
        reportA=cat(1,reportA,reportA);
        Ix                    = size(reportA, 1); % Ix = number of rows
        reportA(:,1)          = (1:Ix)';          % arrange first column 1 2 3..
        reportB=cat(1,reportB,reportB);
        reportBage=cat(1,reportBage,reportBage);
        reportBmoisture=cat(1,reportBmoisture,reportBmoisture); %#ok<NASGU>
        reportBinitial_moisture= cat(1,reportBinitial_moisture,reportBinitial_moisture);
        reportBfactA= cat(1,reportBfactA,reportBfactA);
        reportBmassfraction=cat(1,reportBmassfraction,reportBmassfraction); %#ok<NASGU>
        doublings=doublings+1;
        aditiontime=aditiontime/2;
    end
    
    Doubling_Factor(k,1)=1/(2^doublings);
    
    %Parameters to calculate collision frequency and FB voidage
    [colifreq,Hexp] = function_frequency(dianm_old,u0,airdensity,airviscosity,partdensity,Area,partvolumen_real,FColl);
    interval=1/colifreq; %time step
    x_interval(k)=interval;
    time=time+interval; %Real time
    process_time(k)=time;%Saving the process time
    timerate=timerate+interval; %Time counting for continous droplet addition based on real process time
    x1=size(reportA,1);%Number of entities
    wet=floor(frac*x1);
    
    %==============
    %Drying of deposited droplets
    %==============
    if drying==1 %Secado = 1 significa q tiene secado
        t_dry=interval;
    else
        t_dry=0;
    end
    
    places2=find(reportB==1);
    reportBage(places2)=reportBage(places2)+(t_dry);
    dry_positions=find(reportB==0);
    dry_positions83=find(reportB==83);
    reportBmoisture = reportBinitial_moisture - (reportBfactA.*reportBage); %[Kg]
    reportBmoisture(dry_positions)=0;
    reportBmoisture(dry_positions83)=83;
    reportBmassfraction = reportBmoisture./(reportBmoisture+M_solid);
    reportBmassfraction(dry_positions)=0;
    reportBmassfraction(dry_positions83)=83;
   
    %    Vanish droplets that dry out completely
    dry_drops=find( reportBmoisture<0);
    % finding average drying time of the droplets
    AverageDryTime(k)   =   sum (reportBage (dry_drops)) / length(dry_drops);
    %      making positions of dryed droplest zero in other matrices
    reportB(dry_drops)=0;
    reportBage(dry_drops)=0;
    reportBmoisture(dry_drops)=0;
    reportBmassfraction(dry_drops)=0;
    reportBinitial_moisture(dry_drops)=0;
    reportBfactA(dry_drops)=0;
    
    droplet_after_drying(k)=nnz(reportB==1);
        
    if timerate>=aditiontime
        dropl_event=floor(timerate/aditiontime);
        for u=1:dropl_event
            Hexp_matrix(k,1)=Hexp;
            spray_droplets=spray_droplets+1;
            sprayed_droplets_matrix(k)=dropl_event;
            sprayed_droplets_matrix_kind(k,counter_spray_kind)=counter;
            counter_spray_kind=counter_spray_kind+1;
            
            wet_particles=randperm(x1);
            wet_particles1=wet_particles(1:wet);
            wet_positions_matrix=reportB(wet_particles1,:);
            empty_positions=find(wet_positions_matrix==0);
            pickup_positions=randperm(max(size(empty_positions)));
            pickup_positions1=pickup_positions(1:1);
            wet_positions_matrix(empty_positions(pickup_positions1))=1;
            reportB(wet_particles1,:)=wet_positions_matrix(1:end,:);
            [row, col] = ind2sub(size(wet_positions_matrix),empty_positions(pickup_positions1));
            particle=wet_particles1(row);
            reportBfactA(particle,col)=A(counter);
            
            %Initial moisture content of the binder
            reportBinitial_moisture(particle,col)=M_wo(counter); %[m]
            
            wet_positions_matrix=reportBage(wet_particles1,:);
            wet_positions_matrix(empty_positions(pickup_positions1))=0;
            reportBage(wet_particles1,:)=wet_positions_matrix(1:end,:);
            

            reportBmoisture = reportBinitial_moisture - (reportBfactA.*reportBage); %[Kg]
            dry_positions=find(reportB==0);
            dry_positions83=find(reportB==83);
            reportBmoisture(dry_positions)=0;
            reportBmoisture(dry_positions83)=83;
            reportBmassfraction = reportBmoisture./(reportBmoisture+M_solid);
            reportBmassfraction(dry_positions)=0;
            reportBmassfraction(dry_positions83)=83;
            reportBinitial_moisture(dry_positions)=0;
            reportBinitial_moisture(dry_positions83)=83; 
            reportBfactA(dry_positions)=0;
            reportBfactA(dry_positions83)=83;
            
            droplets_arrive=droplets_arrive+1;
            droplets_arrive_matrix(k)=droplets_arrive;
            droplets_arrive_matrix_kind(k,counter_arrive_kind)=counter;
            counter_arrive_kind=counter_arrive_kind+1;
            droplets_lost=0;
            droplets_lost_matrix(k)=droplets_lost;
            droplets_lost_matrix_kind(k,counter_lost_kind)=0;
            counter_lost_kind=1;
            droplets_elut=0;
            droplets_elut_matrix(k)=droplets_elut;
            droplets_elut_matrix_kind(k,counter_elut_kind)=0;
            counter_elut_kind=1;
            counter=1;
            
        end
        
        timerate=timerate-(dropl_event*aditiontime);
        counter_spray_kind=1;
        counter_arrive_kind=1;
        counter_lost_kind=1;
        counter_elut_kind=1;
        droplets_arrive=0;
        droplets_lost=0;
        droplets_elut=0;
        
    end
    
    droplet_added(k)=nnz(reportB==1);
    [numelereport,numpositions]=size(reportB);
    
    %Two aleatory groups for collison g1 and g2
    
    random_groups = randperm(numelereport)';
    
    pair = floor(numelereport/2);
    g1   = random_groups(1:pair);
    g2   = random_groups(pair+1:2*pair);
    
    numagglo=0;
    breakage_per_event=0;
    breakage_per_event_dry=0;
    
    for i=1:pair %Begining of the collisions i of the event k
        positions_part1=reportA(g1(i),3);%Positions of agglomerate 1
        positions_part2=reportA(g2(i),3);%Positions of agglomerate 2
        permut_drops1=ceil(positions_part1*rand(1));%Choose randonly one colliding position
        permut_drops2=ceil(positions_part2*rand(1));%Choose randonly one colliding position
        collision_kind=(reportB(g1(i),permut_drops1)+reportB(g2(i),permut_drops2));
        d1=reportA(g1(i),4);
        d2=reportA(g2(i),4);
        grandiametro=(2*d1*d2)/(d1+d2);%Average diameter of the two colliding particles (sphere based) [m]
        m1=reportA(g1(i),5); %[Kg]
        m2=reportA(g2(i),5); %[Kg]
        M_agg=(2*m1*m2)/(m1+m2);%Average mass of the two colliding particles %[Kg]
        Y=rand;
        uc_medio=u0/2; %[m/s]
        uc=normrnd(uc_medio,uc_sigma);
        height_collision=ho; %[m] 
         %first criteria for wet collision
        if collision_kind>0  %collision_kind>0 means wet collision
            w1=reportBmassfraction(g1(i),permut_drops1);
            w2=reportBmassfraction(g2(i),permut_drops2);
            
            %Glass transition Temperature
            Tg1=((1 -w1)*Tgs + const*Tgw*w1)/((1- w1) + const*w1); % [K]
            Tg2=((1 -w2)*Tgs + const*Tgw*w2)/((1- w2) + const*w2); % [K]
            
            % second criteria 
            if Ts>Tg1+20 || Ts>Tg2+20
           
                gt=[Tg1,Tg2];
                Tg=max(gt(Ts>gt+20));
                
                A1 = C*(Ts-Tg);                 %%[K]
                B1 = (B + (Ts-Tg));             %%[K]
                mu=(mug)*(10^(A1/B1));       %% viscosity at wet position [Pa.s]
                
                tgcheck(k,counter_Stokes_before)=Tg;
                tg1check(k,counter_Stokes_before)=Tg1;
                tg2check(k,counter_Stokes_before)=Tg2;
                mucheck(k,counter_Stokes_before)=mu;
                St=(2*M_agg*uc)/(3*pi*mu*((grandiametro/2)^2));
                Stcritical=(1+(1/e))*log(height_collision/ha);
                Stokes_before(k,counter_Stokes_before)=St;
                Stokes_crit_before(k,counter_Stokes_before)=Stcritical;
                counter_Stokes_before=counter_Stokes_before+1;
            
            if St <= Stcritical  %Stokes criterion is fulfilled
                Stokes_after(k,counter_Stokes_after)=St;%Saving the St for a posterior analysis
                Stokes_crit_after(k,counter_Stokes_after)=Stcritical;%Saving the critical St for a posterior analysis
                counter_Stokes_after=counter_Stokes_after+1;
                
                reportA(g2(i),1)=0;%Making zero en ID of the entity that disappers due to coalescence
                reportA(g1(i),2)=reportA(g1(i),2)+reportA(g2(i),2);%Adding the number of particles to the 1st entity
                positions_part_g1=reportA(g1(i),3);%Recalculate the number of positions
                positions_part_g2=reportA(g2(i),3);

                % equivalent diameter of the final agglomerate
                reportA(g1(i),4)=1.291*primpartdiameter*((reportA(g1(i),2)/kf).^(1/Df));
                reportA(g1(i),3)=ceil((pi.*(reportA(g1(i),4)).^2)./Acap);         %to recalculate the agglomerate new positions
                reportA(g1(i),5)=primpartmass*reportA(g1(i),2);
                
                %ReportB Changes
                %%%1. to change the wet position that form a bridge
                
                reportB(g1(i),permut_drops1)                   = 0;
                reportBage(g1(i),permut_drops1)                = 0;
                reportBmoisture(g1(i),permut_drops1)           = 0;
                reportBinitial_moisture(g1(i),permut_drops1)   = 0;
                reportBfactA(g1(i),permut_drops1)              = 0;
                reportBmassfraction(g1(i),permut_drops1)       = 0;
                
                reportB(g2(i),permut_drops2)                   = 0;
                reportBage(g2(i),permut_drops2)                = 0;
                reportBmoisture(g2(i),permut_drops2)           = 0;
                reportBinitial_moisture(g2(i),permut_drops2)   = 0;
                reportBfactA(g2(i),permut_drops2)              = 0;
                reportBmassfraction(g2(i),permut_drops2)       = 0;
                
                %%%2. To set matrix size
                
                temp2              =   size(reportB, 2);
                temp3              =   max(reportA(:,3));
                
                reportB(:,temp2+1:temp3)                   = 83; %add positions to the array
                reportBage(:,temp2+1:temp3)                = 83; %positions with 83 mean no existing positions
                reportBmoisture(:,temp2+1:temp3)           = 83; %However, it is necesary to fill out the matrix array
                reportBinitial_moisture(:,temp2+1:temp3)   = 83;
                reportBfactA(:,temp2+1:temp3)              = 83;
                reportBmassfraction(:,temp2+1:temp3)       = 83;
                
                %%%3. Moveing g2 position in g1 and rearranging
                
                combine_positions_agglomerate=positions_part_g1+positions_part_g2;
                calculated_position_agglomerate=reportA(g1(i),3);
                
                
                %Rearrange all; reportB, reportBage and reportBheight etc... for g1
                if combine_positions_agglomerate<= calculated_position_agglomerate
                    
                    reportB(g1(i),1:calculated_position_agglomerate)= cat(2,reportB(g1(i),1: positions_part_g1),reportB(g2(i),1:positions_part_g2),zeros(1,calculated_position_agglomerate-combine_positions_agglomerate));
                    reportBmoisture(g1(i),1:calculated_position_agglomerate)= cat(2,reportBmoisture(g1(i),1: positions_part_g1),reportBmoisture(g2(i),1:positions_part_g2),zeros(1,calculated_position_agglomerate-combine_positions_agglomerate));
                    reportBage(g1(i),1:calculated_position_agglomerate)= cat(2,reportBage(g1(i),1: positions_part_g1),reportBage(g2(i),1:positions_part_g2),zeros(1,calculated_position_agglomerate-combine_positions_agglomerate));
                    reportBinitial_moisture(g1(i),1:calculated_position_agglomerate)= cat(2,reportBinitial_moisture(g1(i),1: positions_part_g1),reportBinitial_moisture(g2(i),1:positions_part_g2),zeros(1,calculated_position_agglomerate-combine_positions_agglomerate));
                    reportBfactA(g1(i),1:calculated_position_agglomerate)= cat(2,reportBfactA(g1(i),1: positions_part_g1),reportBfactA(g2(i),1:positions_part_g2),zeros(1,calculated_position_agglomerate-combine_positions_agglomerate));
                    reportBmassfraction(g1(i),1:calculated_position_agglomerate)= cat(2,reportBmassfraction(g1(i),1: positions_part_g1),reportBmassfraction(g2(i),1:positions_part_g2),zeros(1,calculated_position_agglomerate-combine_positions_agglomerate));
                    
                else
                    R_B=cat(2,reportB(g1(i),1: positions_part_g1),reportB(g2(i),1:positions_part_g2));
                    R_Bmoisture=cat(2,reportBmoisture(g1(i),1: positions_part_g1),reportBmoisture(g2(i),1:positions_part_g2));
                    R_Bage=cat(2,reportBage(g1(i),1: positions_part_g1),reportBage(g2(i),1:positions_part_g2));
                    R_Bheight_initial=cat(2,reportBinitial_moisture(g1(i),1: positions_part_g1),reportBinitial_moisture(g2(i),1:positions_part_g2));
                    R_BfactA=cat(2,reportBfactA(g1(i),1: positions_part_g1),reportBfactA(g2(i),1:positions_part_g2));
                    R_Bmassfraction=cat(2,reportBmassfraction(g1(i),1: positions_part_g1),reportBmassfraction(g2(i),1:positions_part_g2));
                    
                    
                    delete_drydroplet=find(R_B==0,combine_positions_agglomerate-calculated_position_agglomerate);
                    
                    R_B(delete_drydroplet)=[];
                    R_Bmoisture(delete_drydroplet)=[];
                    R_Bage(delete_drydroplet)=[];
                    R_Bheight_initial(delete_drydroplet)=[];
                    R_BfactA(delete_drydroplet)=[];
                    R_Bmassfraction(delete_drydroplet)=[];
                    
                    
                    reportB(g1(i),1:calculated_position_agglomerate)= R_B;
                    reportBmoisture(g1(i),1:calculated_position_agglomerate)= R_Bmoisture;
                    reportBage(g1(i),1:calculated_position_agglomerate)= R_Bage;
                    reportBinitial_moisture(g1(i),1:calculated_position_agglomerate)=R_Bheight_initial;
                    reportBfactA(g1(i),1:calculated_position_agglomerate)=R_BfactA ;
                    reportBmassfraction(g1(i),1:calculated_position_agglomerate)=R_Bmassfraction;
                    
                end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
                numagglo=numagglo+1;
            end
            end
        end  %del for de los pares de colisiones
    end
    
    counter_Stokes_after=1;
    counter_Stokes_before=1;
    counter_mu=1;
    
    ind                           = reportA(:,1)==0;  % find zeros in first column
    reportA(ind,:)                = [];               % delete those zeros rows
    I1                            = size(reportA, 1); % I1 = number of rows
    reportA(:,1)                  = (1:I1)';          % arrange first column 1 2 3 4..
    reportB(ind,:)                = [];               % delete the same rows from reportB
    reportBage(ind,:)             = [];               % delete the same rows from reportBage
    reportBmoisture(ind,:)        = [];               % delete the same rows from reportBheight
    reportBinitial_moisture(ind,:)  = [];               % delete the same rows from reportBheight_initial
    reportBfactA(ind,:)           = [];               % delete the same rows from reportBfactA
    reportBmassfraction(ind,:)          = [];               % delete the same rows from reportBheight_porosity
    
    totalparticles=sum(reportA(:,2));
    droplet_lost_collison(k)=nnz(reportB==1);
    %=====================================
    %Granulometric analysis
    
    diameter(k)=function_diameter(reportA(:,4));
    dianm_old=diameter(k);
    time_counter_PSD=time_counter_PSD+interval;

    %=====================================
%     Sampling the population for a posterior anaylsis of the evolution of PSD
    if time_counter_PSD>=sec_PSD & save_simulation=='Yes'
        time_counter_PSD=0;
        granulom_numb=granulom_numb+1;
        save([file '-reportA-' num2str(granulom_numb) '.mat'], 'reportA');
    end
    
    time_conter_file=time_conter_file+interval;
    
    %Saving the entire simulation in a file for a posterior analysis
    if time_conter_file>=save_file
        time_screen=time_screen+save_file;
        time_conter_file=0;
        
        if save_simulation=='Yes'
            save(['SIMULATION-' file '.mat']);
        end
        msg=['Progress: ' num2str(time_screen)  ' seconds' ]; disp(msg)
        relative_diameter=(((diameter(k))-(diameter(1)-primpartdiameter))/primpartdiameter)
        
    end
         
    coalescence_num(k,1)=k;
    coalescence_num(k,2)=numagglo;
    breakage_num(k,1)=k;
    breakage_num(k,2)=breakage_per_event;
    breakage_num(k,3)=breakage_per_event_dry;
    

%%     %% Terminating the loop after 600 seconds==============================
    if time_screen >= 600
        break
    end
end
growthrate = 1e6.*(mean(reportA(:,4))-primpartdiameter)./time_screen; disp(growthrate)
toc

