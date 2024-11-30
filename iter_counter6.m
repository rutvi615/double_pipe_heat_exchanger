clc;
load('vfs.mat')
fflow=0.07;
fflow=0.5;
mdotmax = 1.5;
mdotmin = 0.07;
T_exh = 1004;   % Temperature of exhaust


mdot = 0.0693+0.025*(fflow-mdotmin)/(mdotmax-mdotmin);   %% mdot bound constraint so that the hot air outlet T is < 533.15 K for controling NOx emissions


T_CA_in = 300;  % Temprature of the inlet cold air
T_amb = 300;     % temperature of ambient air

%Radiation Constants
sig = 5.672e-8;   % W/(m^2-K^4) 
absorp = 0.9;
Lmb = 0.95;
eps = 0.9;


Cp_FG = 1.37e3; % J/kg-K
mdot_FG = mdot;
h_FG = 150;    % W/(m^2-K)    %convection of flue gas
mdot_air = mdot_FG*0.18;
Cp_air = 1e3;   %J/Kg-K
h_air = 150;  % W/(m^2-K)    % heat transfer coefficient of air inside pipe
h_amb = 100;  % W/(m^2-K)    % heat transfer coefficient from insulation to ambient air
k_ins   = 0.15;





r1 = 0.07; % inner radius of the recuperator
r2 = 0.04; % radius of suface 2

L_ins = 0.05; % thickness of insulation
r_ins = r1+L_ins;
hCyl = 0.6; %height of the cylinder
h = 0.1; %height of each discretized surface
A1 = 2*pi*r1*h; %area of each element in surface 1
A2 = 2*pi*r2*h; %area of each element in surface 2


R_amb = log(r_ins/r1)/(2*pi*k_ins*h) + 1/(h_amb*A1) ;  % resistance to heat transfer from surface 1 to ambient air



% initialize wall and gas temperatures
% T_w = 400*ones(3,6); % wall temperatures (surface, element)
T_w=[900 875 850 825 800 775; 600 580 560 540 520 500];
T_FG = 300*ones(6,1); % temperatures of flue gas 
T_CA = 300*ones(6,1); % temperatures of cold air

JLHS = zeros(12,12);  % Radiosity matrix
JRHS = zeros(12,1);  % Right-hand side vector
JRec = zeros(12,1);  % Radiosity vector

JLHS12 = zeros(12,12);  % Radiosity matrix
JRHS12 = zeros(12,1);  % Right-hand side vector
JRec12 = zeros(12,1);  % Radiosity vector
JLHS22 = zeros(12,12); 
JRHS22 = zeros(12,1);
JRec22 = zeros(6,1);

T_FG(1) = T_exh; 
T_CA (6) = T_CA_in;

T_wall = T_w';
T_wall = T_wall(:);

T_wall_store = T_wall;
T_FG_store = T_FG;
T_CA_store = T_CA;


% calculation of JLHS12 that depends upon area and viewfactors
 for i = 1:6
        for j = 1:12
            if i == j
                JLHS12(i,j) = A1 * sum(viewF(i,1:12)) + A1 * (eps) / (1 - eps);
            else
                JLHS(i,j) = -A1 * viewF(i,j);
            end
        end
    end
    
    for i = 7:12
        for j = 1:12
            if i == j
                JLHS12(i,j) = A2 * sum(viewF(i,1:12)) + A2 * (eps) / (1 - eps);
            else
                
               if (j<7)
                JLHS12(i,j) = -A2*viewF(i,j);       
               else
                JLHS12(i,j) = 0;
        end
            end
        end
    end
for i=7:12
    for j=7:12
    if (i==j)
        JLHS22(i,j) = A2*sum(viewF(i,7:12))+A2*(eps)/(1-eps);
    else
        JLHS22(i,j) = -A2*viewF(i,j);     
    end       
    end
end




% we solve the energy balance problem by the principle of superposition.
% First we look at the contribution due to convection and conduction heat
% transfer, then we address the radiation heat transfer. 

for k=1:50      %50


% flue gas section

for i=1:5
    T_FG(i+1) = T_FG(i) - (1/(mdot_FG*Cp_FG))*(h_FG*A1*(T_FG(i)-T_w(1,i)) + h_FG*A2*(T_FG(i)-T_w(2,i)));
end

% air section CA then HAT

for i=5:-1:1
    T_CA(i) = T_CA(i+1) + (h_air*A2/(mdot_air*Cp_air))*(T_w(2,i+1)-T_CA(i+1));
end



%  [T_wall(1:6) T_FG T_wall(7:12) T_HA T_wall(13:18) T_CA]
% T_FG
% 
% T_HA
% 
% T_CA

% solve for radiosities
% for radiosities calculation, the surfaces are indexed differently 
% 1   7   13
% 2   8   14
% .   .   .
% 6   12  18

% changing the index of wall temperatures

T_wall = T_w';
T_wall = T_wall(:);


for i=1:6
   JRHS12(i) = sig*T_wall(i)^4/((1-eps)/(eps*A1));
end

for i=7:12
    JRHS12(i) = sig*T_wall(i)^4/((1-eps)/(eps*A2));
end


JRec12 = JLHS12\JRHS12;

 JRec(1:12) = JRec12;


% Updating wall temperatures for surfaces 1 2 by using radiation heat flux

for i=1:6
    radFlux=0;
    for j=1:12
    if (i==j)
        radFlux = radFlux+A1*sum(viewF(i,1:12))*JRec(j);
    else
        radFlux = radFlux-A1*viewF(i,j)*JRec(j);
    end
    end
T_wall(i) = -(radFlux-(R_amb/R_amb)-(h_FG*A1*T_FG(i)))/(h_FG*A1+1/R_amb);
end



for i=7:12
    radFlux=0;
    for j=1:12
    if (i==j)
        radFlux = radFlux+A2*sum(viewF(i,1:6))*JRec(j);
    else
        if (i<7)
        radFlux = radFlux-A2*viewF(i,j)*JRec(j);
        else
        radFlux = radFlux+0;
        end
    end       
    end
   radFlux;
T_wall(i) = T_FG(i-7+1)-radFlux/(h_FG*A2);
end





for i=7:12
    JRHS22(i) = sig*T_wall(i)^4/((1-eps)/(eps*A2));
end

JRec22 = JLHS22(7:12,7:12)\JRHS22(7:12);

JRec(7:12) = JRec22
JLHS22;
JRHS22;

for i=13:18
    radFlux=0;
    for j=13:18
    if (i==j)
        radFlux = radFlux+A2*sum(viewF(i,13:18))*JRec(j-6);
    else
        radFlux = radFlux-A2*viewF(i,j)*JRec(j-6);
    end       
    end
T_wall(i-6) = T_CA(i-13+1)+radFlux/(h_air*A2);
end


T_wall_store = [T_wall_store  T_wall];

T_w = reshape(T_wall,[6,2])';

T_wall_store = [T_wall_store T_wall];
T_FG_store = [T_FG_store T_FG];
T_CA_store = [T_CA_store T_CA];




end

% 
  [T_wall(1:6) T_FG T_wall(7:12) T_CA]
 


subplot(3,1,1)
plot(T_FG_store')
subplot(3,1,2)
plot(T_CA_store')



% plot of gas tempratures

xax = 1:6;

plot(xax,T_FG,'color','r','LineWidth',2,'Marker','o','MarkerSize',10)
hold all
plot(xax,T_CA,'color','k','LineWidth',2,'Marker','o','MarkerSize',10)
hold off
xlabel('Element Index','FontSize',20)
ylabel('Temperature (K)','FontSize',20)
legend('Flue Gas', 'Cold Air')

set(gca,'Fontsize',20)


%% plot of furnace surface tempratures

figure

xax = 1:6;

plot(xax,T_wall(1:6),'color','r','LineWidth',2,'Marker','o','MarkerSize',10)
hold all
plot(xax,T_wall(7:12),'color','b','LineWidth',2,'Marker','o','MarkerSize',10)

hold off
xlabel('Element Index','FontSize',20)
ylabel('Temperature (K)','FontSize',20)
legend('Surface 1', 'Surface 2')

set(gca,'Fontsize',20)

%plot of gas tempratures

xax = 1:6;

plot(xax,T_FG,'color','r','LineWidth',2,'Marker','o','MarkerSize',10)
hold all
plot(xax,T_CA,'color','k','LineWidth',2,'Marker','o','MarkerSize',10)
plot(xax,T_wall(1:6),'color','r', 'LineStyle', '--','LineWidth',2,'Marker','o','MarkerSize',10)
plot(xax,T_wall(7:12),'color','b', 'LineStyle', '--','LineWidth',2,'Marker','o','MarkerSize',10)

hold off

% set('XLim',[1 6])
set(gca, 'XTick', [1 2 3 4 5 6])


xlabel('Element Index','FontSize',20)
ylabel('Temperature (K)','FontSize',20)
legend('Flue Gas', 'Cold Air', 'Surface 1', 'Surface 2')

set(gca,'Fontsize',20)




