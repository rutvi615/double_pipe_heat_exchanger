clc;
load('vfs.mat')
viewF


fflow=0.5;
mdotmax = 1.5;
mdotmin = 0.07;
mdot = 0.0693 + 0.025 * (fflow - mdotmin) / (mdotmax - mdotmin);  % mdot bound constraint to keep hot air outlet T < 533.15 K to control NOx emissions

T_CA_in = 300;  % Temperature of the inlet cold air
T_amb = 300;    % Ambient temperature
T_exh = 1004;   % Exhaust temperature
% Radiation constants
sig = 5.672e-8;  % W/(m^2-K^4)
absorp = 0.9;
Lmb = 0.95;
eps = 0.9;

Cp_FG = 1.37e3;  % J/kg-K
mdot_FG = mdot;
h_FG = 150;      % W/(m^2-K)    % Convection of flue gas
mdot_air = mdot_FG * 0.18;
Cp_air = 1e3;    % J/kg-K
h_air = 150;     % W/(m^2-K)    % Heat transfer coefficient of air inside pipe
h_amb = 100;     % W/(m^2-K)    % Heat transfer coefficient from insulation to ambient air
k_ins = 0.15;

r2 = 0.07;        % Radius of surface 2
r1 = 0.04;        % Radius of surface 1
L_ins = 0.05;     % Thickness of insulation
r_ins = r1 + L_ins;
hCyl = 0.6;       % Height of the cylinder
h = 0.1;          % Height of each discretized surface
A1 = 2 * pi * r1 * h;  % Area of each element in surface 1
A2 = 2 * pi * r2 * h ; % Area of each element in surface 2

R_amb = log(r_ins/r1) / (2 * pi * k_ins * h) + 1 / (h_amb * A1);  % Resistance to heat transfer from surface 1 to ambient air

JLHS = zeros(12,12);  % Radiosity matrix
JRHS = zeros(12,1);  % Right-hand side vector
JRec = zeros(12,1);  % Radiosity vector
JLHS22 = zeros(12,12); 
JRHS22 = zeros(12,1);
JRec22 = zeros(6,1);

% Initialize wall and gas temperatures
T_w = [900 875 850 825 800 775; 600 620 640 660 680 700];  % Wall temperatures (surface, element)
T_FG = 300*ones(6,1); % temperatures of flue gas 

T_CA = 300*ones(6,1); % temperatures of cold air
T_FG(1) = T_exh; 
% T_FG = [T_exh; 300; 300];  % Flue gas temperatures
% T_CA = [T_CA_in; 300; 300];  % Cold air temperatures

T_wall = T_w';
T_wall = T_wall(:);

T_wall_store = T_wall;
T_wall_arr=T_wall(8:12)';  
T_FG_store = T_FG;
T_CA_store = T_CA;

    for i = 1:6
        for j = 1:12
            if i == j
                JLHS(i,j) = A1 * sum(viewF(i,1:12)) + A1 * (eps) / (1 - eps)
            else
                JLHS(i,j) = -A1 * viewF(i,j);
            end
        end
    end
    
    for i = 7:12
        for j = 1:12
            if i == j
                JLHS(i,j) = A2 * sum(viewF(i,1:12)) + A2 * (eps) / (1 - eps)
            else
                
               if (j<7)
                JLHS(i,j) = -A2*viewF(i,j);       
               else
                JLHS(i,j) = 0;
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


for k = 1:100  % 50 iterations

    % Update flue gas temperatures
    for i = 1:5
        T_FG(i+1) = T_FG(i) - (1/(mdot_FG*Cp_FG))*(h_FG*A1*(T_FG(i)-T_w(1,i)) + h_FG*A2*(T_FG(i)-T_w(2,i)));
   end

    % Update cold air temperatures
    for i = 1:5
        T_CA(i+1) = T_CA(i) + (1/(mdot_air*Cp_air))*(h_air*A2*(T_w(2,i)-T_CA(i)));
    end

    T_wall = T_w';
T_wall = T_wall(:);

    

    % Calculate JRHS12 for surfaces 1 and 2
    for i = 1:6
        JRHS(i) = sig * T_wall(i)^4 / ((1 - eps) / (eps * A1));
    end
    for i = 7:12
        JRHS(i) = sig * T_wall(i)^4 / ((1 - eps) / (eps * A2));
    end

    % Solve for JRec (radiosity)
    JRec12 = JLHS\JRHS;
    
    JRec(1:12) = JRec12;

    

    % Update wall temperatures for surface 1
%     for i=1:3
%     radFlux=0;
%     for j=1:6
%     if (i==j)
%         radFlux = radFlux+A1*sum(viewF(i,1:6))*JRec(j);
%     else
%         radFlux = radFlux-A1*viewF(i,j)*JRec(j);
%     end
%     end
% T_wall(i) = -(radFlux-(R_amb/R_amb)-(h_FG*A1*T_FG(i)))/(h_FG*A1+1/R_amb)
% end
    for i = 1:6
        radFlux = 0;
        for j = 1:12
            if i == j
                radFlux = radFlux + A1 * sum(viewF(i,1:12)) * JRec(j);
            else
                radFlux = radFlux - A1 * viewF(i,j) * JRec(j);
            end
        end
        T_wall(i) = -(radFlux - (T_amb / R_amb) - (h_FG * A1 * T_FG(i))) / (h_FG * A1 + 1 / R_amb);
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




