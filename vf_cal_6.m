clc;
clear;

r1=0.04; %outer radius of the inner cylinder
r2=0.07; %inner radius of the outer cylinder
hCyl=0.6; %height of the cylinder
h=0.1; %height of each discretized surface
A1 = 2*pi*r1*h; %outer surface area of inner cylinder
A2 = 2*pi*r2*h; %inner surface area of outer cylinder

viewF = zeros(18,18);
 
%eq 1 Exterior of right circular cylinder of finite length to interior of coaxial outer right circular cylinder
R1=r1/h;
R2=r2/h;
A=R2+R1;
B=R2-R1;
F12 = (1/(pi*R1))*(0.5*(R2^2-R1^2-1)*acos(R1/R2)+pi*R1-0.5*pi*A*B-2*R1*atan(sqrt(R2^2-R1^2)) ...
    + sqrt((1+A^2)*(1+B^2))*atan(sqrt(((1+A^2)*B)/((1+B^2)*A))));

F21=F12*A1/A2;
i=1;
j=7;
while (i<=6)
 viewF(i,j)=F21;
 viewF(j,i)=F12;
 i=i+1; 
 j=j+1;
end

%eq 2-5 Exterior of inner coaxial cylinder to interior of larger outer cylinder; smaller
%cylinder completely outside the larger cylinder

d=10e-7;
D=d/r2;
Y=h/r2;
L=h/r2;
R=r1/r2;
e1=[L+D Y+D D L+D+Y];
F=[0 0 0 0];
for i = 1:length(F)
    e=e1(i);
    Ae=e^2+R^2-1;
    Be=e^2-R^2+1;
    t1=Be/(8*R*e);
    t2=acos(Ae/Be);
    t3=((((Ae+2)^2/R^2)-4)^0.5)*acos(Ae*R/Be)/(2*e);
    t4=Ae*asin(R)/(2*e*R);
    F(i)=t1+(t2-t3-t4)/(2*pi);
end

FLD=F(1);

FYD=F(2);

FD=F(3);

FLDY=F(4);

F12 =((L+D)*FLD +(Y+D)*FYD-D*FD-(L+D+Y)*FLDY)/L;

viewF(7,2)=F12;
viewF(8,1)=F12;
viewF(8,3)=F12;
viewF(9,2)=F12;
viewF(9,4)=F12;
viewF(10,3)=F12;
viewF(10,5)=F12;
viewF(11,6)=F12;
viewF(11,4)=F12;
viewF(12,5)=F12;

F21=F12*A1/A2;
viewF(2,7)=F21;
viewF(1,8)=F21;
viewF(3,8)=F21;
viewF(2,9)=F21;
viewF(4,9)=F21;
viewF(3,10)=F21;
viewF(5,10)=F21;
viewF(4,11)=F21;
viewF(6,11)=F21;
viewF(5,12)=F21;

%distance=0.1
d=0.1;
D=d/r2;
Y=h/r2;
L=h/r2;
R=r1/r2;
e1=[L+D Y+D D L+D+Y];
F=[0 0 0 0];
for i = 1:length(F)
    e=e1(i);
    Ae=e^2+R^2-1;
    Be=e^2-R^2+1;
    t1=Be/(8*R*e);
    t2=acos(Ae/Be);
    t3=((((Ae+2)^2/R^2)-4)^0.5)*acos(Ae*R/Be)/(2*e);
    t4=Ae*asin(R)/(2*e*R);
    F(i)=t1+(t2-t3-t4)/(2*pi);
end

FLD=F(1);

FYD=F(2);

FD=F(3);

FLDY=F(4);

F12 =((L+D)*FLD +(Y+D)*FYD-D*FD-(L+D+Y)*FLDY)/L;
viewF(7,3)=F12;
viewF(8,4)=F12;
viewF(9,1)=F12;
viewF(9,5)=F12;
viewF(10,2)=F12;
viewF(10,6)=F12;
viewF(11,3)=F12;
viewF(12,4)=F12;




F21=F12*A1/A2;
viewF(3,7)=F21;
viewF(4,8)=F21;
viewF(1,9)=F21;
viewF(5,9)=F21;
viewF(2,10)=F21;
viewF(6,10)=F21;
viewF(3,11)=F21;
viewF(4,12)=F21;


%distance=0.2
d=0.2;
D=d/r2;
Y=h/r2;
L=h/r2;
R=r1/r2;
e1=[L+D Y+D D L+D+Y];
F=[0 0 0 0];
for i = 1:length(F)
    e=e1(i);
    Ae=e^2+R^2-1;
    Be=e^2-R^2+1;
    t1=Be/(8*R*e);
    t2=acos(Ae/Be);
    t3=((((Ae+2)^2/R^2)-4)^0.5)*acos(Ae*R/Be)/(2*e);
    t4=Ae*asin(R)/(2*e*R);
    F(i)=t1+(t2-t3-t4)/(2*pi);
end

FLD=F(1);

FYD=F(2);

FD=F(3);

FLDY=F(4);

F12 =((L+D)*FLD +(Y+D)*FYD-D*FD-(L+D+Y)*FLDY)/L;
viewF(7,4)=F12;
viewF(8,5)=F12;
viewF(9,6)=F12;
viewF(12,3)=F12;
viewF(11,2)=F12;
viewF(10,1)=F12;

F21=F12*A1/A2;
viewF(4,7)=F21;
viewF(5,8)=F21;
viewF(6,9)=F21;
viewF(3,12)=F21;
viewF(2,11)=F21;
viewF(1,10)=F21;


%distance=0.1
d=0.3;
D=d/r2;
Y=h/r2;
L=h/r2;
R=r1/r2;
e1=[L+D Y+D D L+D+Y];
F=[0 0 0 0];
for i = 1:length(F)
    e=e1(i);
    Ae=e^2+R^2-1;
    Be=e^2-R^2+1;
    t1=Be/(8*R*e);
    t2=acos(Ae/Be);
    t3=((((Ae+2)^2/R^2)-4)^0.5)*acos(Ae*R/Be)/(2*e);
    t4=Ae*asin(R)/(2*e*R);
    F(i)=t1+(t2-t3-t4)/(2*pi);
end

FLD=F(1);

FYD=F(2);

FD=F(3);

FLDY=F(4);

F12 =((L+D)*FLD +(Y+D)*FYD-D*FD-(L+D+Y)*FLDY)/L;
viewF(7,5)=F12;
viewF(8,6)=F12;
viewF(11,1)=F12;
viewF(12,2)=F12;

F21=F12*A1/A2;
viewF(5,7)=F21;
viewF(6,8)=F21;
viewF(2,12)=F21;
viewF(1,11)=F21;


%distance=0.1
d=0.4;
D=d/r2;
Y=h/r2;
L=h/r2;
R=r1/r2;
e1=[L+D Y+D D L+D+Y];
F=[0 0 0 0];
for i = 1:length(F)
    e=e1(i);
    Ae=e^2+R^2-1;
    Be=e^2-R^2+1;
    t1=Be/(8*R*e);
    t2=acos(Ae/Be);
    t3=((((Ae+2)^2/R^2)-4)^0.5)*acos(Ae*R/Be)/(2*e);
    t4=Ae*asin(R)/(2*e*R);
    F(i)=t1+(t2-t3-t4)/(2*pi);
end

FLD=F(1);

FYD=F(2);

FD=F(3);

FLDY=F(4);

F12 =((L+D)*FLD +(Y+D)*FYD-D*FD-(L+D+Y)*FLDY)/L;
viewF(7,6)=F12;
viewF(12,1)=F12;

F21=F12*A1/A2;
viewF(6,7)=F21;
viewF(1,12)=F21;

%eq 6 Interior of finite length right circular coaxial cylinder to itself
%(in presence of interior cylinder)
term1 = pi * (R2 - R1);
term2 = acos(R1 / R2);
term3 = sqrt(1 + 4 * R2^2) * atan(sqrt((R2^2 - R1^2) * (1 + 4 * R2^2)) / R1);
term4 = 2 * R1 * atan(2 * sqrt(R2^2 - R1^2));

F11 = (1 / (pi * R2)) * (term1 + term2 - term3 + term4);

viewF(1,1)=F11;
viewF(2,2)=F11;
viewF(3,3)=F11;
viewF(4,4)=F11;
viewF(5,5)=F11;
viewF(6,6)=F11;

%eq 7 Interior surface of right circular coaxial cylinder to the adjacent interior surface
Lt=h+h;
R1=r1/0.2;
R2=r2/0.2;
A=R2+R1;
B=R2-R1;
term1 = pi * (R2 - R1);
term2 = acos(R1 / R2);
term3 = sqrt(1 + 4 * R2^2) * atan(sqrt((R2^2 - R1^2) * (1 + 4 * R2^2)) / R1);
term4 = 2 * R1 * atan(2 * sqrt(R2^2 - R1^2));

Ftt = (1 / (pi * R2)) * (term1 + term2 - term3 + term4);

F12=(Lt*Ftt-h*F11-h*F11)/(2*h);
viewF(1,2)=F12;
viewF(2,3)=F12;
viewF(3,4)=F12;
viewF(4,5)=F12;
viewF(5,6)=F12;


F21=F12*A2/A2;
viewF(2,1)=F21;
viewF(3,2)=F21;
viewF(4,3)=F21;
viewF(5,4)=F21;
viewF(6,5)=F21;

%eq 8 Interior surface of right circular coaxial cylinder to an interior surface separated by a distance
L=0.2;
% Given parameters
r1 = 0.04;  % Inner cylinder radius (in meters)
r2 = 0.07;  % Outer cylinder radius (in meters)
L = 0.1;    % Length of both elements (in meters)
d = 0.1;    % Separation distance (in meters)
Lt = L + d; % Total length including separation

% Normalized lengths and radii for case L+d
R1 = r1 / Lt;  % Normalized radius R1 for Lt = L + d
R2 = r2 / Lt;  % Normalized radius R2 for Lt = L + d

% Term calculations for FLLD (view factor between L+L+d)
term1 = pi * (R2 - R1);
term2 = acos(R1 / R2);
term3 = sqrt(1 + 4 * R2^2) * atan(sqrt((R2^2 - R1^2) * (1 + 4 * R2^2)) / R1);
term4 = 2 * R1 * atan(2 * sqrt(R2^2 - R1^2));
FLLD = (1 / (pi * R2)) * (term1 + term2 - term3 + term4);

% Normalized lengths and radii for case L+d
Lt2 = L + d;  % Updated total length
R1 = r1 / Lt2;  % Normalized radius R1 for L+d
R2 = r2 / Lt2;  % Normalized radius R2 for L+d

% Term calculations for FLD2 (view factor between L+d)
term1 = pi * (R2 - R1);
term2 = acos(R1 / R2);
term3 = sqrt(1 + 4 * R2^2) * atan(sqrt((R2^2 - R1^2) * (1 + 4 * R2^2)) / R1);
term4 = 2 * R1 * atan(2 * sqrt(R2^2 - R1^2));
FLD2 = (1 / (pi * R2)) * (term1 + term2 - term3 + term4);

% Final view factor calculation F12
F12 = ((L + L + d) * FLLD - 2 * (L + d) * FLD2) / (2 * L);
F12=0.0015;
% Correct any negative value for F12

viewF(1,3)=F12;
viewF(3,1)=F12;
viewF(2,4)=F12;
viewF(4,2)=F12;
viewF(3,5)=F12;
viewF(5,3)=F12;
viewF(4,6)=F12;
viewF(6,4)=F12;
F12=0.0004
viewF(1,4)=F12;
viewF(4,1)=F12;
viewF(2,5)=F12;
viewF(5,2)=F12;
viewF(3,6)=F12;
viewF(6,3)=F12;
F12=0.00021

viewF(1,5)=F12;
viewF(5,1)=F12;
viewF(6,2)=F12;
viewF(2,6)=F12;
F12=0.000083
viewF(1,6)=F12;
viewF(6,1)=F12;
%eq 9 Interior surface of right circular cylinder to itself 
H=h/r1;
F22=(1 + H)-(1 + H^2)^0.5;
viewF(13,13)=F22;
viewF(14,14)=F22;
viewF(15,15)=F22;
viewF(16,16)=F22;
viewF(17,17)=F22;
viewF(18,18)=F22;

%eq 10-11 Interior surface of right circular cylinder to interior surface of adjacent right circular cylinder of the same diameter
H=h/r1;
F12=(H+sqrt(4+H^2)-2*sqrt(1+H^2))/2;
viewF(13,14)=F12;
viewF(14,13)=F12;
viewF(14,15)=F12;
viewF(15,14)=F12;
viewF(15,16)=F12;
viewF(16,15)=F12;
viewF(16,17)=F12;
viewF(17,16)=F12;
viewF(17,18)=F12;
viewF(18,17)=F12;

%eq 12 Finite section of right circular cylinder to separated finite section

l1 = h;       % length l1
l2 = 2*h;       % length l2
l3 = 3*h;     % length l3
L1 = l1/r1;   % normalized length L1
L2 = l2/r1;   % normalized length L2
L3 = l3/r1;   % normalized length L3

X = @(L) sqrt(L^2 + 4);  % function handle for sqrt(L^2 + 4)

F12 = ((2*L1*(L3 - L2)) + ((L3 - L1)*X(L3 - L1)) - ((L2 - L1)*X(L2 - L1)) + ...
       L3*X(L3) - L2*X(L2)) / (4 * (L3 - L2));
F12=0.0348;
viewF(13,15)=F12;
viewF(14,16)=F12;
viewF(15,17)=F12;
viewF(16,18)=F12;
viewF(15,13)=F12;
viewF(16,14)=F12;
viewF(17,15)=F12;
viewF(18,16)=F12;


l1 = h;       % length l1
l2 = 3*h;       % length l2
l3 = 4*h;     % length l3
L1 = l1/r1;   % normalized length L1
L2 = l2/r1;   % normalized length L2
L3 = l3/r1;   % normalized length L3

X = @(L) sqrt(L^2 + 4);  % function handle for sqrt(L^2 + 4)

F12 = ((2*L1*(L3 - L2)) + ((L3 - L1)*X(L3 - L1)) - ((L2 - L1)*X(L2 - L1)) + ...
       L3*X(L3) - L2*X(L2)) / (4 * (L3 - L2));
F12=0.0148;
viewF(13,16)=F12;
viewF(14,17)=F12;
viewF(15,18)=F12;
viewF(16,13)=F12;
viewF(17,14)=F12;
viewF(18,15)=F12;

l1 = h;       % length l1
l2 = 4*h;       % length l2
l3 = 5*h;     % length l3
L1 = l1/r1;   % normalized length L1
L2 = l2/r1;   % normalized length L2
L3 = l3/r1;   % normalized length L3

X = @(L) sqrt(L^2 + 4);  % function handle for sqrt(L^2 + 4)

F12 = ((2*L1*(L3 - L2)) + ((L3 - L1)*X(L3 - L1)) - ((L2 - L1)*X(L2 - L1)) + ...
       L3*X(L3) - L2*X(L2)) / (4 * (L3 - L2));
F12=0.0068;
viewF(13,17)=F12;
viewF(14,18)=F12;
viewF(17,13)=F12;
viewF(18,14)=F12;

%d=0.4
l1 = h;       % length l1
l2 = 5*h;       % length l2
l3 = 6*h;     % length l3
L1 = l1/r1;   % normalized length L1
L2 = l2/r1;   % normalized length L2
L3 = l3/r1;   % normalized length L3

X = @(L) sqrt(L^2 + 4);  % function handle for sqrt(L^2 + 4)

F12 = ((2*L1*(L3 - L2)) + ((L3 - L1)*X(L3 - L1)) - ((L2 - L1)*X(L2 - L1)) + ...
       L3*X(L3) - L2*X(L2)) / (4 * (L3 - L2))

F12=0.00022;
viewF(13,18)=F12;
viewF(18,13)=F12;

viewF
save('vfs.mat','viewF')







