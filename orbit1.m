clc;
close all;
clear all;

plotstyle = {'bs','-go','-r^','-c','-mx','-g>','-kd',};

%constants
G = 6.673 * 10^-11;
month = 60*60*24*30;
week = 60*60*24*7;
mass_sun = 1.9891*10^30;
%planetary constants
mass_earth = 5.97219*10^24;
e_earth = 0.017;
a_earth = 149598023*1000;
mass_mars = 6.39*10^23;
e_mars = 0.0935;
a_mars = 2.28*10^8 * 1000;
mass_jupiter = 1.898*10^27;
e_jupiter = 0.0487;
a_jupiter = 7.78*10^8*1000;
mass_saturn = 5.683*10^26;
e_saturn = 0.052;
a_saturn = 14.3353*10^8*1000;
mass_uranus = 8.681*10^25;
e_uranus = 0.04717;
a_uranus = 28.67*10^8*1000;
mass_neptune = 1.024*10^26;
e_neptune = 0.0097;
a_neptune = 45*10^8*1000;
mass_venus = 4.867*10^24;
e_venus = 0.007;
a_venus = 1.08 * 10^8*1000;
%period of earth
P_earth = sqrt((4*pi*pi*a_earth^3)/(G*(mass_sun+mass_earth)));
%plot orbits
[xvenus,yvenus] = orbit(mass_sun,mass_venus,e_venus,a_venus,(15*pi)/180);
[xearth,yearth] =orbit(mass_sun,mass_earth,e_earth,a_earth,0);
[xmars,ymars] = orbit(mass_sun,mass_mars,e_mars,a_mars,(200*pi)/180);
[xjupiter,yjupiter] = orbit(mass_sun,mass_jupiter,e_jupiter,a_jupiter,(25*pi)/180);
%[xsaturn,ysaturn] = orbit(mass_sun,mass_saturn,e_saturn,a_saturn,(315*pi)/180);
%[xuranus,yuranus] = orbit(mass_sun,mass_uranus,e_uranus,a_uranus,(60*pi)/180);
%[xneptune,yneptune] = orbit(mass_sun,mass_neptune,e_neptune,a_neptune,(345*pi)/180);
%x = [xvenus;xearth;xmars;xjupiter;xsaturn;xuranus;xneptune];
%y = [yvenus;yearth;ymars;yjupiter;ysaturn;yuranus;yneptune];
x = [xvenus;xearth;xmars;xjupiter];
y = [yvenus;yearth;ymars;yjupiter];


clf; hold on
ph = zeros(size(x,1));
for i=1:size(x,1)
    plot3(x(:,1),y(:,1),zeros(size(x,1)),'bs');
end
%orbit plotting function
function[xcord,ycord] = orbit(mass_star,mass_planet, eccentricity, semimajor,theta0)
steps = 100;
days = 3;
G = 6.673 * 10^-11;
theta = zeros(steps);
r = zeros(steps);
xcord = zeros(steps);
ycord = zeros(steps);
MeanAnomaly = zeros(steps);
E = zeros(steps);
for i = 1:1:steps
    %initial conditions
    if(i==1)
        theta(i) = theta0;
        r(i) = (semimajor*(1-eccentricity^2))/(1+eccentricity*cos(theta0));
    end
    %iteration 
    
    if(i>0)
    elapsed_time =(i/2)*60*60*24*days;
    MeanAnomaly(i+1) = elapsed_time*sqrt((G*(mass_star+mass_planet))/(semimajor^3));
    for j = 1:1:20
       if(j==1)
          E(i,j) = MeanAnomaly(i+1);
       end
       E(i,j+1) = E(i,j) - (MeanAnomaly(i+1) - E(i,j) + eccentricity*sin(E(i,j)))/(eccentricity*cos(E(i,j))-1) ;
    end
   
    theta(i+1) = 2*atan((sqrt((1+eccentricity)/(1-eccentricity)))*(tan(E(i,20))/2))+theta(1);
    r(i+1) = (semimajor*(1-eccentricity^2))/(1+eccentricity*cos(theta(i+1)));
   
    xcord(i+1) = semimajor*cos(theta(i+1));
    ycord(i+1) = semimajor*sin(theta(i+1));
    end
end
end

