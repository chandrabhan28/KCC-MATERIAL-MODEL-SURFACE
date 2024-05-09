a = [18.11, 0.539, 0.001143]; % values of three constants from the ratio p/3 
fc = 45.2; ft = 3.2;
alpha = 1.15; %fro determing the zeta(4) see page 300;
zeta_four = alpha*fc/( a(1) + 2*alpha*fc/(3*a(2)+2*a(3)*alpha*fc) );  % zeta(4) because of complex relation determined seprately
zeta = [0.5  0.5+1.5*(ft/fc) zeta_four  0.753  1]; % ratio of tensile to compressive merdian
jj=[0 20 40 60]; %equal_division_zeta = zeros(1:80);
% this is for distributing the calculated zeta in 80 equal parts like 20
% parts inbetween z(1) and z(2). This is used to determine the william
% warnke surface at different pressure
for ii=1:4
    zeta_division(1+jj(ii):20*ii) = linspace(zeta(ii),zeta(ii+1),20);
end
%% development of compressive meredian for multiplying with the radius if william warnke surface
number = length(zeta_division);   % number of points along the zeta axis to plot the surface
pressure_z_axis = 10;   % this the end pressure axis (zeta axis)
pressure = linspace(-ft,pressure_z_axis,number); % for plotting the pressure axis or z- axis
delta_sigma = ones(1,number); ii =1;
% determing the number of pressure points between preswure 0 and fc/3
negative = sum(pressure<0); % this is to detemine the number of points below 0 or  negative pressure in preesure variable
result = sum(((pressure< fc/3) & (pressure >0))>0); % determiming the number of points inbetween 0 and fc/3;
zeta_interpolt = linspace(zeta(1),zeta(2),result);
while(ii<number+1)
    if pressure(ii)<fc/3   % for determining the compressive meredian surface of ultinmate strength of concrete
        delta_sigma(ii) = 1.5.*(pressure(ii)+ft); 
        if pressure(ii)<0
            delta_sigma(ii) = delta_sigma(ii) /0.5; ii=ii+1;   % dividing tensile meredian with zeta (tesnile to compressive meredian)
        else
            delta_sigma(ii) =  delta_sigma(ii)./zeta_interpolt(ii-negative) ;ii=ii+1;    % dividing tensile meredian with zeta (tesnile to compressive meredian)
        end
    else 
        delta_sigma(ii) = a(1) + pressure(ii) /(a(2) + a(3)*pressure(ii)); ii=ii+1;
    end
end
%% the surface of KCC model, for unit compressive meredian, later it has to be multipiled with compression meredian to obtain the actual surfcae
costheta = cosd(0:1:60);
for ii=1:80
    % zeta = 0.5; % lookat this while running the complete prograam.
    % zeta_division is used here
    radius(ii,1:61) = (2 * (1-zeta_division(ii)^2).*costheta + (2*zeta_division(ii)-1)*sqrt(4*(1-zeta_division(ii)^2).*costheta.^2 + 5*zeta_division(ii)^2 -4*zeta_division(ii)))./(4*(1-zeta_division(ii)^2).*costheta.^2 + (1-2*zeta_division(ii))^2); 
    thetaa = [0:60, 60:120, 120:180, 180:240, 240:300, 300:360].*(pi/180);
    v(ii,1:366) = [radius(ii,1:61), fliplr(radius(ii,1:61)), radius(ii,1:61), fliplr(radius(ii,1:61)), radius(ii,1:61), fliplr(radius(ii,1:61))];  % converting to a complete 360 degree
    [x(ii,1:366), y(ii,1:366)] = pol2cart(thetaa, delta_sigma(ii).*v(ii,1:366)); % multiply delta_sigma with radius of william warnke surface
end
%plot(x,y,"LineWidth",2);  hold on; daspect([1 1 1]);
% preparing the z for plotting the surface making the matrix dimenison consistent
z =ones(number,366);
pressure = linspace(-ft,pressure_z_axis,number); % for plotting the pressure axis or z- axis
for kk=1:80
    z(kk,1:366) = z(kk,1:366)*pressure(kk);
end
%%
figure;
s = surf(x,y,z);
xlabel('x axis'); ylabel('y-axis'); zlabel('z axis');
s.EdgeColor = 'none'; s.FaceAlpha = 0.5;  daspect([1 1 1]);
%% graph properties

% legend({'maximum envelope'},'FontSize',25);
%     set(gca,'FontSize',15);
%     ylabel('delta sigma (MPa)');
%     xlabel('pressure (MPa)');
%     grid on;
