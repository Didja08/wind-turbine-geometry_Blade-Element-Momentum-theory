% Define wind turbine geometry
R = 40; % rotor radius [m]
N = 3; % number of blades
B = 3; % number of blade sections
r = linspace(0,R,50); % radial coordinates
theta = linspace(0,2*pi,100); % angular coordinates
[Rgrid,Theta] = meshgrid(r,theta); % meshgrid of radial and angular coordinates
[X,Y] = pol2cart(Theta,Rgrid); % convert to Cartesian coordinates
Z = zeros(size(X)); % height of blade sections (assumed to be flat)

% Define flow conditions
U = 10; % wind speed [m/s]
rho = 1.225; % air density [kg/m^3]
mu = 1.7894e-5; % air viscosity [kg/m/s]
M = 0.4;  %Mach number
omega=2.7; 
% Define airfoil properties (e.g. using XFOIL)
alpha = [-15 -10 -5 0 5 10 15]; % angle of attack [deg]
cl = [0.0307 0.2479 0.5432 0.8072 1.0410 1.2591 1.4767]; % lift coefficient
cd = [0.0500 0.0163 0.0101 0.0090 0.0085 0.0082 0.0082]; % drag coefficient

% Initialize variables
chord = zeros(N,B);
twist = zeros(N,B);
r = zeros(N,B+1);
dr = zeros(N,B);
alpha_deg = zeros(N,B);
cl_interp = zeros(N,B);
cd_interp = zeros(N,B);
lift = zeros(N,B);
drag = zeros(N,B);
force = zeros(N,B,3); % force vector [F_x F_y F_z]
Cn = zeros(N,B); % normal force coefficient
Ct = zeros(N,B); % tangential force coefficient
a = zeros(N,B); % axial induction factor
aprime = zeros(N,B); % tangential induction factor
alpha_rad = zeros(N,B); % angle of attack in radians
phi = zeros(N,B); % inflow angle
dphi = zeros(N,B); % incremental inflow angle
Fp = zeros(N,B); % pressure force
Fv = zeros(N,B); % viscous force
Fm = zeros(N,B); % momentum force
P = zeros(size(X)); % sound pressure field

% Loop over blade sections
for j = 1:B
    % Define chord and twist
    for i = 1:N
        chord(i,j) = R/N/7; % chord length
        twist(i,j) = 12-6*(r(j)/(R/N)); % twist angle
    end
    
    % Define radial positions
    for i = 1:N
        r(i,j+1) = (j-1)*R/N + r(j); % radial position
        dr(i,j) = r(i,j+1) - r(i,j); % incremental radius
    end
    
    % Interpolate airfoil properties
    for i = 1:N
        alpha_deg(i,j) = atan2d(U*(1-a(i,j)),(r(i,j)+0.5*chord(i,j))*omega); % angle of attack
        cl_interp(i,j) = interp1(alpha,cl,alpha_deg(i,j),'linear','extrap');                                                                            cd_interp(i,j) = interp1(alpha,cd,alpha_deg(i,j),'linear','extrap');
end

% Calculate lift and drag
for i = 1:N
    alpha_rad(i,j) = alpha_deg(i,j)*pi/180; % angle of attack in radians
    lift(i,j) = 0.5*rho*U^2*chord(i,j)*cl_interp(i,j); % lift force
    drag(i,j) = 0.5*rho*U^2*chord(i,j)*cd_interp(i,j); % drag force
    force(i,j,:) = [-drag(i,j)*sin(alpha_rad(i,j)) lift(i,j) drag(i,j)*cos(alpha_rad(i,j))]; % force vector
    Cn(i,j) = lift(i,j)/(0.5*rho*U^2*chord(i,j)); % normal force coefficient
    Ct(i,j) = drag(i,j)/(0.5*rho*U^2*chord(i,j)); % tangential force coefficient
end

% Calculate axial and tangential induction factors
for i = 1:N
    a(i,j) = 1/(4*sin(theta(j)/2)^2/(Cn(i,j)/Ct(i,j))+1); % axial induction factor
    aprime(i,j) = 4*sin(theta(j)/2)*cos(theta(j)/2)/(Cn(i,j)/Ct(i,j)*sin(theta(j)/2)^2+4*sin(theta(j)/2)*cos(theta(j)/2)); % tangential induction factor
end

% Calculate incremental inflow angle
for i = 1:N
    phi(i,j) = atan2(U*(1-a(i,j)),(r(i,j)+0.5*chord(i,j))*(1+aprime(i,j))*omega); % inflow angle
    if j == 1
        dphi(i,j) = 0;
    else
        dphi(i,j) = phi(i,j) - phi(i,j-1);
    end
end

% Calculate pressure force
for i = 1:N
    Fp(i,j) = 0.5*rho*U^2*chord(i,j)*(Cn(i,j)*cos(phi(i,j))+Ct(i,j)*sin(phi(i,j))); % pressure force
end

% Calculate viscous force
for i = 1:N
Re = rho*U*chord(i,j)/mu; % Reynolds number
if Re <= 1e5
Cf = 1.328/sqrt(Re); % laminar flow
else
Cf = 0.455/((log10(Re))^2.58*(1+0.144*M^2)^0.65); % turbulent flow
end
tau_w = 0.5*rho*U^2.*Cf*chord(i,j); % wall shear stress
Fv(i,j) = tau_w*chord(i,j); % viscous force
end

% Calculate momentum force
for i = 1:N
Fm(i,j) = -0.5*rho*U*dphi(i,j)*r(i,j)*chord(i,j)*(1-a(i,j))^2; % momentum force
end

%P = zeros(length(r), length(theta));

% Calculate sound pressure field
for i = 1:N
for k = 1:length(theta)
r_ = r(i,j) + (0.5*dr(i,j));
theta_ = theta(k);
x_ = r*cos(theta_);
y_ = r*sin(theta_);
z_ = 0;
P(i,k) = P(i,k) + 1/(4*pi*rho)*Fp(i,j)/norm([X-Y-Z]-[x_ y_ z_]);
P(i,k) = P(i,k) + 1/(4*pi*rho)*Fv(i,j)/norm([X-Y-Z]-[x_ y_ z_]);
P(i,k) = P(i,k) + 1/(4*pi*rho)*Fm(i,j)/norm([X-Y-Z]-[x_ y_ z_]);
end
end

% Plot pressure field
figure
contourf(X,Y,P)
xlabel('x [m]')
ylabel('y [m]')
title('Sound pressure field [Pa]')
colorbar

% Plot blade geometry
figure
for i = 1:N
plot3(X(i,:),Y(i,:),Z(i,:),'k')
hold on
end
axis equal
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
title('Blade geometry')
view(3)
end


