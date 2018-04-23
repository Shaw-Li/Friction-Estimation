%% Skid pad
% Maneuver in circles of a nonlinear simple vehicle with Pacejka tire model.
%
%% Code start
%

import VehicleDynamicsLateral.*

%% Model and parameters
% Simulation
%

T = 30;                     % Total simulation time [s]
resol = 300;                % Resolution
TSPAN = 0:T/resol:T;        % Time span [s]

%%

TireModel = TirePacejka();

%%
% Vehicle

System = VehicleSimpleNonlinear();

% Defining vehicle parameters
System.mF0 = 700;
System.mR0 = 600;
System.IT = 10000;
System.lT = 3.5;
System.nF = 2;
System.nR = 2;
System.wT = 2;
System.muy = .7;
System.deltaf = 30*pi/180;
System.Fxf = 0;
System.Fxr = @RearVelControl;

System.tire = TireModel;
simulator = Simulator(System, TSPAN);
simulator.V0 = 8.333;

%% Simulation
%

simulator.Simulate();

%% Results
%

% Retrieving states
XT = simulator.XT;
YT = simulator.YT;
PSI = simulator.PSI;
VEL = simulator.VEL;
ALPHAT = simulator.ALPHAT;
dPSI = simulator.dPSI;

mu_est = zeros(300,2);
mu_true = zeros(300,2);
for i=1:300
    mu_true(i,:) = System.GetCoefficient(VEL(i), ALPHAT(i), dPSI(i));
    mu_est(i,:) = System.GetCoefficientEst(VEL(i), ALPHAT(i), PSI(i), dPSI(i), (VEL(i+1)-VEL(i))*10, (ALPHAT(i+1)-ALPHAT(i))*10);
end

% Plot the slip angle as a function of time
figure(1)
hold on ; grid on ; box on
plot(TSPAN,ALPHAT)
ylim([0 .1]);
xlabel('time [s]')
ylabel('Vehicle slip angle [rad]')
title('Slip Angle for Constant Steering')


% Define the moving average window
window = 7;

% Plot the estimated and true value for the front Tyres
figure(2)
hold on ; grid on ; box on
plot(TSPAN(1:300),abs(mu_est(:,1)), '.')
plot(TSPAN(1:(300-window)), abs(movingAverage(mu_est(:,1), window)), 'g-')
plot(TSPAN(1:300), mu_true(:,1), 'b--')
ylim([-1 1]);
legend('Estimated','True')
xlabel('time [s]')
ylabel('Estimated Friction Coefficient Front')
title('Comparison of Estimated vs True Friction Coefficient (Front Axle)')


% Plot the estimated and true value for the Rear Tyres
figure(3)
hold on ; grid on ; box on
plot(TSPAN(1:300),abs(mu_est(:,2)), '.')
plot(TSPAN(1:(300-window)), abs(movingAverage(mu_est(:,2), window)), 'g-')
plot(TSPAN(1:300), mu_true(:,2), 'b--')
ylim([-1 1]);
legend('Estimated', 'True')
xlabel('time [s]')
ylabel('Estimated Friction Coefficient Rear')
title('Comparison of Estimated vs True Friction Coefficient (Rear Axle)')
