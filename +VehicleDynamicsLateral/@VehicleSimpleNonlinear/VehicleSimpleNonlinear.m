classdef VehicleSimpleNonlinear < VehicleDynamicsLateral.VehicleSimple
    % VehicleSimpleNonlinear Nonlinear simple vehicle model.
    %
    % It inherits properties from VehicleSimple.

    methods
        % Constructor
        function self = VehicleSimpleNonlinear()
            self.mF0 = 700;
            self.mR0 = 600;
            self.IT = 10000;
            self.lT = 3.5;
            self.nF = 2;
            self.nR = 2;
            self.wT = 2;
            self.muy = .8;
            self.deltaf = 0;
            self.Fxf = 0;
            self.Fxr = 0;
        end

        %% Model
        
        function dx = Model(self, t, states,tspan)
            % Data
            m = self.mT;
            I = self.IT;
            a = self.a;
            b = self.b;
            nF = self.nF;
            nR = self.nR;
            muy = self.muy;


            g = 9.81;                 % Gravity [m/s^2]

            FzF = self.mF0 * g;       % Vertical load @ F [N]
            FzR = self.mR0 * g;       % Vertical load @ R [N]

            % Estados
            X = states(1);
            Y = states(2);
            PSI = states(3);
            v = states(4);
            ALPHAT = states(5);
            dPSI = states(6);


            if isa(self.deltaf,'function_handle')
                DELTA = self.deltaf([X;Y;PSI;v;ALPHAT;dPSI],t);
            elseif length(self.deltaf)>1
                DELTA = interp1(tspan,self.deltaf,t);
            else
                DELTA = self.deltaf;
            end


            % Slip angles
            ALPHAF = atan2((v * sin(ALPHAT) + a * dPSI), (v * cos(ALPHAT))) - DELTA; % Dianteiro
            ALPHAR = atan2((v * sin(ALPHAT) - b * dPSI), (v * cos(ALPHAT)));         % Traseiro

            % Longitudinal forces
            if isa(self.Fxf,'function_handle')
                FxF = self.Fxf([X;Y;PSI;v;ALPHAT;dPSI],t);
            elseif length(self.Fxf)>1
                FxF = interp1(tspan,self.Fxf,t);
            else
                FxF = self.Fxf;
            end

            if isa(self.Fxr,'function_handle')
                FxR = self.Fxr([X;Y;PSI;v;ALPHAT;dPSI],t);
            elseif length(self.Fxr)>1
                FxR = interp1(tspan,self.Fxr,t);
            else
                FxR = self.Fxr;
            end

            % Characteristic curve
            FyF = nF * self.tire.Characteristic(ALPHAF, FzF/nF, muy);
            FyR = nR * self.tire.Characteristic(ALPHAR, FzR/nR, muy);

            % Equations of motion
            dx(1,1) = v * cos(ALPHAT + PSI); % X
            dx(2,1) = v * sin(ALPHAT + PSI); % Y
            dx(3,1) = dPSI; % dPSI
            dx(4,1) = (FxF * cos(ALPHAT - DELTA) + FxR * cos(ALPHAT) + FyF * sin(ALPHAT - DELTA) + FyR * sin(ALPHAT))/(m);
            dx(5,1) = ( - FxF * sin(ALPHAT - DELTA) - FxR * sin(ALPHAT) + FyF * cos(ALPHAT - DELTA) + FyR * cos(ALPHAT) - m * v * dPSI) / (m * v);
            dx(6,1) = (FxF * a * sin(DELTA) + FyF * a * cos(DELTA) - FyR * b) / I;

        end
        
        %% Friction Estimation using Tyre model
        
        function mu_estimate = GetCoefficient(self, v, ALPHAT, dPSI)
         
            a = self.a;
            b = self.b;
            nF = self.nF;
            nR = self.nR;
            muy = self.muy;


            g = 9.81;                 % Gravity [m/s^2]

            FzF = self.mF0 * g;       % Vertical load @ F [N]
            FzR = self.mR0 * g;       % Vertical load @ R [N]

            % Estados
           
            if isa(self.deltaf,'function_handle')
                DELTA = self.deltaf([X;Y;PSI;v;ALPHAT;dPSI],t);
            elseif length(self.deltaf)>1
                DELTA = interp1(tspan,self.deltaf,t);
            else
                DELTA = self.deltaf;
            end


            % Slip angles
            ALPHAF = atan2((v * sin(ALPHAT) + a * dPSI), (v * cos(ALPHAT))) - DELTA; % Dianteiro
            ALPHAR = atan2((v * sin(ALPHAT) - b * dPSI), (v * cos(ALPHAT)));         % Traseiro

       
            % Characteristic curve
            FyF = nF * self.tire.Characteristic(ALPHAF, FzF/nF, muy);
            FyR = nR * self.tire.Characteristic(ALPHAR, FzR/nR, muy);

            mu_estimate = [FyF/FzF, FyR/FzR];

        end
        
        %% Friction Estimation using Accelaration Values
       
        function mu_estimate = GetCoefficientEst(self, v, ALPHAT, PSI, dPSI, dvt, dalphat)
            % Data
            m = self.mT;
            
            g = 9.81;                 % Gravity [m/s^2]

            FzF = self.mF0 * g;       % Vertical load @ F [N]
            FzR = self.mR0 * g;       % Vertical load @ R [N]

            % Estados
           
            if isa(self.deltaf,'function_handle')
                DELTA = self.deltaf([X;Y;PSI;v;ALPHAT;dPSI],t);
            elseif length(self.deltaf)>1
                DELTA = interp1(tspan,self.deltaf,t);
            else
                DELTA = self.deltaf;
            end


            % Longitudinal forces
            if isa(self.Fxf,'function_handle')
                FxF = self.Fxf([0;0;0;v;0;0],0);
            elseif length(self.Fxf)>1
                FxF = interp1(tspan,self.Fxf,t);
            else
                FxF = self.Fxf;
            end

            if isa(self.Fxr,'function_handle')
                FxR = self.Fxr([0;0;0;v;0;0],0);
            elseif length(self.Fxr)>1
                FxR = interp1(tspan,self.Fxr,t);
            else
                FxR = self.Fxr;
            end

            % Calculating Lateral forces from Lateral and Longitudanal
            % accelaration values
            ax = cos(ALPHAT)*dvt - v*sin(ALPHAT)*dalphat;
            ay = sin(ALPHAT)*dvt + v*cos(ALPHAT)*dalphat;
            
            at = [ax;ay];
            
            A = [-sin(PSI+DELTA), -sin(PSI); 
                  cos(PSI+DELTA), cos(PSI)];
            
            B = [cos(PSI+DELTA), cos(PSI); 
                 sin(PSI+DELTA), sin(PSI)];
            
            
            Fy_est = m*A\at - A\B*[FxF; FxR];
            
            FyF = Fy_est(1);
            FyR = Fy_est(2);

            mu_estimate = [FyF/FzF, FyR/FzR];

        end
        
    end
end

%% See Also
%
% <../../index.html Home>
%
