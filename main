function main(indx)
% TEP4240 - SYSTEM SIMULATION
% PROJECT
% Developed by: R. Ensalzado, 2016
% MATLAB Version: R2015a

if nargin < 1
    indx = 1;
end

% || Reservoir data
% Pws : reservoir average pressure [kPa]
% Pwf : bottom well pressure       [kPa]
% Pwh : wellhead pressure          [kPa]
% Jo  : Oil productivity index     [Sm3/d/kPa]
% Qsc : Well oil production        [Sm3/d]
% Tw  : process temperature        [C]

Pwh = 6000;
Pws = 20000;
Psm = (Pwh + Pws)/2;
Tw = 45;
Jo = 0.2;
Qo = Pws*Jo;

% || Other constants
% g : gravity acceleration
g = 9.81;

% || Pump data
% Reference characteristic from a ESP (electro submersible pump) with water
% Discharge pressure, in kPa, can be calculated as H*SG*g  
% - H : head     [m]
% - Q : capacity [m3/h]
% - S : centrifugal stages
% - F : flow multiplier

S = 30;
F = 3;
H = S*[24.1; 24.1; 23.8; 23.1; 22.7; 22.3; 21.9; 21.9; 22.0; 22.3; 22.4; ...
         22.5; 22.5; 21.9; 21.1; 19.6; 17.5; 14.2];   
Q = F*[0.0; 10.8; 21.7; 32.5; 43.3; 54.2; 65.0; 75.8; 86.7; 97.5; 108.3; ...
        119.2; 130.0; 140.8; 151.7; 162.5; 173.3; 184.2];

PPW = cESPObj;
PPW.CurveHQ = [H Q];
PPW.BEP = [S*21.1 F*151.7];
PPW.N =  3600;
PPW.Nc = 3600;

% || Fluid models
% Oil phase is model using black oil approach, and tuning the model
% so real fluid properties are calculated.
% - _oil : heavy oil
% - _dil : diluent

SG_oil = 0.9692;
Tsm_oil = [37.8 50.0];
vmu_oil = [906 329];

SG_dil = 0.8720;
Tsm_dil = [38 54];
vmu_dil = [9.8 7.5];

BL0 = cInjectionObj;
BL0.Qo = 1;
BL0.Qd = 1;
BL0.Qw = 0;
BL0.Tb = Tw;
BL0.mrmu = [vmu_oil; vmu_oil]; 
BL0.mrT = [Tsm_oil; Tsm_oil];
BL0.drho = SG_oil*1e3;
BL0.orho = SG_oil*1e3;
BL0.BlendingReference

BL1 = cInjectionObj;
BL1.Qo = 0;
BL1.Qd = 1;
BL1.Qw = 0;
BL1.Tb = Tw;
BL1.mrmu = [vmu_oil; vmu_dil]; 
BL1.mrT = [Tsm_oil; Tsm_dil];
BL1.drho = SG_dil*1e3;
BL1.orho = SG_oil*1e3;
BL1.BlendingReference

BO1 = cBOObj;
BO1.P = conversionP(Psm, 'kPa -> psi');
BO1.T = conversionT(Tw, 'c -> f');
BO1.Tbosc = 851 - 460;
BO1.Tbgsc = 196.47 - 460;
BO1.gSG = 0.64;
BO1.wSG = 1.00;
BO1.oSG = SG_oil;
BO1.GOR = 130;
BO1.WC = 0;
BO1.TuneViscosity(BL0.bmrmu, BL0.bmrT*1.8 + 32)

BO2 = cBOObj;
BO2.P = conversionP(Psm, 'kPa -> psi');
BO2.T = conversionT(Tw, 'c -> f');
BO2.Tbosc = 851 - 460;
BO2.Tbgsc = 196.47 - 460;
BO2.gSG = 0.64;
BO2.wSG = 1.00;
BO2.oSG = SG_oil;
BO2.GOR = 130;
BO2.WC = 0;
BO2.TuneViscosity(BL0.bmrmu, BL0.bmrT*1.8 + 32)

% || Piping data
% Geometrical data to characterize the piping sections
% - di : internal diameter [m]
% - pr : roughness         [m]
% - Ap : sectional area    [m]
% - L  : piping length     [m]

di1 = 5*0.0254;
pr1 = 1.524e-5;
Ap1 = pi*di1^2/4;
L1 = 500;

di2 = 4*0.0254;
pr2 = 1.524e-5;
Ap2 = pi*di2^2/4;
L2 = 500;

% || System solution

out = ode45(@system, [0 1000], 0);

% || Relevant physical variables
% - Qp0 : oil production "consumed" by the source resistance
% - Qp1 : oil production to wellhead
% - Qp2 : oil production to wellhead + diluent

n = length(out.x);
Qp1sc = zeros(1, n);

BO1.P = conversionP(Psm, 'kPa -> psi');
BO1.T = conversionT(Tw, 'c -> f');
MT = BO1.qsc2ac;

rho1 = conversionrho(BO1.orho, 'lb/ft3 -> kg/m3');
vmu1 = BO1.omu*1e-3;
Ip1 = rho1*L1/Ap1;

Qp1 = out.y/Ip1;
Qp2 = out.y/Ip1 + Sfd(out.x);

for j = 1:n
    k = MT\[0; Qp1(j); 0];
    Qp1sc(j) = k(2)*3600*24; 
end

Qp0 = Qo - Qp1sc;

Pwf = fRsw(Qp0).*Qp0;
Ppm = H(1)*g*rho2/1e3 - fRep(Qp2, [rho2/1e3 vmu2]).*Qp2;
dP1 = fRpi(Qp1, 1, [rho1 vmu1]).*Qp1; 
dP2 = fRpi(Qp2, 2, [rho2 vmu2]).*Qp2; 

% || Plotting
% plot 1
% - fluid flow, piping section 1
% - fluid flow, piping section 2
% plot 2
% - pressure (real effort) from ideal sources
% plot 3
% - Pressure drop in piping sections
% - Pressure drop in pump effort source

figure(1)
plot(out.x, Qp1*24*3600, 'ob-',...
     out.x, Qp2*24*3600, 'or-')
title('Flowrate in piping sections')
xlabel('time [s]'); ylabel('flow [m3/d]')
legend('piping section 1', 'piping section 2', 'Location', 'southeast')
grid on

figure(2)
subplot(1, 2, 1)
plot(out.x, Ppm, 'ob-')
title('Real effort from pump (effort source) - Thevenin approximation')
xlabel('time [s]'); ylabel('Pressure [kPa]')
grid on

subplot(1, 2, 2)
plot(out.x, Pwf, 'or-')
title('Effort from well (flow source) - Norton approximation')
xlabel('time [s]'); ylabel('Pressure [kPa]')
grid on

figure(3)
plot(out.x, dP1, 'ob-', out.x, dP2, 'or-', ...
     out.x, fRep(Qp2, [rho2/1e3 vmu2]).*Qp2, 'ok-')
title('Pressure drop related to piping resistance')
xlabel('time [s]'); ylabel('Pressure [kPa]')
legend('piping section 1', 'piping section 2', 'ESP','Location', 'southeast')
grid on

ExcelTable = [out.x' (Qp1*24*3600)' (Qp2*24*3600)' Ppm' Pwf' dP1' dP2'...
             (fRep(Qp2, [rho2/1e3 vmu2]).*Qp2)'];
writetable(table(ExcelTable),'TEP4240-B-RP-01 01 PLOTS.xlsx','Sheet',1)

    function dydt = system(t, y)
        % ode system
        
        % Fluid properties (piping 1)
        % rho : density             [kg/m3]
        % vmu : dynamic viscosity   [Pa.s]
 
        BO1.P = conversionP(Psm, 'kPa -> psi');
        BO1.T = conversionT(Tw, 'c -> f');
        rho1 = conversionrho(BO1.orho, 'lb/ft3 -> kg/m3');
        vmu1 = BO1.omu*1e-3;
        Ip1 = rho1*L1/Ap1;
        
        % Physical system properties
        % Qp0 : oil production "consumed"
        % Qp1 : oil production to wellhead 
        % Qp2 : oil production to wellhead + diluent
        
        BO1.P = 14.7;
        BO1.T = 60;
        
        Qp1sc = BO1.qsc2ac\[0; y(1)/Ip1; 0];
        Qp1sc = Qp1sc(2)*24*3600;
        
        [Qd, dQd] = Sfd(t);
        Qp0 = Qo - Qp1sc;
        Qp1 = y(1)/Ip1;
        Qp2 = Qp1 + Qd;
        Pwf = fRsw(Qp0)*Qp0;
        
        % Mixing production + diluent
        % Qp1 : production to wellhead 
        % Qd  : diluent rate
        
        BL1.Qo = Qp1;
        BL1.Qd = Qd;
        BL1.BlendingReference
        BO2.TuneViscosity(BL1.bmrmu, BL1.bmrT*1.8 + 32)
        
        % Fluid properties (piping 1)
        % rho : density             [kg/m3]
        % vmu : dynamic viscosity   [Pa.s]
        
        rho2 = conversionrho(BO2.orho, 'lb/ft3 -> kg/m3');
        vmu2 = BO2.omu*1e-3;
        Ip2 = rho2*L2/Ap2;
        Ppm = H(1)*g*rho2/1e3;
        
        % System equations
        % fRpi : resistance from piping
        % fRsw : resistance from well source
        % fRep : resistance from pump
                
        switch indx
            case 1
                % with pump
                num = Pwf ...
                    - fRpi(Qp1, 1, [rho1 vmu1])*Qp1 ...
                    - fRpi(Qp2, 2, [rho2 vmu2])*Qp2 ...
                    - Pwh - (rho1*L1*g + rho2*L2*g)/1e3 ...
                    - dQd*Ip2 + (Ppm - fRep(Qp2, [rho2/1e3 vmu2])*Qp2);
            case 2
                % no pump
                num = Pwf ...
                    - fRpi(Qp1, 1, [rho1 vmu1])*Qp1 ...
                    - fRpi(Qp2, 2, [rho2 vmu2])*Qp2 ...
                    - Pwh - (rho1*g*L1 + rho2*g*L2)/1e3 ...
                    - dQd*Ip2;
        end
        
        den = 1 + Ip2/Ip1;
        dydt = num/den;
    end

    function Rsw = fRsw(Qsc)
        % flow source resistance
        % inflow perfomance index (IPR) - Well Productivity Index (Well PI)

        Rsw = Pws./Qsc - (Pws./(Qsc + (Qsc == 0)) - 1/Jo);
    end

    function Rep = fRep(Qo, props)
        % effort source resistance
        % centrifugal pump performance
        
        n = length(Qo);
        Rep = zeros(1, n);
        
        PPW.SG = props(1);
        PPW.vnu = props(2)*1e3;
        PPW.ViscosityAdjustment
        PPW.FrequencyAdjustment
        
        for i = 1:n
            PPW.Q = conversionF(Qo, 'm3/s -> m3/h');
            Rep(1, i) = (H(1) - PPW.H)*g*PPW.SG/(Qo + (Qo == 0)); 
        end
    end

    function Rpi = fRpi(Qp, indx, props)
        % piping resistance
        % hydraulic losses in piping system
        
        % variables
        % Ap : piping crosssection area
        % di : piping internal diameter
        % Re : Reynolds number
        % er : piping roughness
        % L  : piping length
        
        rho = props(1);
        vmu = props(2);
        
        switch indx
            case 1
                Ap = Ap1; 
                di = di1;
                Re = rho*Qp*di1/(vmu*Ap1);
                er = pr1;
                L = L1;
            case 2
                Ap = Ap2; 
                di = di2;
                Re = rho*Qp*di2/(vmu*Ap2);
                er = pr2;
                L = L2;
        end
        fd = colebrook(Re, er/di);
        Rpi = fd.*(rho/1e3).*Qp*L/(2*di*Ap^2);
    end
    
    function [Qd, dQd] = Sfd(~)
        % flow source
        % diluent injection
        
        % developer's note: 
        % The rate could change in time, with a particular flowrate change
        % for simulating startup.
        
        Qd = 300/(24*3600);
        dQd = 0;
    end

    function ff = colebrook(Re, edr)
        % auxiliary function
        % fiction factor in piping - Colebrook and White (1931) correlation
        
        if Re == 0
            ff = 1;
        elseif Re <= 2000
            ff = 64./Re;
        else
            fff = @(f) (2/log(10))*log(edr/3.7 + 2.51/(Re*sqrt(f))) + 1/sqrt(f);
            dff = @(f) -((2/log(10)*(2.51*0.5/Re)*(edr/3.7 + 2.51/(Re*sqrt(f)))^(-1)*f^(-1.5)) + 0.5*f^(-1.5));
            
            eff = 1;
            ff0 = 1e-3;
            
            while eff > 1e-8
                ff = ff0 - fff(ff0)/dff(ff0);
                eff = abs(ff - ff0);
                ff0 = ff;
            end
        end
    end

    function nT = conversionT(oT, type)
        % auxiliary function
        % conversionT. Temperature
        switch lower(type)
            case 'f -> c'
                nT = (oT - 32)/1.8;
            case 'c -> f'
                nT = oT*1.8 + 32;
        end
    end

    function nP = conversionP(oP, type)
        % auxiliary function
        % conversationP. Pressure
        switch lower(type)
            case 'kpa -> psi'
                nP = oP*14.7/101.325;
            case 'psi -> kpa'
                nP = oP*101.325/14.7;
        end
    end

    function nrho = conversionrho(orho, type)
        % auxiliary function
        % conversionrho. Density
        switch lower(type)
            case 'lb/ft3 -> kg/m3'
                nrho = orho*16.0186;
            case 'kg/m3 -> lb/ft3'
                nrho = orho/16.0186;
        end
    end

    function nF = conversionF(oF, type)
        % auxiliary function
        % conversionF. Volumetric flowrate
        switch lower(type)
            case 'bpd -> m3/h'
                nF = oF*0.1589873/24;
            case 'm3/h -> bpd'
                nF = oF*24/0.1589873;
            case 'm3/h -> m3/s'
                nF = oF/3600;
            case 'm3/s -> m3/h'
                nF = oF*3600;
        end
    end

end
