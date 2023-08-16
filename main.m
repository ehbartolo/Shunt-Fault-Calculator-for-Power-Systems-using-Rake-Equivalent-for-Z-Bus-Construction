function Out = MainFunction(In)

%%%%%%%%%%%%%%%%%% Assumptions : %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% > bus 1 is connected to a generator
% > 'From Bus' of Lines & Two winding  is always less than their 'To Bus'
% > Generator Impedance Set to Xd Initially
% > for the face of the clock (1*30 = 30 deg is how much low side lags high side
% > for connection of transformer: 1-Wye to ground, 2-Delta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%%%%%%%%%%%%% Hard Code Data for Testing Purposes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BaseMVA = 100; 
BusNo = 9;                                                                                     
BasekV = 138; 
FaultBus = 4;
FaultType = 4; % 1-Three phase Fault, 2 - SLGF, 3 - Line to Line, 4- Double line to ground
Zf = 0;  

% ThreeTransformers = [8 3 6 138 69 13.8 150 45 45 0.148i 0.21i 0.369i 1 1 2 0 1; 3 2 1 67 13.8 2.4 15 6 6 0.0665i 0.0469i 0.0151i 1 1 2 0 1];
% TwoTransformers = [4 5 67 13.8 20 0.08i 2 1 1; 6 7 13.8 2.4 6.5 0.0585i 1 2 1];
% Utilities = [9 138 875.3 617.2 100000 100000];
% Generators = [1 2.4 5 0.1 0.2 0.32 0.1];
% Lines = [8 9 138 8.26i 8.26i 31.46i; 3 4 69 2.6i 2.6i 7.8i]
% Motors = [5 13.8 10 0.1 0.2 0.32 0.1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get input from user %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
BaseMVA = input('What is the base MVA?: ');
BusNo = input('What is the Bus Number for Base kV? ');
BasekV = input(strcat('What is the Base Line-Line kV for bus no.', int2str(BusNo), ' ? : ') );
FaultBus = input('What is the faulted Bus Number? : ');
FaultType = input('Type of fault(1:3P, 2:LG, 3:LL, 4:DLG)? : ');
Zf = input('What is the fault impedance(Ohms)?: ');

FaultPhase = 1;
if FaultType == 2
    FaultPhase = input('What is the FAULTED Phase(1:PhaseA, 2:PhaseB, 3:PhaseC)?: ');
elseif FaultType == 3 || FaultType == 4
    FaultPhase = input('What is the UNFAULTED Phase(1:PhaseA, 2:PhaseB, 3:PhaseC)?: ');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Read Data from Excel File %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ThreeTransformers = GetData('FaultStudyData.xls','3transformer');
TwoTransformers = GetData('FaultStudyData.xls','2transformer');
Utilities = GetData('FaultStudyData.xls','utility');
Lines = GetData('FaultStudyData.xls','line');
Generators = GetData('FaultStudyData.xls','gen');
Motors = GetData('FaultStudyData.xls','motor');


% Get Number of Buses in the power system
Maxbus = GetMaxBus(TwoTransformers ,ThreeTransformers,Lines); 
OrigBus = Maxbus;       % used to preserve initial number of buses of bus voltage calculation

% Get Base Values:
BaseValues = GetBaseValues(TwoTransformers,ThreeTransformers,Lines,Maxbus,BusNo,BasekV,BaseMVA);  % Result: Bus Number, Base kV, Base kA, Base Impedance, TAG

% Get p.u. of fault impedance:
Zf = Zf/BaseValues(FaultBus,4);

% Get Impedance Table
[PosZpu,NegZpu,ZerZpu,Maxbus,phaseshift] = GetImpedanceTable(BaseValues,ThreeTransformers,TwoTransformers,Utilities,Generators,Lines,Motors,Maxbus);

% Get Phase shift caused by transformer configuration:
phase = CalculatePhaseShift(phaseshift,Maxbus,FaultBus);

% Get Z bus 
PosZpu = Arrange2(PosZpu,Maxbus);   % Arrange First to avoid a line/transformer that connects to 2 new buses in Zbus Algorithm
NegZpu = Arrange2(NegZpu,Maxbus);
ZerZpu = Arrange2(ZerZpu,Maxbus);
[ZbusposR,ZbusposX] = GetZbus(PosZpu,Maxbus);
[ZbusnegR,ZbusnegX] = GetZbus(NegZpu,Maxbus);
[ZbuszerR,ZbuszerX] = GetZbus(ZerZpu,Maxbus);

% Calculate Fault Bus Voltage in p.u. and line Current in p.u. & kA
[Ifpu,IfkA,VbuspuA,VbuspuB,VbuspuC,IlinepuA,IlinepuB,IlinepuC,IlineAkA,IlineBkA,IlineCkA]  = CalculateFault(ZbusposR,ZbusposX,ZbusnegR,ZbusnegX,ZbuszerR,ZbuszerX,FaultBus,FaultType,BaseValues,Zf,PosZpu,NegZpu,ZerZpu,Maxbus,phase,FaultPhase);

% Display Results:
% I. Base Values:
fprintf('\n\n      ### Table of Base Values ###\n');
fprintf('Bus Number   Base kV   Base kA   Base Ohm \n');
fprintf('----------   -------  ---------  --------\n');
disp(BaseValues);
fprintf('\n\n');

% II. Impedance Table:
fprintf('### Table of POSITIVE Impedance in per unit ###\n');
fprintf(' From Bus   To Bus   Resistance   Reactance \n');
fprintf(' --------   ------   ----------   --------- \n');
disp(Arrange(PosZpu(:,1:4)));
% fprintf('\n\n');

fprintf('### Table of NEGATIVE Impedance in per unit ###\n');
fprintf(' From Bus   To Bus   Resistance   Reactance \n');
fprintf(' --------   ------   ----------   --------- \n');
disp(Arrange(NegZpu(:,1:4)));
% fprintf('\n\n');

fprintf('### Table of ZERO Impedance in per unit ###\n');
fprintf(' From Bus   To Bus   Resistance   Reactance \n');
fprintf(' --------   ------   ----------   --------- \n');
disp(Arrange(ZerZpu(:,1:4)));
fprintf('\n\n');

% III. Zbus Matrix
fprintf('                        ### POSITIVE Zbus Matrix(Magnitude)###\n');
disp(abs(ZbusposX*1i+ZbusposR));
fprintf('\n');
fprintf('                        ### NEGATIVE Zbus Matrix(Magnitude) ###\n');
disp(abs(ZbusnegX*1i+ZbusnegR))
fprintf('\n');
fprintf('                        ### ZERO Zbus Matrix(Magnitude) ###\n');
disp(abs(ZbuszerX*1i+ZbuszerR))
fprintf('\n\n');

% IV. Fault Curent at Faulted Bus
fprintf('Fault Current in Bus %d: %.2e p.u. or %.2e kA \n\n',FaultBus,abs(Ifpu),abs(IfkA))

% V: Bus Voltages in per unit
VbuspuA = VbuspuA(1:OrigBus,:);
VbuspuA = [BaseValues(:,1),abs(VbuspuA),abs(VbuspuA.*BaseValues(:,2)),angle(VbuspuA)*180/pi];
fprintf('### Table of Bus Voltages for Phase A ###\n');
fprintf('Bus Number      P.U       kV    Angle_deg \n');
disp(VbuspuA);
fprintf('\n');

VbuspuB = VbuspuB(1:OrigBus,:);
VbuspuB = [BaseValues(:,1),abs(VbuspuB),abs(VbuspuB.*BaseValues(:,2)),angle(VbuspuB)*180/pi];
fprintf('### Table of Bus Voltages for Phase B ###\n');
fprintf('Bus Number      P.U       kV    Angle_deg \n');
disp(VbuspuB);
fprintf('\n');

VbuspuC = VbuspuC(1:OrigBus,:);
VbuspuC = [BaseValues(:,1),abs(VbuspuC),abs(VbuspuC.*BaseValues(:,2)),angle(VbuspuC)*180/pi];
fprintf('### Table of Bus Voltages for Phase C ###\n');
fprintf('Bus Number      P.U       kV    Angle_deg \n');
disp(VbuspuC);

fprintf('\n\n');

Ia = [real(IlinepuA(:,1:2)), abs(IlinepuA(:,3)),  abs(IlineAkA(:,3)),angle(IlinepuA(:,3))*(180/pi)];
Ib = [real(IlinepuB(:,1:2)), abs(IlinepuB(:,3)),  abs(IlineBkA(:,3)),angle(IlinepuB(:,3))*(180/pi)];
Ic = [real(IlinepuC(:,1:2)), abs(IlinepuC(:,3)),  abs(IlineCkA(:,3)),angle(IlinepuC(:,3))*(180/pi)];

% VI: Line Currents
fprintf('### Table of Line Currents with respect to From_Bus  ###\n');
fprintf('###                     Phase A:                     ###\n');
fprintf('From Bus     To Bus     P.U.      kA      Angle_deg\n');
fprintf('--------     ------     -----     -----   ----------\n');
disp(Ia);
fprintf('\n');
fprintf('###                     Phase B:                     ###\n');
fprintf('From Bus     To Bus     P.U.      kA      Angle_deg\n');
fprintf('--------     ------     -----    ------   ---------\n');
disp(Ib);
fprintf('\n');
fprintf('###                     Phase C:                     ###\n');
fprintf('From Bus     To Bus     P.U.      kA      Angle_deg\n');
fprintf('--------     ------     -----     -----   ---------\n');
disp(Ic);
fprintf('\n');

end

% Converts xlsread data to array of complex numbers
function myVal = GetData(Workbook,Worksheet)

    [num,txt,raw] = xlsread(Workbook,Worksheet);
    raw(1,:) = [];  % delete header of excel
    rawsize = size(raw,2);
    numsize = size(num,2);
    if numsize < rawsize
       num = [num, zeros(size(num,1), rawsize-numsize)];    %Append zeros if num doesn't read all column
    end
    num(isnan(num)) = 0;        % Zero out NAN Values
    raw = str2double(raw);       % Convert complex text to complex number
    raw(isnan(raw)) = 0;        % Zero out NAN Values
    
    myVal = raw + num;
    
    % Note if this gets error then delete succeding rows below the data in
    % excel... this means a spaces are read by xlsread !!!!!
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Get Base Values of Power System %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bVal = GetBaseValues(Data2Trans,Data3Trans,DataLine,NumBus,InitialBus,InitialKV,SBase)
bVal = zeros(NumBus,5);
NxtBus = InitialBus;
RefBus = InitialBus;
RefkV = InitialKV;        % kV of reference Bus for Transformation of Base Values
BasekV = InitialKV;

%Initialize first Bus
bVal(NxtBus,1) = NxtBus;
bVal(NxtBus,2) = BasekV;
bVal(NxtBus,3) = SBase/ (sqrt(3) * BasekV);
bVal(NxtBus,4) = (BasekV)^2/SBase;

Countbus = 1;
while Countbus < NumBus
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%% Get Base Values of Buses connected to Reference Bus %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Rows_to_Delete = [];       % Collects to be deleted Rows of Line and transformers
    % Two winding Transformer
    for j = 1:size(Data2Trans,1)
        for y = 1:2
            if RefBus == Data2Trans(j,y)
                if y == 1
                    y2 = 2;
                elseif y == 2
                    y2 = 1;
                end
                NxtBus = Data2Trans(j,y2);
                if bVal(NxtBus,1) == 0
                    bVal(NxtBus,1) = NxtBus;
                    BasekV = RefkV * (Data2Trans(j,y2+2)/Data2Trans(j,y+2));
                    bVal(NxtBus,2) = BasekV;
                    bVal(NxtBus,3) = SBase/ (sqrt(3) * BasekV);
                    bVal(NxtBus,4) = (BasekV)^2/SBase;
                    Countbus = Countbus +1;
                    Rows_to_Delete = [Rows_to_Delete, j]; % collects to be deleted rows
                end
            end
        end
    end
    Data2Trans(Rows_to_Delete,:) = []; % delete rows of transformer that already have base values
    
    if Countbus == NumBus
        break;					% Exit while loop if Countbus = Maximum Bus
    end
    
    Rows_to_Delete = [];				% collects to be deleted rows
    % Three winding Transformer
    for j = 1:size(Data3Trans,1)
        for y = 1:3
            if RefBus == Data3Trans(j,y)
                for y2 = 1:3
                    if y2 ~= y % not equal
                        NxtBus = Data3Trans(j,y2);
                        if bVal(NxtBus,1) == 0
                            bVal(NxtBus,1) = NxtBus;
                            BasekV = RefkV * (Data3Trans(j,y2 + 3)/Data3Trans(j,y +3)) ;
                            bVal(NxtBus,2) = BasekV;
                            bVal(NxtBus,3) = SBase/ (sqrt(3) * BasekV);
                            bVal(NxtBus,4) = (BasekV)^2/SBase;
                            Countbus = Countbus +1;
                            Rows_to_Delete = [Rows_to_Delete, j]; % collects to be deleted rows
                        end
                    end
                end
            end
        end
    end
    Data3Trans(Rows_to_Delete,:) = [];    % delete rows of transformer that already have base values
    
    if Countbus == NumBus
        break;					% Exit while loop if Countbus = Maximum Bus
    end
    
    Rows_to_Delete = [];				% collects to be deleted rows
    % Lines
    for j = 1:size(DataLine,1)
        for y = 1:2
            if RefBus == DataLine(j,y)
                for y2 = 1:2
                    if y2 ~= y % not equal
                        NxtBus = DataLine(j,y2);
                        if bVal(NxtBus,1) == 0
                            bVal(NxtBus,1) = NxtBus;
                            BasekV = RefkV;
                            bVal(NxtBus,2) = BasekV;
                            bVal(NxtBus,3) = SBase/ (sqrt(3) * BasekV);
                            bVal(NxtBus,4) = (BasekV)^2/SBase;
                            Countbus = Countbus +1;
                            Rows_to_Delete = [Rows_to_Delete, j]; % collects to be deleted rows
                        end
                    end
                end
            end
        end
    end
    DataLine(Rows_to_Delete,:) = [];   % delete rows of transformer that already have base values
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get Next Reference Bus %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    found = 0;				% '1' if found a valid reference bus
    for j = 1:NumBus
        if found == 1
            break;
        end
        if bVal(j,5) == 0
            % Find Reference Bus at the two winding transformer
            for y = 1:size(Data2Trans,1)
                if found == 1
                    break;
                end
                for y2 = 1:2
                    if Data2Trans(y,y2) == bVal(j,1)
                        found = 1;
                        RefBus = bVal(j,1);
                        RefkV = bVal(j,2);
                        bVal(RefBus,5) = 1;  % Tag as Used Reference Bus
                        break;
                    end
                end
            end
            % Find RefBus at the three winding transformer
            for z = 1:size(Data3Trans,1)
                if found == 1
                    break;
                end
                for z2 = 1:3
                    if Data3Trans(z,z2) == bVal(j,1)
                        found = 1;
                        RefBus = bVal(j,1);
                        RefkV = bVal(j,2);
                        bVal(RefBus,5) = 1;  % Tag as Used Reference Bus
                        break;
                    end
                end
            end
            % Find RefBus at the Line
            for z = 1:size(DataLine,1)
                if found == 1
                    break;
                end
                for z2 = 1:2
                    if DataLine(z,z2) == bVal(j,1)
                        found = 1;
                        RefBus = bVal(j,1);
                        RefkV = bVal(j,2);
                        bVal(RefBus,5) = 1;  % Tag as Used Reference Bus
                        break;
                    end
                end
            end
            
        end
    end % Loop for finding RefBus
end % end of while loop

bVal(:,5) = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Get Maximum bus of the power system %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function NumBus = GetMaxBus(Data2Trans,Data3Trans,DataLine)

NumBus = 1;

% Two winding transformer
for j = 1:size(Data2Trans,1)
    for y = 1:2
        if NumBus < Data2Trans(j,y)
            NumBus = Data2Trans(j,y);
        end
    end
end
 
% three winding transformer
for j = 1:size(Data3Trans,1)
     for y = 1:3
         if NumBus < Data3Trans(j,y)
             NumBus = Data3Trans(j,y);
         end
     end
end
 
% lines
for j = 1:size(DataLine,1)
     for y = 1:2
         if NumBus < DataLine(j,y)
             NumBus = DataLine(j,y);
         end
     end
end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Get Impedance Table %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pZ,nZ,zZ,upMaxBus,phShift] = GetImpedanceTable(bVal,thtrans,twtrans,util,gen,line,motor,numbus)
%% Get Z Table: From Bus - To Bus - Resistance pu - Reactance pu - Type(1-Generator,0-branch) - Phase shift%%%%% 
phShift = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Get Per Unit Values of impedance of utility    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Uz1 = util(:,2).*util(:,2)./util(:,3);                                      % positive and negative impedance of utility (X1 = X2 = (kV)^2 / 3PF_MVA)
Uz0 = 3*util(:,2).*util(:,2)./util(:,4) - 2*Uz1; 							% zero impedance X1+X2+X0 = 3*(kV)^2 / SLGF_MVA

Ur1 = Uz1./sqrt(ones(size(util,1),1)+(util(:,5).*util(:,5)));				% R^2 = Z^2/(1+ (X/R)^2)
Ur0 = Uz0./sqrt(ones(size(util,1),1)+(util(:,6).*util(:,6)));

Ux1 = Uz1.*util(:,5)./sqrt(ones(size(util,1),1)+(util(:,5).*util(:,5)));	% X^2 = Z^2 * (X/R)^2 /(1+ (X/R)^2)
Ux0 = Uz0.*util(:,6)./sqrt(ones(size(util,1),1)+(util(:,6).*util(:,6)));


PosUtil = [zeros(size(util,1),1),util(:,1),  Ur1, Ux1,zeros(size(util,1),1)];
NegUtil = PosUtil;
ZerUtil = [zeros(size(util,1),1), util(:,1),  Ur0, Ux0,zeros(size(util,1),1)];

% Get per unit Values
for j = 1:size(PosUtil,1)
	for y = 1:size(bVal,1)
		if PosUtil(j,2) == bVal(y,1)
			PosUtil(j,3:4) = (1/bVal(y,4))*PosUtil(j,3:4);
            NegUtil(j,3:4) = (1/bVal(y,4))*NegUtil(j,3:4);
            ZerUtil(j,3:4) = (1/bVal(y,4))*ZerUtil(j,3:4);
			break;
		end
	end
end

phShift = [phShift;PosUtil(:,1:2),zeros(size(PosUtil,1),1)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     Get Per Unit Values of Impedance of Line   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PosLine = [sort(line(:,1:2),2), real(line(:,4)), imag(line(:,4)), zeros(size(line,1),1)]; % Sort so that FromBus < To Bus
NegLine = [sort(line(:,1:2),2), real(line(:,5)), imag(line(:,5)), zeros(size(line,1),1)];
ZerLine = [sort(line(:,1:2),2), real(line(:,6)), imag(line(:,6)), zeros(size(line,1),1)];

% Get per unit Values
for j = 1:size(PosLine,1)
	for y = 1:size(bVal,1)
		if PosLine(j,1) == bVal(y,1)
			PosLine(j,3:4) = (1/bVal(y,4))*PosLine(j,3:4);
            NegLine(j,3:4) = (1/bVal(y,4))*NegLine(j,3:4);
            ZerLine(j,3:4) = (1/bVal(y,4))*ZerLine(j,3:4);
			break;
		end
	end
end

phShift = [phShift;PosLine(:,1:2),zeros(size(PosLine,1),1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%     Get Per Unit Values of Impedance of Generator   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generator Reactance Set to Steady State Value Initially
PosGen = [zeros(size(gen,1),1),gen(:,1), zeros(size(gen,1),1), gen(:,6).*gen(:,2).*gen(:,2)./gen(:,3),ones(size(gen,1),1)];   % Xpu*(kV^2)/MVA ohms
NegGen = PosGen;
ZerGen = [zeros(size(gen,1),1),gen(:,1), zeros(size(gen,1),1), gen(:,7).*gen(:,2).*gen(:,2)./gen(:,3),ones(size(gen,1),1)];

% Get per unit Values
for j = 1:size(PosGen,1)
	for y = 1:size(bVal,1)
		if PosGen(j,2) == bVal(y,1)
			PosGen(j,3:4) = (1/bVal(y,4)).*PosGen(j,3:4);
            NegGen(j,3:4) = (1/bVal(y,4)).*NegGen(j,3:4);
            ZerGen(j,3:4) = (1/bVal(y,4)).*ZerGen(j,3:4);
			break;
		end
	end
end

phShift = [phShift;PosGen(:,1:2),zeros(size(PosGen,1),1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%     Get Per Unit Values of Impedance of Motor   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Motor Reactance Set to Steady State Value Initially

if size(motor,1) > 0
    PosMotor = [zeros(size(motor,1),1),motor(:,1), zeros(size(motor,1),1), motor(:,6).*motor(:,2).*motor(:,2)./motor(:,3),zeros(size(motor,1),1)];   % Xpu*(kV^2)/MVA ohms
    NegMotor = PosMotor;
    ZerMotor = [zeros(size(motor,1),1),motor(:,1), zeros(size(motor,1),1), motor(:,7).*motor(:,2).*motor(:,2)./motor(:,3),zeros(size(motor,1),1)];

    % Get per unit Values
    for j = 1:size(PosMotor,1)
        for y = 1:size(bVal,1)
            if PosMotor(j,2) == bVal(y,1)
                PosMotor(j,3:4) = (1/bVal(y,4)).*PosMotor(j,3:4);
                NegMotor(j,3:4) = (1/bVal(y,4)).*NegMotor(j,3:4);
                ZerMotor(j,3:4) = (1/bVal(y,4)).*ZerMotor(j,3:4);
                break;
            end
        end
    end

    phShift = [phShift;PosMotor(:,1:2),zeros(size(PosMotor,1),1)];
else
    PosMotor = [];
    NegMotor = [];
    ZerMotor = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     Get Per Unit Values of 3 winding transformer    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%     Assumptions: per unit values are referred to High side and Z1=Z2=Z3  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert to ohmic values.... ZHX,ZHY,ZXY = Zpu*(kV_high_side)^2/MVA_old
thtransX = [thtrans(:,1:3), imag(thtrans(:,10)).*thtrans(:,4).*thtrans(:,4)./thtrans(:,7), imag(thtrans(:,11)).*thtrans(:,4).*thtrans(:,4)./thtrans(:,8), imag(thtrans(:,12)).*thtrans(:,4).*thtrans(:,4)./thtrans(:,9) ];
thtransR = [thtrans(:,1:3), real(thtrans(:,10)).*thtrans(:,4).*thtrans(:,4)./thtrans(:,7), real(thtrans(:,11)).*thtrans(:,4).*thtrans(:,4)./thtrans(:,8), real(thtrans(:,12)).*thtrans(:,4).*thtrans(:,4)./thtrans(:,9) ];

% Convert to ZH = 0.5(Zhx+Zhy-Zxy), ZX = 0.5(Zhx-Zhy+Zxy), ZY = 0.5(-Zhx+Zhy+Zxy) 
thtransX = [thtransX(:,1:3), 0.5*(thtransX(:,4) + thtransX(:,5) - thtransX(:,6)), 0.5*(thtransX(:,4) - thtransX(:,5) + thtransX(:,6)), 0.5*(-thtransX(:,4) + thtransX(:,5) + thtransX(:,6))];
thtransR = [thtransR(:,1:3), 0.5*(thtransR(:,4) + thtransR(:,5) - thtransR(:,6)), 0.5*(thtransR(:,4) - thtransR(:,5) + thtransR(:,6)), 0.5*(-thtransR(:,4) + thtransR(:,5) + thtransR(:,6))];

% Convert to per unit
for j = 1:size(thtransX,1)
    for y = 1:size(bVal,1)  
        if thtransX(j,1) == bVal(y,1)        
           thtransX(j,4:6) = (1/bVal(y,4))*thtransX(j,4:6);
           thtransR(j,4:6) = (1/bVal(y,4))*thtransR(j,4:6);
           break;
        end            
    end
end

% Insert new bus to mutual connection of three winding transformer
Pos3Trans = zeros(3*size(thtransX,1),5);
Zer3Trans = zeros(3*size(thtransX,1),5); % Sixt Column will be assigned to connection of transfomer

for j = 1:size(thtransX,1)
	numbus = numbus +1;
	for y = 1:3
		Pos3Trans(y+3*(j-1),1:5) = [thtransX(j,y), numbus, thtransR(j,y+3), thtransX(j,y+3),0];
        Zer3Trans(y+3*(j-1),1:5) = [thtransX(j,y), numbus, thtransR(j,y+3), thtransX(j,y+3),0];
        if thtrans(j,y+12) == 2     % if Delta Connected Ground the from Bus
            Zer3Trans(y+3*(j-1),1) = 0;
        end
        % Determine phase shift
        if y ~= 1
            phShift = [phShift;numbus,thtransX(j,y),thtrans(j,y+14)];
        else
             phShift = [phShift;numbus,thtransX(j,y),0];
        end
	end
end

Neg3Trans = Pos3Trans;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     Get Per Unit Values of 2 winding transformer    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     Assumption: per unit values are referred to High side   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ohmic Values = Zpu*(kV^2)/MVA
Pos2Trans = [sort(twtrans(:,1:2),2), real(twtrans(:,6)).*twtrans(:,3).*twtrans(:,3)./twtrans(:,5), imag(twtrans(:,6)).*twtrans(:,3).*twtrans(:,3)./twtrans(:,5), zeros(size(twtrans,1),1)];
% Get per unit values
for j = 1:size(Pos2Trans,1)
    for y = 1:size(bVal,1)
        if Pos2Trans(j,1) == bVal(y,1)
            Pos2Trans(j,3:4) = (1/bVal(y,4))*Pos2Trans(j,3:4);
            break;
        end
    end
end
Zer2Trans = Pos2Trans;
for j = 1:size(Zer2Trans,1)
    for k = 1:2
        if twtrans(j,k+6) == 2  % Ground the delta connection of transfomer
           Zer2Trans(j,k) = 0;
        end
    end
end
phShift = [phShift;twtrans(:,1:2),twtrans(:,9)];
%Sort From and To Bus of Impedance:
Zer2Trans = [sort(Zer2Trans(:,1:2),2), Zer2Trans(:,3:5)];
Neg2Trans = Pos2Trans;

% pZ = PosLine;
% nZ = NegLine;
% zZ = ZerLine;
pZ = [PosGen;PosUtil;PosMotor;PosLine;Pos3Trans;Pos2Trans];
nZ = [NegGen;NegUtil;NegMotor;NegLine;Neg3Trans;Neg2Trans];
zZ = [ZerGen;ZerUtil;ZerMotor;ZerLine;Zer3Trans;Zer2Trans];

pZ = Arrange(pZ);
nZ = Arrange(nZ);
zZ = Arrange(zZ);

upMaxBus = numbus;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Arrange Imp. Table such that From Bus < To Bus %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a = Arrange(Zt)
[Y,I] = sort(Zt(:,1));  % Y contains the sorted values of column 1
B = Zt(I,:);            % I contains the row number of the element at the original Matrix
int = 1;
fin = 1;

M = [];
for j = 1:size(B,1)
    if B(j,1) == B(int,1)
        fin = j;
    else
        A = B(int:fin,:);
        [Y,I] = sort(A(:,2));
        M = [M; A(I,:)];
        int = j;
        fin = j;
    end
end
A = B(int:fin,:);       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Y,I] = sort(A(:,2));   %%% For Last Iteration %%%%%%%%%%%%
M = [M; A(I,:)];        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = M;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%  Arrangement based on connection: to Avoid a line/transformer that connects to 2 new buses  %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a2 = Arrange2(Zt,NumBus)
    % Group A: Generators
    % Group B: Impedance not connected to ground
    % Group C: Impedance connected to ground
    GroupA  = [];
    GroupB  = [];
    GroupC  = [];
    GroupBp = [];
    
    for j = 1:size(Zt,1)
        if Zt(j,5) == 1
            GroupA = [GroupA;j];
        elseif Zt(j,1) == 0
            GroupC = [GroupC;j];
        else
            GroupB = [GroupB;j];
        end
    end
    a2 = Zt([GroupA;GroupC;GroupB],:);
    oldbus = zeros(1,NumBus);
    for i = 1:size(GroupA,1)
        oldbus(1,Zt(GroupA(i,1),2)) = 1;
    end
    for i = 1:size(GroupC,1)
        oldbus(1,Zt(GroupC(i,1),2)) = 1;
    end
    for Count = 1:size(GroupB,1)
        for i = 1:size(GroupB,1)
            if GroupB(i,1) > 0
                if oldbus(1,Zt(GroupB(i,1),1)) == 1 || oldbus(1,Zt(GroupB(i,1),2)) == 1
                    oldbus(1,Zt(GroupB(i,1),1)) = 1;
                    oldbus(1,Zt(GroupB(i,1),2)) = 1;
                    GroupBp = [GroupBp;GroupB(i,1)];
                    GroupB(i,1) = 0;
                    break;
                end
            end
        end
    end
    GroupA;
    GroupB = GroupBp;
    GroupC;
    a2 = Zt([GroupA;GroupC;GroupB],:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Z Bus Matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [r,x] = GetZbus(zTable,nBus)
    
    oldBus = zeros(2,nBus); % value of olBus(1,n) will be '1' if bus n is already included in the Zbus, row assignment
    rBus = [zTable(1,3)];
    xBus = [zTable(1,4)];
    
    if zTable(1,1)>0
        oldBus(1,zTable(1,1)) = 1;
        oldBus(2,zTable(1,1)) = 1;
    end
    if zTable(1,2)>0
        oldBus(1,zTable(1,2)) = 1;
        oldBus(2,zTable(1,2)) = 1;
    end
    
    for m = 2:size(zTable,1)
       j = zTable(m,1); % From Bus
       k = zTable(m,2); % To Bus
       % Type 1 Operation (Generator to new bus)
       if zTable(m,1) == 0 & oldBus(1,k) == 0           % (m,5): 1-for generator, 0-for others // oldBus(1,n): 1-oldbus, 0-newbus
           rBus = [rBus,zeros(size(rBus,1),1);zeros(1,size(rBus,2)),zTable(m,3)];
           %fprintf('Type 1\n');
           xBus = [xBus,zeros(size(xBus,1),1);zeros(1,size(xBus,2)),zTable(m,4)];
           oldBus(1,k) = 1;                  % Tag as an old Bus
           oldBus(2,k) = size(rBus,1);       % monitor position of the bus number of the original Impedance Table relative to Zbus Matrix
           %oldBus
       % Type 2 Operation (Generator to old Bus)
       elseif zTable(m,1) == 0 & oldBus(1,k) == 1
          rBus = [rBus, rBus(:,oldBus(2,k)); rBus(oldBus(2,k),:), rBus(oldBus(2,k),oldBus(2,k))+zTable(m,3)];
          xBus = [xBus, xBus(:,oldBus(2,k)); xBus(oldBus(2,k),:), xBus(oldBus(2,k),oldBus(2,k))+zTable(m,4)]; 
          %fprintf('Type 2: oldBus = %d\n',oldBus(2,k));
          rBus = KronReduction(rBus);
          xBus = KronReduction(xBus);
       else 
           if j ~= 0 % Don't include non-generator that is connected to ground.....
                % Type 4 Operation (old to old)
                if (oldBus(1,j) == 1) & (oldBus(1,k) == 1)
                     j = oldBus(2,j);
                     k = oldBus(2,k);
                     rV = rBus(j,j) + rBus(k,k) - 2*rBus(j,k) + zTable(m,3);       % Zv = Zjj + Zkk - 2Zjk + Zb                                                             
                     xV = xBus(j,j) + xBus(k,k) - 2*xBus(j,k) + zTable(m,4);       %
                     rBus = [rBus,rBus(:,j)-rBus(:,k); rBus(j,:)-rBus(k,:), rV];    
                     %fprintf('Type 4: oldBuses: %d,%d\n',j,k);
                     xBus = [xBus,xBus(:,j)-xBus(:,k); xBus(j,:)-xBus(k,:), xV];

                     rBus = KronReduction(rBus);
                     xBus = KronReduction(xBus);
                % Type 3 Operation (old to new)
                else
                    if oldBus(1,j) == 0
                       addBus = j;      %addBus is the new Bus to be included in the old bus
                       oBus = k;        %oBus is an old bus
                       oBus = oldBus(2,k);
                    else
                       addBus = k;
                       oBus = j;
                       oBus = oldBus(2,j);
                    end                  
                    rBus = [rBus, rBus(:,oBus); rBus(oBus,:), rBus(oBus,oBus)+zTable(m,3)];
                    %fprintf('Type 3: oldBus: %d\n',oBus);
                    xBus = [xBus, xBus(:,oBus); xBus(oBus,:), xBus(oBus,oBus)+zTable(m,4)];
                    oldBus(1,addBus) = 1;
                    oldBus(2,addBus) = size(rBus,1);
                end
           else
           end
       end
    end
    
    z =  transpose(oldBus(2,:));
    
    r = zeros(nBus,nBus);
    x = zeros(nBus,nBus);
    
    % Re-Arrange Values of Original Matrix
    for j = 1:nBus
        for k = 1:nBus
            if (z(j,1) ~= 0) & (z(k,1) ~= 0)
                r(j,k) = rBus(z(j,1),z(k,1));
            end
        end
    end
    for j = 1:nBus
        for k = 1:nBus
            if (z(j,1) ~= 0) & (z(k,1) ~= 0)
                x(j,k) = xBus(z(j,1),z(k,1));
            end
        end
    end

end

function k = KronReduction(Matrix1)
    
    % Matrix' = X2*X3*(1/X4)
    X1 = Matrix1(1:(size(Matrix1,1)-1),1:(size(Matrix1,2)-1));
    X2 = Matrix1(1:(size(Matrix1,1)-1), size(Matrix1,2));
    X3 = Matrix1(size(Matrix1,1),1:(size(Matrix1,2)-1));
    X4 = Matrix1(size(Matrix1,1), size(Matrix1,2));
    
    if X4 == 0 
        k = X1;
    else
        k = X1-X2*(1/X4)*X3;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% FaultType: 1-Three Phase Fault, 2- Single Line to Ground Fault, 3- Line to Line , 4- Double Line to Ground%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ifp,IfA,VbpA,VbpB,VbpC,IlpA,IlpB,IlpC,IlaA,IlaB,IlaC] = CalculateFault(RbusPos, XbusPos, RbusNeg, XbusNeg, RbusZer, XbusZer, FaultedBusNo, FaultType, BaseMatrix, ZFault, PosImpTab, NegImpTab,ZerImpTab,MaximumBus,FOC,faultphase)
   
   PosImpTab = Arrange(PosImpTab);
   NegImpTab = Arrange(NegImpTab);
   ZerImpTab = Arrange(ZerImpTab);
   A = exp(sqrt(-1)*120*pi/180);
   A2 = exp(sqrt(-1)*240*pi/180);
   
   n = FaultedBusNo;
   nBus = MaximumBus;
   IlpA = zeros(size(PosImpTab,1),3); % Line current in per unit phase A
   IlpB = zeros(size(PosImpTab,1),3); % Line current in per unit phase B
   IlpC = zeros(size(PosImpTab,1),3); % Line current in per unit phase C
   
   IlpA(:,1:2) = PosImpTab(:,1:2);
   IlpB(:,1:2) = PosImpTab(:,1:2);
   IlpC(:,1:2) = PosImpTab(:,1:2);
   
   Ilp1 = IlpA;                       % Positive Sequence Line current in per unit 
   Ilp2 = IlpA;                       % Negative Sequence Line current in per unit 
   Ilp0 = IlpA;                       % Zero Sequence Line current in per unit 
   
   
   if FaultType  == 1 % Three Phase Fault.... Balance System (Ia0 = Ia2 = 0)
        Ifp = 1/(RbusPos(n,n)+XbusPos(n,n)*1i+ZFault) ; %Ik = 1/Zkk
        IfA = Ifp*BaseMatrix(n,3);
        VbpA = ones(nBus,1) - (RbusPos(1:nBus,n) + XbusPos(1:nBus,n)*1i) .* Ifp; % Vn = 1-Znk/Zkk
  
       
        for j = 1:size(PosImpTab,1)
            if PosImpTab(j,1) == 0 
                Vm  = 1;
            else
                Vm = VbpA(PosImpTab(j,1),1);
            end
            Vn = VbpA(PosImpTab(j,2),1);
            % Get phase shift:
            % Maximum phase shift of the two buses. Phase Shift = e^(j*FOC*30*pi/180) = e^(j*FOC*pi/6)  %%%%%%%%%%
            if PosImpTab(j,1) == 0 
                ang1  = 0;
            else
                ang1 = FOC(PosImpTab(j,1),1);
            end
            if PosImpTab(j,2) == 0 
                ang2  = 0;
            else
                ang2 = FOC(PosImpTab(j,2),1);
            end
            nFOC = max(ang1,ang2);
             
            IlpA(j,3) = ( (Vm-Vn)/complex(PosImpTab(j,3),PosImpTab(j,4)) ) * exp(sqrt(-1)*nFOC*(-1)*pi/6); % Imn = (Vm-Vn)/zmn * (phase shift <0) 
        end
        VbpA = VbpA .* exp(sqrt(-1)*FOC*(-1)*pi/6); %**********************************
        VbpB = VbpA*A2;                             % Apply phase shift to bus voltages
        VbpC = VbpA*A;                              %**********************************
        
        IlpB(:,3) = IlpA(:,3)*A2; %Ilp = IlpB<240 deg)
        IlpC(:,3) = IlpA(:,3)*A; %Ilp = IlpB<120 deg)
   elseif FaultType == 2 % Single Line to Ground Fault (Ia0 = Ia1 = Ia2)
       Ia0 = 1/(RbusPos(n,n)+XbusPos(n,n)*1i + RbusNeg(n,n)+XbusNeg(n,n)*1i+ RbusZer(n,n)+XbusZer(n,n)*1i + 3*(ZFault)); % Ia0 = Ia1 = Ia2 = 1/(Z1+Z2+Z0+3Zf)

           
       if faultphase == 2                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           Ia1 = Ia0/A2;                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           Ia2 = Ia0/A;                         %%%%%%
       elseif faultphase == 3                   %%%%%% For Change in Symmetry
           Ia1 = Ia0/A;                         %%%%%%
           Ia2 = Ia0/A2;                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       else                                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           Ia1 = Ia0;                           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           Ia2 = Ia0;                           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       end                                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       Ifp = 3*Ia0;                 % Ia = Ia1 + Ia2 + Ia0
       IfA = Ifp*BaseMatrix(n,3);
       
       Vbp1 = ones(nBus,1) - (RbusPos(1:nBus,n) + XbusPos(1:nBus,n)*1i) .* Ia1;% Va1 = 1-Z1*Ia1
       Vbp2 = -(RbusNeg(1:nBus,n) + XbusNeg(1:nBus,n)*1i) .* Ia2;              % Va2 = -Z2*Ia2
       Vbp0 = -(RbusZer(1:nBus,n) + XbusZer(1:nBus,n)*1i) .* Ia0;              % Va0 = -Z0*Ia0
       

       k = 1;
       % Compute for line current
       for j = 1:size(PosImpTab,1)
            % Solve for positive sequence line current
            if PosImpTab(j,1) == 0 
                Vm  = 1;
            else
                Vm = Vbp1(PosImpTab(j,1),1);
            end
            if PosImpTab(j,2) == 0 
                Vn  = 1;
            else
                Vn = Vbp1(PosImpTab(j,2),1);
            end
            % Get phase shift:
            % Maximum phase shift of the two buses. Phase Shift = e^(j*FOC*30*pi/180) = e^(j*FOC*pi/6)  %%%%%%%%%%
            if PosImpTab(j,1) == 0 
                ang1  = 0;
            else
                ang1 = FOC(PosImpTab(j,1),1);
            end
            if PosImpTab(j,2) == 0 
                ang2  = 0;
            else
                ang2 = FOC(PosImpTab(j,2),1);
            end
            nFOC = max(ang1,ang2);
            Ilp1(j,3) =( (Vm-Vn)/complex(PosImpTab(j,3),PosImpTab(j,4)) ) * exp(sqrt(-1)*(-1)*nFOC*pi/6); % Imn = ((Vm-Vn)/zmn)(phase shift<0)
            
            % end of positive line current
            % Solve for negative sequence line current
            if NegImpTab(j,1) == 0 
                Vm  = 0;
            else
                Vm = Vbp2(NegImpTab(j,1),1);
            end
            if NegImpTab(j,2) == 0 
                Vn  = 0;
            else
                Vn = Vbp2(NegImpTab(j,2),1);
            end
            
            % Get phase shift of negative line current:
            % Maximum phase shift of the two buses. Phase Shift = e^(j*FOC*30*pi/180) = e^(j*FOC*pi/6)  %%%%%%%%%%
            if NegImpTab(j,1) == 0 
                ang1  = 0;
            else
                ang1 = FOC(NegImpTab(j,1),1);
            end
            if NegImpTab(j,2) == 0 
                ang2  = 0;
            else
                ang2 = FOC(NegImpTab(j,2),1);
            end
            nFOC = max(ang1,ang2);
            Ilp2(j,3) =( (Vm-Vn)/complex(NegImpTab(j,3),NegImpTab(j,4)) ) * exp(sqrt(-1)*(+1)*nFOC*pi/6); % Imn = ((Vm-Vn)/zmn)(phase shift>0)      
            
            zmn0 = 0;
            % finds if there is a corresponding zero sequence impedance & line current...
            while k <=size(ZerImpTab,1)
                if ZerImpTab(k,1) > PosImpTab(j,1)
                   break; 
                elseif (ZerImpTab(k,1) == PosImpTab(j,1)) & (ZerImpTab(k,2) == PosImpTab(j,2))
                    zmn0 = complex(ZerImpTab(k,3),ZerImpTab(k,4));
                    break;
                else
                    k = k +1;
                end
            end
            if zmn0 ~= 0
                % Solve for zero sequence line current
                if PosImpTab(j,1) == 0 
                    Vm  = 0;
                else
                    Vm = Vbp0(PosImpTab(j,1),1);
                end
                if PosImpTab(j,2) == 0 
                    Vn  = 0;
                else
                    Vn = Vbp0(PosImpTab(j,2),1);
                end     
                Ilp0(j,3) = (Vm-Vn)/zmn0; % Imn = (Vm-Vn)/zmn
            else
                Ilp0(j,3) = 0;
            end
                
       end
       IlpA(:,3) = Ilp0(:,3) + Ilp1(:,3) + Ilp2(:,3) ;
       IlpB(:,3) = Ilp0(:,3) + A2* Ilp1(:,3) + A*Ilp2(:,3);
       IlpC(:,3) = Ilp0(:,3) + A*Ilp1(:,3) + A2*Ilp2(:,3) ;
       
       
       Vbp1 = Vbp1 .*exp(sqrt(-1)*(-1)*FOC*pi/6);   %
       Vbp2 = Vbp2 .*exp(sqrt(-1)*(+1)*FOC*pi/6);   %
       VbpA = Vbp1 + Vbp2 + Vbp0;                   % Apply phase shift in bus voltages
       VbpB = Vbp1*A2 + Vbp2*A + Vbp0;              %
       VbpC = Vbp1*A + Vbp2*A2 + Vbp0;              %
       
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Line to Line Fault %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   elseif FaultType == 3 % Line to Line Fault ( Ia0 = 0, Ia1 = -Ia2 )
       Ia0 = 0;
       Ia1 = 1/(RbusPos(n,n)+XbusPos(n,n)*1i + RbusNeg(n,n)+XbusNeg(n,n)*1i+ (ZFault)); % Ia0 = Ia1 = Ia2 = 1/(Z1+Z2+Zf)
       Ia2 = (-1) * Ia1;
       
       %%%%%%%%%%%%%%%%%% For Change in Symmetry %%%%%%%%%%%%%%%%%%%%%%%%
       if faultphase == 2
           Ia1 = Ia1/A2;
           Ia2 = Ia2/A;
           Ifp = Ia0 + A*Ia1 + A2*Ia2;                 % Since IfpB = 0 then Ifault = IfC = Ia0+ a*Ia1 + a2*Ia2
       elseif faultphase == 3
           Ia1 = Ia1/A;
           Ia2 = Ia2/A2;
           Ifp = Ia0 + Ia1 + Ia2;                      % Since IfpC = 0 then Ifault = IfA = Ia0+ Ia1 + Ia2
       else
           Ifp = Ia0 + A2*Ia1 + A*Ia2;                 % Since IfpA = 0 then Ifault = IfB = Ia0+ a2*Ia1 + a*Ia2
       end
       %%%%%%%%%%%%%%%%%% For Change in Symmetry %%%%%%%%%%%%%%%%%%%%%%%%
       
       
       IfA = Ifp*BaseMatrix(n,3);
       
       Vbp1 = ones(nBus,1) - (RbusPos(1:nBus,n) + XbusPos(1:nBus,n)*1i) .* Ia1;% Va1 = 1-Zjk1*Ia1
       Vbp2 = -(RbusNeg(1:nBus,n) + XbusNeg(1:nBus,n)*1i) .* Ia2;              % Va2 = -Zjk2*Ia2
       Vbp0 = zeros(nBus,1);                                                   % Va0 = 0
       VbpA = Vbp1 + Vbp2 + Vbp0;
       VbpB = Vbp1*A2 + Vbp2*A + Vbp0;
       VbpC = Vbp1*A + Vbp2*A2 + Vbp0;
       
       k = 1;
       % Compute for line current
       for j = 1:size(PosImpTab,1)
            % Solve for positive sequence line current
            if PosImpTab(j,1) == 0 
                Vm  = 1;
            else
                Vm = Vbp1(PosImpTab(j,1),1);
            end
            if PosImpTab(j,2) == 0 
                Vn  = 1;
            else
                Vn = Vbp1(PosImpTab(j,2),1);
            end
            % Get phase shift:
            % Maximum phase shift of the two buses. Phase Shift = e^(j*FOC*30*pi/180) = e^(j*FOC*pi/6)  %%%%%%%%%%
            if PosImpTab(j,1) == 0 
                ang1  = 0;
            else
                ang1 = FOC(PosImpTab(j,1),1);
            end
            if PosImpTab(j,2) == 0 
                ang2  = 0;
            else
                ang2 = FOC(PosImpTab(j,2),1);
            end
            nFOC = max(ang1,ang2);
            Ilp1(j,3) =( (Vm-Vn)/complex(PosImpTab(j,3),PosImpTab(j,4)) ) * exp(sqrt(-1)*(-1)*nFOC*pi/6); % Imn = ((Vm-Vn)/zmn)(phase shift<0)
            
            % end of positive line current
            % Solve for negative sequence line current
            if NegImpTab(j,1) == 0 
                Vm  = 0;
            else
                Vm = Vbp2(NegImpTab(j,1),1);
            end
            if NegImpTab(j,2) == 0 
                Vn  = 0;
            else
                Vn = Vbp2(NegImpTab(j,2),1);
            end
            
            % Get phase shift of negative line current:
            % Maximum phase shift of the two buses. Phase Shift = e^(j*FOC*30*pi/180) = e^(j*FOC*pi/6)  %%%%%%%%%%
            if NegImpTab(j,1) == 0 
                ang1  = 0;
            else
                ang1 = FOC(NegImpTab(j,1),1);
            end
            if NegImpTab(j,2) == 0 
                ang2  = 0;
            else
                ang2 = FOC(NegImpTab(j,2),1);
            end
            nFOC = max(ang1,ang2);
            Ilp2(j,3) =( (Vm-Vn)/complex(NegImpTab(j,3),NegImpTab(j,4)) ) * exp(sqrt(-1)*(+1)*nFOC*pi/6); % Imn = ((Vm-Vn)/zmn)(phase shift>0)      
            
            zmn0 = 0;
            % finds if there is a corresponding zero sequence impedance & line current...
            while k <=size(ZerImpTab,1)
                if ZerImpTab(k,1) > PosImpTab(j,1)
                   break; 
                elseif (ZerImpTab(k,1) == PosImpTab(j,1)) & (ZerImpTab(k,2) == PosImpTab(j,2))
                    zmn0 = complex(ZerImpTab(k,3),ZerImpTab(k,4));
                    break;
                else
                    k = k +1;
                end
            end
            if zmn0 ~= 0
                % Solve for zero sequence line current
                if PosImpTab(j,1) == 0 
                    Vm  = 0;
                else
                    Vm = Vbp0(PosImpTab(j,1),1);
                end
                if PosImpTab(j,2) == 0 
                    Vn  = 0;
                else
                    Vn = Vbp0(PosImpTab(j,2),1);
                end     
                Ilp0(j,3) = (Vm-Vn)/zmn0; % Imn = (Vm-Vn)/zmn
            else
                Ilp0(j,3) = 0;
            end
                
       end
       IlpA(:,3) = Ilp0(:,3) + Ilp1(:,3) + Ilp2(:,3) ;
       IlpB(:,3) = Ilp0(:,3) + A2* Ilp1(:,3) + A*Ilp2(:,3);
       IlpC(:,3) = Ilp0(:,3) + A*Ilp1(:,3) + A2*Ilp2(:,3) ;
       
       Vbp1 = Vbp1 .*exp(sqrt(-1)*(-1)*FOC*pi/6);   %
       Vbp2 = Vbp2 .*exp(sqrt(-1)*(+1)*FOC*pi/6);   %
       VbpA = Vbp1 + Vbp2 + Vbp0;                   % Apply phase shift in bus voltages
       VbpB = Vbp1*A2 + Vbp2*A + Vbp0;              %
       VbpC = Vbp1*A + Vbp2*A2 + Vbp0;              %
       
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Double Line to Ground Fault %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   elseif FaultType == 4 % Double Line to Ground Fault
       Z0T = RbusZer(n,n)+XbusZer(n,n)*1i + ZFault; % Z0T = Z0 + Zf + 3 (ZG = 0)
       Z1T = RbusPos(n,n)+XbusPos(n,n)*1i + ZFault; % Z1T = Z1 + Zf
       Z2T = RbusNeg(n,n)+XbusNeg(n,n)*1i + ZFault; % Z1T = Z1 + Zf
       
       Ia1 = 1/(Z1T + ((Z0T*Z2T)/(Z0T+Z2T)));      %Ia1 = 1/(Z1T + (Z0T//Z2T) ) 
       Ia0 = (-1)*(Z2T/ (Z0T+Z2T)) * Ia1;           %Ia0 = -Z2T/(Z0T+Z2T)
       Ia2 = (-1)*(Z0T/ (Z0T+Z2T)) * Ia1;           %Ia2 = -Z0T/(Z0T+Z2T)
       
      %%%%%%%%%%%%%%%%%% For Change in Symmetry %%%%%%%%%%%%%%%%%%%%%%%%
       if faultphase == 2
           Ia1 = Ia1/A2;
           Ia2 = Ia2/A;
           Ifp = Ia0 + A*Ia1 + A2*Ia2;                 % Since IfpB = 0 then Ifault = IfC = Ia0+ a*Ia1 + a2*Ia2
       elseif faultphase == 3
           Ia1 = Ia1/A;
           Ia2 = Ia2/A2;
           Ifp = Ia0 + Ia1 + Ia2;                      % Since IfpC = 0 then Ifault = IfA = Ia0+ Ia1 + Ia2
       else
           Ifp = Ia0 + A2*Ia1 + A*Ia2;                 % Since IfpA = 0 then Ifault = IfB = Ia0+ a2*Ia1 + a*Ia2
       end
       %%%%%%%%%%%%%%%%%% For Change in Symmetry %%%%%%%%%%%%%%%%%%%%%%%%       
       
       
       IfA = Ifp*BaseMatrix(n,3);
       
       Vbp1 = ones(nBus,1) - (RbusPos(1:nBus,n) + XbusPos(1:nBus,n)*1i) .* Ia1;% Va1 = 1-Zjk1*Ia1
       Vbp2 = -(RbusNeg(1:nBus,n) + XbusNeg(1:nBus,n)*1i) .* Ia2;              % Va2 = -Zjk2*Ia2
       Vbp0 = -(RbusZer(1:nBus,n) + XbusZer(1:nBus,n)*1i) .* Ia0;              % Va0 = -Zjk0*Ia0
       
       k = 1;
       % Compute for line current
       for j = 1:size(PosImpTab,1)
            % Solve for positive sequence line current
            if PosImpTab(j,1) == 0 
                Vm  = 1;
            else
                Vm = Vbp1(PosImpTab(j,1),1);
            end
            if PosImpTab(j,2) == 0 
                Vn  = 1;
            else
                Vn = Vbp1(PosImpTab(j,2),1);
            end
            % Get phase shift:
            % Maximum phase shift of the two buses. Phase Shift = e^(j*FOC*30*pi/180) = e^(j*FOC*pi/6)  %%%%%%%%%%
            if PosImpTab(j,1) == 0 
                ang1  = 0;
            else
                ang1 = FOC(PosImpTab(j,1),1);
            end
            if PosImpTab(j,2) == 0 
                ang2  = 0;
            else
                ang2 = FOC(PosImpTab(j,2),1);
            end
            nFOC = max(ang1,ang2);
            Ilp1(j,3) =( (Vm-Vn)/complex(PosImpTab(j,3),PosImpTab(j,4)) ) * exp(sqrt(-1)*(-1)*nFOC*pi/6); % Imn = ((Vm-Vn)/zmn)(phase shift<0)
            
            % end of positive line current
            % Solve for negative sequence line current
            if NegImpTab(j,1) == 0 
                Vm  = 0;
            else
                Vm = Vbp2(NegImpTab(j,1),1);
            end
            if NegImpTab(j,2) == 0 
                Vn  = 0;
            else
                Vn = Vbp2(NegImpTab(j,2),1);
            end
            
            % Get phase shift of negative line current:
            % Maximum phase shift of the two buses. Phase Shift = e^(j*FOC*30*pi/180) = e^(j*FOC*pi/6)  %%%%%%%%%%
            if NegImpTab(j,1) == 0 
                ang1  = 0;
            else
                ang1 = FOC(NegImpTab(j,1),1);
            end
            if NegImpTab(j,2) == 0 
                ang2  = 0;
            else
                ang2 = FOC(NegImpTab(j,2),1);
            end
            nFOC = max(ang1,ang2);
            Ilp2(j,3) =( (Vm-Vn)/complex(NegImpTab(j,3),NegImpTab(j,4)) ) * exp(sqrt(-1)*(+1)*nFOC*pi/6); % Imn = ((Vm-Vn)/zmn)(phase shift>0)      
            
            zmn0 = 0;
            % finds if there is a corresponding zero sequence impedance & line current...
            while k <=size(ZerImpTab,1)
                if ZerImpTab(k,1) > PosImpTab(j,1)
                   break; 
                elseif (ZerImpTab(k,1) == PosImpTab(j,1)) & (ZerImpTab(k,2) == PosImpTab(j,2))
                    zmn0 = complex(ZerImpTab(k,3),ZerImpTab(k,4));
                    break;
                else
                    k = k +1;
                end
            end
            if zmn0 ~= 0
                % Solve for zero sequence line current
                if PosImpTab(j,1) == 0 
                    Vm  = 0;
                else
                    Vm = Vbp0(PosImpTab(j,1),1);
                end
                if PosImpTab(j,2) == 0 
                    Vn  = 0;
                else
                    Vn = Vbp0(PosImpTab(j,2),1);
                end     
                Ilp0(j,3) = (Vm-Vn)/zmn0; % Imn = (Vm-Vn)/zmn
            else
                Ilp0(j,3) = 0;
            end
                
       end
       IlpA(:,3) = Ilp0(:,3) + Ilp1(:,3) + Ilp2(:,3) ;
       IlpB(:,3) = Ilp0(:,3) + A2* Ilp1(:,3) + A*Ilp2(:,3);
       IlpC(:,3) = Ilp0(:,3) + A*Ilp1(:,3) + A2*Ilp2(:,3) ;
       
       Vbp1 = Vbp1 .*exp(sqrt(-1)*(-1)*FOC*pi/6);   %
       Vbp2 = Vbp2 .*exp(sqrt(-1)*(+1)*FOC*pi/6);   %
       VbpA = Vbp1 + Vbp2 + Vbp0;                   % Apply phase shift in bus voltages
       VbpB = Vbp1*A2 + Vbp2*A + Vbp0;              %
       VbpC = Vbp1*A + Vbp2*A2 + Vbp0;              %
   end
   
   IlaA = IlpA;
   IlaB = IlpB;
   IlaC = IlpC;
   
   for i = 1:size(IlaA,1)
       if IlaA(i,1) == 0
           IlaA(i,3) = IlaA(i,3)*BaseMatrix(IlaA(i,2),3);
           IlaB(i,3) = IlaB(i,3)*BaseMatrix(IlaB(i,2),3);
           IlaC(i,3) = IlaC(i,3)*BaseMatrix(IlaC(i,2),3);
       else
           IlaA(i,3) = IlaA(i,3)*BaseMatrix(IlaA(i,1),3);
           IlaB(i,3) = IlaB(i,3)*BaseMatrix(IlaB(i,1),3);
           IlaC(i,3) = IlaC(i,3)*BaseMatrix(IlaC(i,1),3);
       end
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Calculates the phase shift(lag) for each bus  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function phs = CalculatePhaseShift(phashift,nBus,fBus)

intbus = fBus;
y = 1;
phs = zeros(nBus,3);
phs(intbus,1) =1;
phs(intbus,2) =1;

% phs: (0 or 1) used as initial bus, (0 or 1) has computed phase shift,FOC
while y < nBus
    rowdel = []; % delete rows to avoid redundancy
    for i = 1: size(phashift,1)
        % column 1 of phashift is the high side and column 2 of phashift is the low %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
        if (phashift(i,1) == intbus) && phashift(i,2) > 0
            if(phs(phashift(i,2),2) ~= 1)
                phs(phashift(i,2),3) = phs(intbus,3) + phashift(i,3);
                rowdel = [rowdel;i];
                phs(phashift(i,2),2) = 1;
                y = y +1;
            end
        elseif (phashift(i,2) == intbus) && phashift(i,1) > 0
            if(phs(phashift(i,1),2) ~= 1)
                phs(phashift(i,1),3) = phs(intbus,3) - phashift(i,3); %verified
                rowdel = [rowdel;i];
                phs(phashift(i,1),2) = 1;
                y = y+1;
            end
        end

    end
    
    phashift(rowdel,:) = [];
    
    found = 0; % 1 if found another initial bus
    for z = 1:size(phashift,1)
        for z2 = 1:2
            if phashift(z,z2) > 0
                if phs(phashift(z,z2),1) == 0 && phs(phashift(z,z2),2) == 1
                    found = 1;
                    intbus = phashift(z,z2);
                    phs(phashift(z,z2),1) = 1;
                    break;
                end
            end
        end
        if found == 1
            break;
        end
    end
   
end
phs = phs(:,3);
end
