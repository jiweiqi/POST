%%
%{

This program is for postprocess pressure trace data from Tsinghua RCM,
The original pressure data filename should be like:20140113_Test_18.txt,

Updating at

https://github.com/jiweiqi/rcm-pressure-post

%}

clear;clc;close all;

%% INPUT INITIAL PARAMETERS
processible = 0;            % If the current run could not be processed for some reason,please change this flag to 0, and highlight the reason in reamrk. 
dataFilename = '20140817_Test_20.txt';
rcm.initial.P = 0.7699;     % bar
rcm.initial.T = 298.25;     % K
rcm.bool_sampling = 1;      % 0 is corresponding to without sampling  
rcm.sampling.delay = 14;	% ms
rcm.config = '10';          % mm
Remark = '';                % Remark about the present experiment.
rcm.sampling.length = 2;    % Scale in [ms].default is 2ms
rcm.ign_stage = 1;          % One-stage ignition or Two-stage ignition.
rcm.TDC_lower_bool =1;      % Set 0 when no visuable drop of pressure befor ignition
rcm.phi = '0.4';            

%% look up TDC(TOP DEAD CENTER) & ignition point
rcm.lookup.TDC.begin_lower = 0.480 *10^5;       
                                            %Input LOWER limit for looking up TDC 
rcm.lookup.TDC.begin_upper = 0.4901*10^5;   
                                            %Input UPPER limit for looking up TDC 
rcm.lookup.TDC.eff_begin_lower = rcm.lookup.TDC.begin_lower;   
                                            %Input LOWER limit for looking up beginning of Effective pressure
rcm.lookup.TDC.eff_begin_upper = rcm.lookup.TDC.begin_upper;  
                                            %Input UPPER limit for looking up End of Effective
rcm.lookup.TDC.finish_upper = 0.5058*10^5;      
                                            %Input UPPER limit for looking up End of Effective pressure
                                            %Also used as LOWER limit for looking up total ignition point
                                            
rcm.lookup.ign.total_lower = rcm.lookup.TDC.finish_upper;    
rcm.lookup.ign_upper = 10^5-100;      
                                            %Default 100, Specify it when multi peaks for dp/dt exist
   
                                            %Set this start point [only if first stage happens]
                                            %Input LOWER limit for looking up Total ign
rcm.lookup.ign.first_lower = 0.489 *10^5;   %Input LOWER limit for looking up Fisrt ign [only if first stage happens]
rcm.lookup.ign.first_upper = 0.4897*10^5;   %Input UPPER limit for looking up Fisrt ign [only if first stage happens]

%% Fuel Components
rcm.fuel = 'Isobutanol';
Molar_SUM = 0.4 + 6*(1+7.26);
rcm.Molar = struct(...
    'O2',6/Molar_SUM , ...
    'N2',6*7.26*1/2.8/Molar_SUM , ...
    'Ar',6*7.26*1.8/2.8/Molar_SUM , ...
    'H2',0, ...
    'CO2',0, ...
    'CO',0,...
    'H2O',0, ...
    'nHeptane',0, ...
    'Toluene',0, ...
    'Isooctane',0,...
    'MCH',0, ...
    'CH4',0,...
    'Isobutanol',0.4/Molar_SUM ,...
    'DMH26',0);

%% Coefficients for Pressure Trancedure and Filting
rcm.P.AMP = 20;     %bar/V
rcm.P.Filter.hrrCOF = 3000; %Hz
rcm.P.Filter.NyquistFreq = 10^5/2;

rcm.sampling.P_AMP = 0.5;	%bar/V
rcm.sampling.Filter.hrrCOF = 3000;
rcm.sampling.Filter.NyquistFreq = 10^5/2;

%% Primary Settings and Checkings
FontName = 'Times New Roman';
FontSize = 15;
FontWeight = 'Bold';
LineWidth = 1.5;
no_data_sign = '/';

rcm.Molar.Dilution_ratio = (rcm.Molar.N2 + rcm.Molar.Ar + rcm.Molar.CO2 )./rcm.Molar.O2;

Mole_sum = rcm.Molar.O2 + rcm.Molar.N2 + rcm.Molar.Ar + rcm.Molar.H2 + rcm.Molar.CO2 + ...
    rcm.Molar.CO + rcm.Molar.H2O + ...
        rcm.Molar.nHeptane + rcm.Molar.Toluene + rcm.Molar.Isooctane ...
        + rcm.Molar.MCH + rcm.Molar.Isobutanol +rcm.Molar.CH4 +rcm.Molar.DMH26; 

assert( abs(Mole_sum-1)< 0.001 , ...
        ['The Sum of Molar Fraction is Not Unity!' 10 ...
        'Please Check Molar Fraction!'] );
    
assert( rcm.bool_sampling ==0 | rcm.bool_sampling == 1,  ...
    ['Please input 0 or 1 for sampling switch !' 10 ...
    'Input 0 for without Sampling.' 10 ...
    'Input 1 if there is Sampling.']);

assert( rcm.ign_stage == 1 | rcm.ign_stage == 2, ...
    ['Please input 1 or 2 for Ignition Stage !' 10 ...
    'Input 1 for One-stage Ignition.' 10 ...
    'Input 2 for two-stage ignition.'] );

%% Folder£ºAssign and Creat Directory
currentFolder = pwd;

rcm.date = dataFilename(1:8);

dataFolder = [currentFolder, '\data\',rcm.date];

figFolder = [currentFolder, '\fig\',rcm.date];

rcmFolder = [currentFolder, '\rcm\',rcm.date];

xlsFolder = [currentFolder, '\xls'];

rcm.dataPath = [dataFolder, '\', dataFilename];

sampling_figFolder = [currentFolder, '\fig\samplingFig\',rcm.date];

if exist( figFolder, 'dir' )==0
    mkdir( figFolder )
end

if exist( rcmFolder, 'dir')==0
    mkdir( rcmFolder )
end

if exist( xlsFolder, 'dir' )==0
    mkdir(xlsFolder)
end

if exist( sampling_figFolder,'dir' )==0
    mkdir( sampling_figFolder )
end

%% Get date and time
listing = dir(rcm.dataPath);

% check we got a single entry corresponding to the file
assert(numel(listing) == 1, 'No such file: %s', rcm.dataPath);

rcm.raw.data = dlmread(rcm.dataPath);

if dataFilename(16) == '.'
    ID_map_str = [dataFilename(1:8),'0', dataFilename(15)];
    rcm.dayID = ['0', dataFilename(15)];
else
    ID_map_str = [dataFilename(1:8),dataFilename(15:16)];
    rcm.dayID = dataFilename(15:16);
end

rcm.ID = [rcm.date rcm.dayID];

%% Raw chamber pressure data
rcm.raw.length = size(rcm.raw.data, 1);

rcm.P.raw.time = 1/rcm.raw.length: 1/rcm.raw.length: 1;

rcm.P.raw.time = rcm.P.raw.time';

N_COLUMN = 1;            % Which column to be filted
                          % 1st column refers to the chamber pressure
rcm.P.post = ( (rcm.raw.data(:,N_COLUMN)-...
    sum(rcm.raw.data(1:20000,N_COLUMN)/20000))*rcm.P.AMP ) + rcm.initial.P;
                          % P.post refers to the chamber pressure calibrated with initial pressure
                          
%% Filter Chamber Pressure

% Plz comment on this section
                          % dp_dCA refers to rwa pressure trace
                          % cadArray refers to NUMBER of samples
                          % hrrCOF refers to filting frequence
cadArray= rcm.raw.data(:,N_COLUMN);
hrrCOF=rcm.P.Filter.hrrCOF;
NyquistFreq = rcm.P.Filter.NyquistFreq;
cof = hrrCOF/NyquistFreq;
dp_dCA= rcm.P.post;
fft_of_dpdca = fft (dp_dCA);
aaa = (0:1:(length(cadArray)-1)) * 0.8326 / (length(cadArray)/2*cof);
aaa = aaa';
filter_func = exp (-1 * aaa.^2);
rev_filter_func = filter_func (length(filter_func):-1:1);
filter_func = filter_func + rev_filter_func;
filter_limit = ones(size(filter_func));
filter_func = min (filter_func, filter_limit);
dp_dCA = ifft (fft_of_dpdca .* filter_func);
dp_dCA = real (dp_dCA);
rcm.P.post =dp_dCA;
                          % Diff P
rcm.P.dP = diff(rcm.P.post)./diff(rcm.P.raw.time);

%% LOOK UP TDC
                               % Beginning of Effective Pressure
[rcm.lookup.TDC.eff_begin_P, rcm.lookup.TDC.eff_begin_index] = ...
    max( rcm.P.post(rcm.lookup.TDC.eff_begin_lower+1: rcm.lookup.TDC.eff_begin_upper) );
                               
rcm.lookup.TDC.eff_begin_index = rcm.lookup.TDC.eff_begin_index +  rcm.lookup.TDC.eff_begin_lower;
                                % Beginning of Ignition delay
[rcm.lookup.TDC.begin_P, rcm.lookup.TDC.begin_index] = ...
    max( rcm.P.post(rcm.lookup.TDC.begin_lower+1:rcm.lookup.TDC.begin_upper) );

rcm.lookup.TDC.zero_index = rcm.lookup.TDC.begin_index +  rcm.lookup.TDC.begin_lower;
                
rcm.lookup.TDC.zero_time = rcm.lookup.TDC.zero_index/10^2; 
                                % Scale in [ms]

%% LOOK UP IGN
rcm.lookup.TDC.finish_lower = rcm.lookup.TDC.eff_begin_index;

if rcm.TDC_lower_bool == 0
    
    rcm.lookup.TDC.finish_index = rcm.lookup.TDC.finish_upper;

else
    
    [rcm.lookup.TDC.finish_P, rcm.lookup.TDC.finish_index] = ...
        min( rcm.P.post(rcm.lookup.TDC.finish_lower+1:rcm.lookup.TDC.finish_upper) );
    
    rcm.lookup.TDC.finish_index = rcm.lookup.TDC.finish_index + rcm.lookup.TDC.finish_lower;
end

rcm.P_TDC = rcm.P.post(rcm.lookup.TDC.finish_index);
                            % Pressure corresponding to the TDC point

rcm.P_eff = sum( rcm.P.post(rcm.lookup.TDC.eff_begin_index:rcm.lookup.TDC.finish_index) )...
                            /( rcm.lookup.TDC.finish_index - rcm.lookup.TDC.eff_begin_index +1 );
                            % Cal P_eff
rcm.P_ratio = rcm.P_eff / rcm.initial.P;
rcm.P_TDC_ratio = rcm.P_TDC / rcm.initial.P;
                            % Look up first ignition delay
                            % Mark first ignition delay in future version
[rcm.lookup.ign.first_flag, rcm.lookup.ign.first_index] = max( ...
    rcm.P.dP( rcm.lookup.ign.first_lower+1: rcm.lookup.ign.first_upper ) );

rcm.lookup.ign.first_index = rcm.lookup.ign.first_index + rcm.lookup.ign.first_lower;

rcm.ign.first_time = ( rcm.lookup.ign.first_index - rcm.lookup.TDC.zero_index +1)/100;

[rcm.lookup.ign.total_flag, rcm.lookup.ign.total_index] = max(...
    rcm.P.dP( rcm.lookup.ign.total_lower+1:rcm.lookup.ign_upper));

rcm.lookup.ign.total_index = rcm.lookup.ign.total_index + rcm.lookup.ign.total_lower;

rcm.ign.total_time = (rcm.lookup.ign.total_index - rcm.lookup.TDC.zero_index+1)/100;

rcm.ign.second_time = rcm.ign.total_time - rcm.ign.first_time;

rcm.ign_time = [rcm.ign.first_time,rcm.ign.second_time, rcm.ign.total_time];


%% Cal Effective Temperature % Volume Compression ratio
dT = 0.01;
dP_sum = 0;
                  %Universal gas constant
R=8.3145;        %J/(mol-K)
T = rcm.initial.T;
while dP_sum < log( rcm.P_ratio )
    %% Cal Cp   
    if( T>1000 )
        %Cp/R of N2  MW=28.01348
        N2=2.95257626+1.39690057E-3.*T-4.92631691E-7.*T.^2+7.86010367E-11.*T.^3-4.60755321E-15.*T.^4;
        %Cp/R of O2   MW=31.9988
        O2=3.66096083+6.56365523E-4.*T-1.41149485E-7.*T.^2+2.05797658E-11.*T.^3-1.29913248E-15.*T.^4;
        CO2=4.63659493+2.74131991E-3.*T-9.95828531E-7.*T.^2+1.60373011E-10.*T.^3-9.16103468E-15.*T.^4;
        H2O=2.67703787+2.97318329E-3.*T-7.73769690E-7.*T.^2+9.44336689E-11.*T.^3-4.26900959E-15.*T.^4;
        CO=3.04848583+1.35172818E-3.*T-4.85794075E-7.*T.^2+7.88536486E-11.*T.^3-4.69807489E-15.*T.^4;
        H2=3.33727920E+00-4.94024731E-05.*T+4.99456778E-07.*T.^2-1.79566394E-10.*T.^3+2.00255376E-14.*T.^4;
        CH4=7.48514950E-02+ 1.33909467E-02.*T-5.73285809E-06.*T^2 + 1.22292535E-09.*T^3 -1.01815230E-13.*T^4;
    else
        N2=3.53100528-1.23660987E-4.*T-5.02999437E-7.*T.^2+2.43530612E-09.*T.^3-1.40881235E-12.*T.^4;
        O2=3.78245636-2.99673415E-3.*T+9.84730200E-6.*T.^2-9.68129508E-09.*T.^3+3.24372836E-12.*T.^4;
        CO2=2.35677352+8.98459677E-3.*T-7.12356269E-6.*T.^2+2.45919022E-09.*T.^3-1.43699548E-13.*T.^4;
        H2O=4.19864056-2.03643410E-3.*T+6.52040211E-6.*T.^2-5.48797062E-09.*T.^3+1.77197817E-12.*T.^4;
        CO=3.57953347-6.10353680E-4.*T+1.01681433E-6.*T.^2+9.07005884E-10.*T.^3-9.04424499E-13.*T.^4;
        H2=2.34433112E+00+7.98052075E-03.*T-1.94781510E-05.*T.^2+2.01572094E-08.*T.^3-7.37611761E-12.*T.^4;
        CH4=5.14987613E+00 -1.36709788E-02.*T + 4.91800599E-05.*T^2 -4.84743026E-08.*T^3+ 1.66693956E-11.*T^4;
    end
    
    if (T>1391)
        n_Heptane = 2.22148969E+01+3.47675750E-02.*T-1.18407129E-05.*T.^2+1.83298478E-09.*T.^3-1.06130266E-13.*T.^4; %last number should be 4.
        MCH = 2.14785343E+01+3.32215917E-02.*T-1.14861934E-05.*T.^2+1.79638933E-09.*T.^3-1.04761864E-13.*T.^4; %last number should be 4.
    else
        n_Heptane =-1.26836187E+00+8.54355820E-02.*T-5.25346786E-05.*T.^2+1.62945721E-08.*T.^3-2.02394925E-12.*T.^4;
        MCH = -8.09426478E+00+1.00736150E-01.*T--7.00859796E-05.*T.^2+2.48687934E-08.*T.^3-3.59166681E-12.*T.^4; %last number should be 4.
    end

    if (T>1396)
        Toluene =   1.69989044e+01+2.20052825e-02.*T-7.58020314e-06.*T.^2+1.18267678e-09.*T.^3-6.88594262e-14.*T.^4;
        Isooctane = 2.71373590E+01+3.79004890E-02.*T-1.29437358E-05.*T.^2+2.00760372E-09.*T.^3-1.16400580E-13.*T.^4;
    else
        Toluene   =-5.45944164e+00+7.71089443e-02.*T-5.99305765e-05.*T.^2+2.40404364e-08.*T.^3-3.92116250e-12.*T.^4;
        Isooctane =-4.20868893E+00+1.11440581E-01.*T-7.91346582E-05.*T.^2+2.92406242E-08.*T.^3-4.43743191E-12.*T.^4;
    end
    
    if (T>1394)
        Isobutanol = 1.47423131e+01 + 2.19843767e-02.*T -7.51584192e-06.*T.^2 +1.16633393e-09.*T.^3 -6.76421638e-14.*T.^4;
    else
        Isobutanol =-8.37465362e-01 + 5.76520639e-02.*T - 3.90215462e-05.*T.^2 +1.40131231e-08.*T.^3-2.11159111e-12.*T.^4;
    end
    
    if (T>1393)
        DMH26 = 2.86146686E+01 + 4.41140780E-02.*T -1.51388773E-05.*T.^2 + 2.35538159E-09.*T.^3 -1.36850135E-13.*T.^4;
    else
        DMH26 = -2.69759199E+00 + 1.13719824E-01.*T -7.40347807E-05.*T.^2 +2.50653320E-08.*T.^3 -3.52291338E-12.*T.^4;
    end

    %Cp/R of Ar   MW=40
    Ar=2.5;
    
    Cp=(H2*rcm.Molar.H2 + N2*rcm.Molar.N2 + O2*rcm.Molar.O2 + Ar*rcm.Molar.Ar +...
        CO2*rcm.Molar.CO2 + H2O*rcm.Molar.H2O + CO*rcm.Molar.CO + ...
        n_Heptane*rcm.Molar.nHeptane + Toluene * rcm.Molar.Toluene +...
        Isooctane*rcm.Molar.Isooctane + MCH*rcm.Molar.MCH +CH4*rcm.Molar.CH4 +...
        Isobutanol*rcm.Molar.Isobutanol +DMH26*rcm.Molar.DMH26 )*R;
    Cv = Cp - R;
    gamma = Cp / Cv;
    
    dP_sum = dP_sum + (gamma/(gamma-1)/T)*dT;
    T = T + dT;
end
rcm.T_eff = T; 

                  %Cal Volume Compression ratio
CR_V_sum = 0;
T = rcm.initial.T;
N = 50000;
dT = (rcm.T_eff - rcm.initial.T)/N;
for i=1:N
    %% Cal Cp   
    if( T>1000 )
        %Cp/R of N2  MW=28.01348
        N2=2.95257626+1.39690057E-3.*T-4.92631691E-7.*T.^2+7.86010367E-11.*T.^3-4.60755321E-15.*T.^4;
        %Cp/R of O2   MW=31.9988
        O2=3.66096083+6.56365523E-4.*T-1.41149485E-7.*T.^2+2.05797658E-11.*T.^3-1.29913248E-15.*T.^4;
        CO2=4.63659493+2.74131991E-3.*T-9.95828531E-7.*T.^2+1.60373011E-10.*T.^3-9.16103468E-15.*T.^4;
        H2O=2.67703787+2.97318329E-3.*T-7.73769690E-7.*T.^2+9.44336689E-11.*T.^3-4.26900959E-15.*T.^4;
        CO=3.04848583+1.35172818E-3.*T-4.85794075E-7.*T.^2+7.88536486E-11.*T.^3-4.69807489E-15.*T.^4;
        H2=3.33727920E+00-4.94024731E-05.*T+4.99456778E-07.*T.^2-1.79566394E-10.*T.^3+2.00255376E-14.*T.^4;
        CH4=7.48514950E-02+ 1.33909467E-02.*T-5.73285809E-06.*T^2 + 1.22292535E-09.*T^3 -1.01815230E-13.*T^4;
    else
        N2=3.53100528-1.23660987E-4.*T-5.02999437E-7.*T.^2+2.43530612E-09.*T.^3-1.40881235E-12.*T.^4;
        O2=3.78245636-2.99673415E-3.*T+9.84730200E-6.*T.^2-9.68129508E-09.*T.^3+3.24372836E-12.*T.^4;
        CO2=2.35677352+8.98459677E-3.*T-7.12356269E-6.*T.^2+2.45919022E-09.*T.^3-1.43699548E-13.*T.^4;
        H2O=4.19864056-2.03643410E-3.*T+6.52040211E-6.*T.^2-5.48797062E-09.*T.^3+1.77197817E-12.*T.^4;
        CO=3.57953347-6.10353680E-4.*T+1.01681433E-6.*T.^2+9.07005884E-10.*T.^3-9.04424499E-13.*T.^4;
        H2=2.34433112E+00+7.98052075E-03.*T-1.94781510E-05.*T.^2+2.01572094E-08.*T.^3-7.37611761E-12.*T.^4;
        CH4=5.14987613E+00 -1.36709788E-02.*T + 4.91800599E-05.*T^2 -4.84743026E-08.*T^3+ 1.66693956E-11.*T^4;
    end
    
    if (T>1391)
        n_Heptane = 2.22148969E+01+3.47675750E-02.*T-1.18407129E-05.*T.^2+1.83298478E-09.*T.^3-1.06130266E-13.*T.^4; %last number should be 4.
        MCH = 2.14785343E+01+3.32215917E-02.*T-1.14861934E-05.*T.^2+1.79638933E-09.*T.^3-1.04761864E-13.*T.^4; %last number should be 4.
    else
        n_Heptane =-1.26836187E+00+8.54355820E-02.*T-5.25346786E-05.*T.^2+1.62945721E-08.*T.^3-2.02394925E-12.*T.^4;
        MCH = -8.09426478E+00+1.00736150E-01.*T--7.00859796E-05.*T.^2+2.48687934E-08.*T.^3-3.59166681E-12.*T.^4; %last number should be 4.
    end

    if (T>1396)
        Toluene =   1.69989044e+01+2.20052825e-02.*T-7.58020314e-06.*T.^2+1.18267678e-09.*T.^3-6.88594262e-14.*T.^4;
        Isooctane = 2.71373590E+01+3.79004890E-02.*T-1.29437358E-05.*T.^2+2.00760372E-09.*T.^3-1.16400580E-13.*T.^4;
    else
        Toluene   =-5.45944164e+00+7.71089443e-02.*T-5.99305765e-05.*T.^2+2.40404364e-08.*T.^3-3.92116250e-12.*T.^4;
        Isooctane =-4.20868893E+00+1.11440581E-01.*T-7.91346582E-05.*T.^2+2.92406242E-08.*T.^3-4.43743191E-12.*T.^4;
    end
    
    if (T>1394)
        Isobutanol = 1.47423131e+01 + 2.19843767e-02.*T -7.51584192e-06.*T.^2 +1.16633393e-09.*T.^3 -6.76421638e-14.*T.^4;
    else
        Isobutanol =-8.37465362e-01 + 5.76520639e-02.*T - 3.90215462e-05.*T.^2 +1.40131231e-08.*T.^3-2.11159111e-12.*T.^4;
    end
    
    if (T>1393)
        DMH26 = 2.86146686E+01 + 4.41140780E-02.*T -1.51388773E-05.*T.^2 + 2.35538159E-09.*T.^3 -1.36850135E-13.*T.^4;
    else
        DMH26 = -2.69759199E+00 + 1.13719824E-01.*T -7.40347807E-05.*T.^2 +2.50653320E-08.*T.^3 -3.52291338E-12.*T.^4;
    end

    %Cp/R of Ar   MW=40
    Ar=2.5;
    
    Cp=(H2*rcm.Molar.H2 + N2*rcm.Molar.N2 + O2*rcm.Molar.O2 + Ar*rcm.Molar.Ar +...
        CO2*rcm.Molar.CO2 + H2O*rcm.Molar.H2O + CO*rcm.Molar.CO + ...
        n_Heptane*rcm.Molar.nHeptane + Toluene * rcm.Molar.Toluene +...
        Isooctane*rcm.Molar.Isooctane + MCH*rcm.Molar.MCH +CH4*rcm.Molar.CH4 +...
        Isobutanol*rcm.Molar.Isobutanol +DMH26*rcm.Molar.DMH26 )*R;
    Cv = Cp - R;
    gamma = Cp / Cv;
    
    CR_V_sum = CR_V_sum + 1/(gamma-1)/T*dT;
    T = T + dT;
end
ln_CR_V = CR_V_sum;
rcm.V_CR = exp(ln_CR_V);

%% Cal TDC Temperature % Volume Compression ratio
dT = 0.01;
dP_sum = 0;
                  %Universal gas constant
R=8.3145;        %J/(mol-K)
T = rcm.initial.T;
while dP_sum < log( rcm.P_TDC_ratio )
    %% Cal Cp   
    if( T>1000 )
        %Cp/R of N2  MW=28.01348
        N2=2.95257626+1.39690057E-3.*T-4.92631691E-7.*T.^2+7.86010367E-11.*T.^3-4.60755321E-15.*T.^4;
        %Cp/R of O2   MW=31.9988
        O2=3.66096083+6.56365523E-4.*T-1.41149485E-7.*T.^2+2.05797658E-11.*T.^3-1.29913248E-15.*T.^4;
        CO2=4.63659493+2.74131991E-3.*T-9.95828531E-7.*T.^2+1.60373011E-10.*T.^3-9.16103468E-15.*T.^4;
        H2O=2.67703787+2.97318329E-3.*T-7.73769690E-7.*T.^2+9.44336689E-11.*T.^3-4.26900959E-15.*T.^4;
        CO=3.04848583+1.35172818E-3.*T-4.85794075E-7.*T.^2+7.88536486E-11.*T.^3-4.69807489E-15.*T.^4;
        H2=3.33727920E+00-4.94024731E-05.*T+4.99456778E-07.*T.^2-1.79566394E-10.*T.^3+2.00255376E-14.*T.^4;
        CH4=7.48514950E-02+ 1.33909467E-02.*T-5.73285809E-06.*T^2 + 1.22292535E-09.*T^3 -1.01815230E-13.*T^4;
    else
        N2=3.53100528-1.23660987E-4.*T-5.02999437E-7.*T.^2+2.43530612E-09.*T.^3-1.40881235E-12.*T.^4;
        O2=3.78245636-2.99673415E-3.*T+9.84730200E-6.*T.^2-9.68129508E-09.*T.^3+3.24372836E-12.*T.^4;
        CO2=2.35677352+8.98459677E-3.*T-7.12356269E-6.*T.^2+2.45919022E-09.*T.^3-1.43699548E-13.*T.^4;
        H2O=4.19864056-2.03643410E-3.*T+6.52040211E-6.*T.^2-5.48797062E-09.*T.^3+1.77197817E-12.*T.^4;
        CO=3.57953347-6.10353680E-4.*T+1.01681433E-6.*T.^2+9.07005884E-10.*T.^3-9.04424499E-13.*T.^4;
        H2=2.34433112E+00+7.98052075E-03.*T-1.94781510E-05.*T.^2+2.01572094E-08.*T.^3-7.37611761E-12.*T.^4;
        CH4=5.14987613E+00 -1.36709788E-02.*T + 4.91800599E-05.*T^2 -4.84743026E-08.*T^3+ 1.66693956E-11.*T^4;
    end
    
    if (T>1391)
        n_Heptane = 2.22148969E+01+3.47675750E-02.*T-1.18407129E-05.*T.^2+1.83298478E-09.*T.^3-1.06130266E-13.*T.^4; %last number should be 4.
        MCH = 2.14785343E+01+3.32215917E-02.*T-1.14861934E-05.*T.^2+1.79638933E-09.*T.^3-1.04761864E-13.*T.^4; %last number should be 4.
    else
        n_Heptane =-1.26836187E+00+8.54355820E-02.*T-5.25346786E-05.*T.^2+1.62945721E-08.*T.^3-2.02394925E-12.*T.^4;
        MCH = -8.09426478E+00+1.00736150E-01.*T--7.00859796E-05.*T.^2+2.48687934E-08.*T.^3-3.59166681E-12.*T.^4; %last number should be 4.
    end

    if (T>1396)
        Toluene =   1.69989044e+01+2.20052825e-02.*T-7.58020314e-06.*T.^2+1.18267678e-09.*T.^3-6.88594262e-14.*T.^4;
        Isooctane = 2.71373590E+01+3.79004890E-02.*T-1.29437358E-05.*T.^2+2.00760372E-09.*T.^3-1.16400580E-13.*T.^4;
    else
        Toluene   =-5.45944164e+00+7.71089443e-02.*T-5.99305765e-05.*T.^2+2.40404364e-08.*T.^3-3.92116250e-12.*T.^4;
        Isooctane =-4.20868893E+00+1.11440581E-01.*T-7.91346582E-05.*T.^2+2.92406242E-08.*T.^3-4.43743191E-12.*T.^4;
    end
    
    if (T>1394)
        Isobutanol = 1.47423131e+01 + 2.19843767e-02.*T -7.51584192e-06.*T.^2 +1.16633393e-09.*T.^3 -6.76421638e-14.*T.^4;
    else
        Isobutanol =-8.37465362e-01 + 5.76520639e-02.*T - 3.90215462e-05.*T.^2 +1.40131231e-08.*T.^3-2.11159111e-12.*T.^4;
    end
    
    if (T>1393)
        DMH26 = 2.86146686E+01 + 4.41140780E-02.*T -1.51388773E-05.*T.^2 + 2.35538159E-09.*T.^3 -1.36850135E-13.*T.^4;
    else
        DMH26 = -2.69759199E+00 + 1.13719824E-01.*T -7.40347807E-05.*T.^2 +2.50653320E-08.*T.^3 -3.52291338E-12.*T.^4;
    end

    %Cp/R of Ar   MW=40
    Ar=2.5;
    
    Cp=(H2*rcm.Molar.H2 + N2*rcm.Molar.N2 + O2*rcm.Molar.O2 + Ar*rcm.Molar.Ar +...
        CO2*rcm.Molar.CO2 + H2O*rcm.Molar.H2O + CO*rcm.Molar.CO + ...
        n_Heptane*rcm.Molar.nHeptane + Toluene * rcm.Molar.Toluene +...
        Isooctane*rcm.Molar.Isooctane + MCH*rcm.Molar.MCH +CH4*rcm.Molar.CH4 +...
        Isobutanol*rcm.Molar.Isobutanol +DMH26*rcm.Molar.DMH26 )*R;
    Cv = Cp - R;
    gamma = Cp / Cv;
    
    dP_sum = dP_sum + (gamma/(gamma-1)/T)*dT;
    T = T + dT;
end
rcm.T_TDC_eff = T; 

                  %Cal Volume Compression ratio
CR_V_sum = 0;
T = rcm.initial.T;
N = 50000;
dT = (rcm.T_TDC_eff - rcm.initial.T)/N;
for i=1:N
    %% Cal Cp   
    if( T>1000 )
        %Cp/R of N2  MW=28.01348
        N2=2.95257626+1.39690057E-3.*T-4.92631691E-7.*T.^2+7.86010367E-11.*T.^3-4.60755321E-15.*T.^4;
        %Cp/R of O2   MW=31.9988
        O2=3.66096083+6.56365523E-4.*T-1.41149485E-7.*T.^2+2.05797658E-11.*T.^3-1.29913248E-15.*T.^4;
        CO2=4.63659493+2.74131991E-3.*T-9.95828531E-7.*T.^2+1.60373011E-10.*T.^3-9.16103468E-15.*T.^4;
        H2O=2.67703787+2.97318329E-3.*T-7.73769690E-7.*T.^2+9.44336689E-11.*T.^3-4.26900959E-15.*T.^4;
        CO=3.04848583+1.35172818E-3.*T-4.85794075E-7.*T.^2+7.88536486E-11.*T.^3-4.69807489E-15.*T.^4;
        H2=3.33727920E+00-4.94024731E-05.*T+4.99456778E-07.*T.^2-1.79566394E-10.*T.^3+2.00255376E-14.*T.^4;
        CH4=7.48514950E-02+ 1.33909467E-02.*T-5.73285809E-06.*T^2 + 1.22292535E-09.*T^3 -1.01815230E-13.*T^4;
    else
        N2=3.53100528-1.23660987E-4.*T-5.02999437E-7.*T.^2+2.43530612E-09.*T.^3-1.40881235E-12.*T.^4;
        O2=3.78245636-2.99673415E-3.*T+9.84730200E-6.*T.^2-9.68129508E-09.*T.^3+3.24372836E-12.*T.^4;
        CO2=2.35677352+8.98459677E-3.*T-7.12356269E-6.*T.^2+2.45919022E-09.*T.^3-1.43699548E-13.*T.^4;
        H2O=4.19864056-2.03643410E-3.*T+6.52040211E-6.*T.^2-5.48797062E-09.*T.^3+1.77197817E-12.*T.^4;
        CO=3.57953347-6.10353680E-4.*T+1.01681433E-6.*T.^2+9.07005884E-10.*T.^3-9.04424499E-13.*T.^4;
        H2=2.34433112E+00+7.98052075E-03.*T-1.94781510E-05.*T.^2+2.01572094E-08.*T.^3-7.37611761E-12.*T.^4;
        CH4=5.14987613E+00 -1.36709788E-02.*T + 4.91800599E-05.*T^2 -4.84743026E-08.*T^3+ 1.66693956E-11.*T^4;
    end
    
    if (T>1391)
        n_Heptane = 2.22148969E+01+3.47675750E-02.*T-1.18407129E-05.*T.^2+1.83298478E-09.*T.^3-1.06130266E-13.*T.^4; %last number should be 4.
        MCH = 2.14785343E+01+3.32215917E-02.*T-1.14861934E-05.*T.^2+1.79638933E-09.*T.^3-1.04761864E-13.*T.^4; %last number should be 4.
    else
        n_Heptane =-1.26836187E+00+8.54355820E-02.*T-5.25346786E-05.*T.^2+1.62945721E-08.*T.^3-2.02394925E-12.*T.^4;
        MCH = -8.09426478E+00+1.00736150E-01.*T--7.00859796E-05.*T.^2+2.48687934E-08.*T.^3-3.59166681E-12.*T.^4; %last number should be 4.
    end

    if (T>1396)
        Toluene =   1.69989044e+01+2.20052825e-02.*T-7.58020314e-06.*T.^2+1.18267678e-09.*T.^3-6.88594262e-14.*T.^4;
        Isooctane = 2.71373590E+01+3.79004890E-02.*T-1.29437358E-05.*T.^2+2.00760372E-09.*T.^3-1.16400580E-13.*T.^4;
    else
        Toluene   =-5.45944164e+00+7.71089443e-02.*T-5.99305765e-05.*T.^2+2.40404364e-08.*T.^3-3.92116250e-12.*T.^4;
        Isooctane =-4.20868893E+00+1.11440581E-01.*T-7.91346582E-05.*T.^2+2.92406242E-08.*T.^3-4.43743191E-12.*T.^4;
    end
    
    if (T>1394)
        Isobutanol = 1.47423131e+01 + 2.19843767e-02.*T -7.51584192e-06.*T.^2 +1.16633393e-09.*T.^3 -6.76421638e-14.*T.^4;
    else
        Isobutanol =-8.37465362e-01 + 5.76520639e-02.*T - 3.90215462e-05.*T.^2 +1.40131231e-08.*T.^3-2.11159111e-12.*T.^4;
    end
    
    if (T>1393)
        DMH26 = 2.86146686E+01 + 4.41140780E-02.*T -1.51388773E-05.*T.^2 + 2.35538159E-09.*T.^3 -1.36850135E-13.*T.^4;
    else
        DMH26 = -2.69759199E+00 + 1.13719824E-01.*T -7.40347807E-05.*T.^2 +2.50653320E-08.*T.^3 -3.52291338E-12.*T.^4;
    end

    %Cp/R of Ar   MW=40
    Ar=2.5;
    
    Cp=(H2*rcm.Molar.H2 + N2*rcm.Molar.N2 + O2*rcm.Molar.O2 + Ar*rcm.Molar.Ar +...
        CO2*rcm.Molar.CO2 + H2O*rcm.Molar.H2O + CO*rcm.Molar.CO + ...
        n_Heptane*rcm.Molar.nHeptane + Toluene * rcm.Molar.Toluene +...
        Isooctane*rcm.Molar.Isooctane + MCH*rcm.Molar.MCH +CH4*rcm.Molar.CH4 +...
        Isobutanol*rcm.Molar.Isobutanol +DMH26*rcm.Molar.DMH26 )*R;
    Cv = Cp - R;
    gamma = Cp / Cv;
    
    CR_V_sum = CR_V_sum + 1/(gamma-1)/T*dT;
    T = T + dT;
end
ln_CR_V = CR_V_sum;
rcm.TDC_V_CR = exp(ln_CR_V);


%% PLOT Chamber Pressure Trace and Differential Trace & Save Fig

figure(1);

[AX, H1, H2] = plotyy(  rcm.P.raw.time,                         rcm.P.post,...
                        rcm.P.raw.time(1:rcm.raw.length-1),     rcm.P.dP    );

set(get(AX(1), 'Ylabel'), 'String', 'Pressure(bar)', 'FontSize', FontSize, ...
    'FontWeight', FontWeight);
set(get(AX(2), 'Ylabel'), 'String', 'dP', 'FontSize', FontSize, ...
    'FontWeight', FontWeight);

set(AX(1), 'fontsize', FontSize);
set(AX(2), 'fontsize', FontSize);

set(H1, 'Color', 'r', 'LineWidth', LineWidth, 'LineStyle', '-');
set(H2, 'Color', 'b', 'LineWidth', LineWidth, 'LineStyle', '--');

xlabel('Time[s]', 'FontSize', FontSize, 'FontWeight', FontWeight);

legend ( 'Presure','dP','Location','NorthWest');

title(rcm.ID, 'FontSize', FontSize, 'FontWeight', FontWeight);

hold all;

plot( [rcm.lookup.TDC.zero_time rcm.lookup.TDC.zero_time]./1E3,...
    [min( rcm.P.post ) max(rcm.P.post)] ,'-.');
                                        % Plot Ignition delay start
hold all;

plot( [rcm.lookup.TDC.zero_time+rcm.ign.total_time ...
    rcm.lookup.TDC.zero_time+rcm.ign.total_time]./1E3,...
    [min( rcm.P.post ) max(rcm.P.post)] ,'-.k');
                                        % Plot ignition
hold all;

plot( [0 1],[rcm.P.post(rcm.lookup.TDC.eff_begin_index) ...
    rcm.P.post(rcm.lookup.TDC.eff_begin_index)] ,'-.');
                                        % Plot P_eff start     
hold all;

plot( [0 1],[rcm.P.post(rcm.lookup.TDC.finish_index) ...
    rcm.P.post(rcm.lookup.TDC.finish_index)] ,'-.');

hold all;

str = [ 'Peff = ' num2str(rcm.P_eff ,'%.2f')  '  Teff = ' num2str(rcm.T_eff,'%.2f') ...
    '  Ign = ' num2str(rcm.ign.total_time,'%.2f') ...
    '  CR = ' num2str([rcm.V_CR rcm.P.post(rcm.lookup.TDC.eff_begin_index) ], '%.2f\n')...
    ' - ' num2str( rcm.P.post(rcm.lookup.TDC.finish_index),'%.2f' ) ' = ' ...
    num2str( rcm.P.post(rcm.lookup.TDC.eff_begin_index) - rcm.P.post(rcm.lookup.TDC.finish_index),'%.2f' )] ;
text('Position',[rcm.lookup.TDC.zero_time./1E3 max(rcm.P.post)], ...
    'String',str,...
    'FontName', 'Times New Roman','FontSize', 8);
                                        % Plot P_eff end
FileFig1 = ['RCM_',rcm.ID,'.fig']; 

saveas(figure(1), fullfile(figFolder,FileFig1));
saveas(figure(1), fullfile([currentFolder, '\fig'],FileFig1));

while ( rcm.bool_sampling == 1 )
%% Sampling Data Process
N_COLUMN = 2;

rcm.sampling.P = (rcm.raw.data(:,N_COLUMN)-...
    sum(rcm.raw.data(1:20000,N_COLUMN)/20000))*rcm.sampling.P_AMP;

%% Filter Sampling Pressure
    
cadArray= rcm.raw.data(:,N_COLUMN);
hrrCOF=rcm.sampling.Filter.hrrCOF;
NyquistFreq = rcm.sampling.Filter.NyquistFreq;
cof = hrrCOF/NyquistFreq;
dp_dCA = rcm.sampling.P;
fft_of_dpdca = fft (dp_dCA);
aaa = (0:1:(length(cadArray)-1)) * 0.8326 / (length(cadArray)/2*cof);
aaa = aaa';
filter_func = exp (-1 * aaa.^2);
rev_filter_func = filter_func (length(filter_func):-1:1);
filter_func = filter_func + rev_filter_func;
filter_limit = ones(size(filter_func));
filter_func = min (filter_func, filter_limit);
dp_dCA = ifft (fft_of_dpdca .* filter_func);
dp_dCA = real (dp_dCA);
rcm.sampling.P =dp_dCA;

%% Cal sampling time
                  % Cal effective P_sampling, corrsponding to sampling amount
N_COLUMN = 3;
                  % PULSE is 5 V TTL
rcm.sampling.start_index = find(rcm.raw.data(:,N_COLUMN) > 0.2, 1 );

rcm.sampling.finish_index = find(rcm.raw.data(:,N_COLUMN) > 0.2, 1 );

rcm.sampling.middle_index = (rcm.sampling.start_index + rcm.sampling.finish_index)./2;

rcm.sampling.P_before = sum( rcm.sampling.P(1000:21000) )/(21000-1000);

rcm.sampling.P_after = sum( rcm.sampling.P(100000-21000:100000-1000) )/20000;

rcm.sampling.deltaP = rcm.sampling.P_after - rcm.sampling.P_before;

sampling_P_mid = (rcm.sampling.P_before + rcm.sampling.P_after)/2;

rcm.sampling.index = find( rcm.sampling.P(20000:80000) > sampling_P_mid ,1)...
    +20000;

rcm.sampling.time = rcm.P.raw.time(rcm.sampling.index)*10^3 ...
    - rcm.lookup.TDC.zero_time;

rcm.sampling.time_pulse = rcm.P.raw.time(rcm.sampling.middle_index)*10^3 ...
    - rcm.lookup.TDC.zero_time;

rcm.sampling.normal_time = rcm.sampling.time./rcm.ign_time(3);

%% PLOT Sampling Pressure Trace and Chamber Pressure Trace IF Sampling & Save Fig

figure(2);

[AX, H1, H2] = plotyy(  rcm.P.raw.time,     rcm.P.post,...
                        rcm.P.raw.time,     rcm.sampling.P  );

set(get(AX(1), 'Ylabel'), 'String', 'Chamber Pressure(bar)', ...
    'FontSize', FontSize, 'FontWeight', FontWeight, 'Color', 'k');
set(get(AX(2), 'Ylabel'), 'String', 'Sampling Pressure(bar)', ...
    'FontSize', FontSize, 'FontWeight', FontWeight, 'Color', 'b');

set(AX(1), 'fontsize', FontSize);
set(AX(2), 'fontsize', FontSize);

set(H1, 'Color', 'r', 'LineWidth', LineWidth, 'LineStyle', '-');
set(H2, 'Color', 'b', 'LineWidth', LineWidth, 'LineStyle', '--');

xlabel('Time[s]', 'FontSize', FontSize, 'FontWeight', FontWeight);

legend ( 'Chamber Pressure','Sampling Pressure','Location','NorthWest');

title([rcm.ID,'Sampling'], 'FontSize', FontSize, 'FontWeight', FontWeight);

hold all;

plot( [rcm.lookup.TDC.zero_time+rcm.sampling.time ...
    rcm.lookup.TDC.zero_time+rcm.sampling.time]./1E3,...
    [min( rcm.P.post ) max(rcm.P.post)] );

hold all;

plot( rcm.P.raw.time, rcm.raw.data(:,3).*max(rcm.P.post)./5 );

hold all;

plot( [rcm.lookup.TDC.zero_time rcm.lookup.TDC.zero_time]./1E3,...
    [min( rcm.P.post ) max(rcm.P.post)] ,'-.');

hold all;

plot( [rcm.lookup.TDC.zero_time+rcm.ign.total_time ...
    rcm.lookup.TDC.zero_time+rcm.ign.total_time]./1E3,...
    [min( rcm.P.post ) max(rcm.P.post)] ,'-.k');

hold all;

str = [ 'S NT = ' num2str(rcm.sampling.normal_time ,'%.2f') ...
    '  S dP = ' num2str(rcm.sampling.deltaP,'%.2f') ...
    '  S time = ' num2str(rcm.sampling.time,'%.2f') ];
text('Position',[rcm.lookup.TDC.zero_time./1E3 max(rcm.P.post)], ...
    'String',str,...
    'FontName', 'Times New Roman','FontSize', 8);

FileFig2 = ['RCM_Sampling_',rcm.ID,'.fig']; 

saveas(figure(2), fullfile(sampling_figFolder,FileFig2));
saveas(figure(2), fullfile([currentFolder, '\fig'],FileFig2));

break;

end

%% Save Data to Mat File
if processible == 1     % The current run could be processed.

rcm.P.raw_time = rcm.P.raw.time*1000;

RCM.ChamberPressure.shift = [rcm.P.raw_time-rcm.lookup.TDC.zero_time, rcm.P.post];

RCM.ChamberPressure.post = [rcm.P.raw_time, rcm.P.post];

RCM.ChamberPressure.raw = [rcm.P.raw_time, rcm.raw.data(:,1)];

RCM.ID = rcm.ID;

RCM.RunNum = rcm.dayID;

RCM.Fuel = rcm.fuel;

RCM.Phi = rcm.phi;

RCM.Dilution = rcm.Molar.Dilution_ratio;

RCM.p_initial = rcm.initial.P;

RCM.T_initial = rcm.initial.T;

RCM.p_eff = rcm.P_eff;

RCM.T_eff = rcm.T_eff;

RCM.TDC.P = rcm.P_TDC;
RCM.TDC.T = rcm.T_TDC_eff;
RCM.TDC.CR = rcm.TDC_V_CR;

if rcm.ign_stage == 1

    RCM.Ign_delay = rcm.ign_time(3);

elseif rcm.ign_stage == 2

    RCM.Ign_delay = [rcm.ign_time(1), rcm.ign_time(3)];

end

RCM.CR = rcm.V_CR;

RCM.ChamberLength = rcm.config;

RCM.Remark = Remark;

RCM.Molar = rcm.Molar;

RCM.ChamberPressure.Filter = rcm.P.Filter;

RCM.ChamberPressure.Filter.AMP = rcm.P.AMP;

if rcm.bool_sampling == 1
    RCM.Sampling = rcm.sampling;
elseif rcm.bool_sampling == 0
    RCM.Sampling = 'No Sampling';
end

RCM.lookup = rcm.lookup;

savemat = ['RCM_', rcm.ID, '.mat' ];

save(fullfile(rcmFolder, savemat), 'RCM');
save(fullfile([currentFolder, '\rcm'], savemat), 'RCM');

%% Print

fprintf( 'Effective Pressure:        P_eff = %.2f bar\n',rcm.P_eff);  

fprintf( 'Effective Temperature      T_eff = %.2f K\n',rcm.T_eff); 

fprintf( 'Volume Compression Ratio    V_CR = %.2f \n',rcm.V_CR); 

fprintf( 'ign delay time:             t_id = %.2f ms\n',RCM.Ign_delay); 
disp('TDC')
fprintf('P = %.2f bar\t T = %.2f K\t TDC_CR=%.2f\n',RCM.TDC.P,RCM.TDC.T,RCM.TDC.CR);
disp('Sampling')
fprintf('S_dP = %.2f bar\t S_t = %.2f ms\t S_NT=%.2f\n',RCM.Sampling.deltaP,RCM.Sampling.time,...
    RCM.Sampling.normal_time);

%% Save Data to Excel
Todayxls = fullfile(xlsFolder, [rcm.date, '.xlsx']);

ROW = str2num(rcm.dayID) + 1;

headers = { 'Date','DayID','Fuel',...
    rcm.fuel, 'O2','N2','Ar','CO2', 'PHI','Dilution Ratio',...
    'P_initial(bar)', 'T_initial', 'P_eff(bar)', 'T_eff(K)',...
    'ign_1st(ms)' ,'ign delay(ms)', 'Chamber length(mm)', 'Chamber Parts(mm)' ...
    'CR' , 'S_dP','S_Time(ms)','S_NT','Remarks' };

while ( rcm.bool_sampling == 0 )
if rcm.ign_stage == 1
    A_row = [ {rcm.date}, {rcm.dayID},{rcm.fuel},... 
        eval(['rcm.Molar.' rcm.fuel]), rcm.Molar.O2, rcm.Molar.N2,...
        rcm.Molar.Ar, rcm.Molar.CO2, rcm.phi, round(rcm.Molar.Dilution_ratio*100)/100, ...
        round(rcm.initial.P*10000)/10000, round(rcm.initial.T*100)/100, round(rcm.P_eff*100)/100, round(rcm.T_eff*100)/100, ...
        ' ', RCM.Ign_delay(1),eval(rcm.config),{rcm.config},...
        round(rcm.V_CR*1000)/1000, ...
        ' ',...
        ' ',...
        ' ', ...
        Remark ]; 

elseif rcm.ign_stage == 2
    A_row = [ {rcm.date}, {rcm.dayID},{rcm.fuel},... 
        eval(['rcm.Molar.' rcm.fuel]), rcm.Molar.O2, rcm.Molar.N2,...
        rcm.Molar.Ar, rcm.Molar.CO2, rcm.phi, round(rcm.Molar.Dilution_ratio*100)/100, ...
        round(rcm.initial.P*10000)/10000, round(rcm.initial.T*100)/100, round(rcm.P_eff*100)/100, round(rcm.T_eff*100)/100, ...
        RCM.Ign_delay(1), RCM.Ign_delay(2),eval(rcm.config),{rcm.config},...
        {rcm.config},round(rcm.V_CR*1000)/1000,...
        ' ',...
        ' ',...
        ' ',...
        Remark ]; 
end

xlRange = [ 'A',num2str(ROW) ];

xlswrite(Todayxls,headers,'Sheet1','A1');

xlswrite(Todayxls,A_row,'Sheet1',xlRange);

break;
end

while ( rcm.bool_sampling == 1 )
if rcm.ign_stage == 1
    A_row = [ {rcm.date}, {rcm.dayID},{rcm.fuel},... 
        eval(['rcm.Molar.' rcm.fuel]), rcm.Molar.O2, rcm.Molar.N2,...
        rcm.Molar.Ar, rcm.Molar.CO2, rcm.phi, round(rcm.Molar.Dilution_ratio*100)/100, ...
        round(rcm.initial.P*10000)/10000, round(rcm.initial.T*100)/100, round(rcm.P_eff*100)/100, round(rcm.T_eff*100)/100, ...
        ' ', RCM.Ign_delay(1),eval(rcm.config),{rcm.config},...
        round(rcm.V_CR*1000)/1000, ...
        rcm.sampling.deltaP,...
        rcm.sampling.time,...
        rcm.sampling.normal_time, ...
        Remark ]; 

elseif rcm.ign_stage == 2
    A_row = [ {rcm.date}, {rcm.dayID},{rcm.fuel},... 
        eval(['rcm.Molar.' rcm.fuel]), rcm.Molar.O2, rcm.Molar.N2,...
        rcm.Molar.Ar, rcm.Molar.CO2, rcm.phi, round(rcm.Molar.Dilution_ratio*100)/100, ...
        round(rcm.initial.P*10000)/10000, round(rcm.initial.T*100)/100, round(rcm.P_eff*100)/100, round(rcm.T_eff*100)/100, ...
        RCM.Ign_delay(1), RCM.Ign_delay(2),eval(rcm.config),{rcm.config},...
        {rcm.config},round(rcm.V_CR*1000)/1000,...
        rcm.sampling.deltaP,...
        rcm.sampling.time,...
        rcm.sampling.normal_time, ...
        Remark ]; 
end

xlRange = [ 'A',num2str(ROW) ];

xlswrite(Todayxls,headers,'Sheet1','A1');

xlswrite(Todayxls,A_row,'Sheet1',xlRange);

break;
end

elseif processible == 0     % The Current Run could not be processed.

rcm.P.raw_time = rcm.P.raw.time*1000;

RCM.ChamberPressure.post = [rcm.P.raw_time, rcm.P.post];

RCM.ChamberPressure.raw = [rcm.P.raw_time, rcm.raw.data(:,1)];

RCM.ID = rcm.ID;

RCM.RunNum = rcm.dayID;

RCM.Fuel = rcm.fuel;

RCM.Phi = rcm.phi;

RCM.Dilution = rcm.Molar.Dilution_ratio;

RCM.p_initial = rcm.initial.P;

RCM.T_initial = rcm.initial.T;

RCM.ChamberLength = rcm.config;

RCM.Remark = Remark;

RCM.Molar = rcm.Molar;

RCM.TDC.P = rcm.P_TDC;
RCM.TDC.T = rcm.T_TDC_eff;
RCM.TDC.CR = rcm.TDC_V_CR;

RCM.ChamberPressure.raw = rcm.P.raw;

RCM.ChamberPressure.post = rcm.P.post;

RCM.ChamberPressure.Filter = rcm.P.Filter;

RCM.ChamberPressure.Filter.AMP = rcm.P.AMP;

savemat = ['RCM_', rcm.ID, '.mat' ];

save(fullfile(rcmFolder, savemat), 'RCM');

fprintf( ['The Current Run Can not be Processed!' 10 'Can not Decide End of Compression!\n']);  

Todayxls = fullfile(xlsFolder, [rcm.date, '.xlsx']);

ROW = str2num(rcm.dayID) + 1;

headers = { 'Date','DayID','Fuel',...
    rcm.fuel, 'O2','N2','Ar','CO2', 'PHI','Dilution Ratio',...
    'P_initial(bar)', 'T_initial', 'P_eff(bar)', 'T_eff(K)',...
    'ign_1st(ms)' ,'ign delay(ms)', 'Chamber length(mm)', 'Chamber Parts(mm)' ...
    'CR' ,'Remarks' };

A_row = [ {rcm.date}, {rcm.dayID},{rcm.fuel},... 
	eval(['rcm.Molar.' rcm.fuel]), rcm.Molar.O2, rcm.Molar.N2,...
	rcm.Molar.Ar, rcm.Molar.CO2, rcm.phi, round(rcm.Molar.Dilution_ratio*100)/100, ...
	round(rcm.initial.P*10000)/10000, round(rcm.initial.T*100)/100, no_data_sign, no_data_sign, ...
	no_data_sign, no_data_sign,eval(rcm.config),{rcm.config},...
	no_data_sign, Remark ]; 

xlRange = [ 'A',num2str(ROW) ];

xlswrite(Todayxls,headers,'Sheet1','A1');

xlswrite(Todayxls,A_row,'Sheet1',xlRange);

end
% End of processible flag.





