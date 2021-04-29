function [results, success, raw] = dcopf_solver(om, mpopt)
% tic
%DCOPF_SOLVER  Solves a DC optimal power flow.
%
%   [RESULTS, SUCCESS, RAW] = DCOPF_SOLVER(OM, MPOPT)
%
%   Inputs are an OPF model object and a MATPOWER options struct.
%
%   Outputs are a RESULTS struct, SUCCESS flag and RAW output struct.
%
%   RESULTS is a MATPOWER case struct (mpc) with the usual baseMVA, bus
%   branch, gen, gencost fields, along with the following additional
%   fields:
%       .order      see 'help ext2int' for details of this field
%       .x          final value of optimization variables (internal order)
%       .f          final objective function value
%       .mu         shadow prices on ...
%           .var
%               .l  lower bounds on variables
%               .u  upper bounds on variables
%           .lin
%               .l  lower bounds on linear constraints
%               .u  upper bounds on linear constraints
%
%   SUCCESS     1 if solver converged successfully, 0 otherwise
%
%   RAW         raw output in form returned by MINOS
%       .xr     final value of optimization variables
%       .pimul  constraint multipliers
%       .info   solver specific termination code
%       .output solver specific output information
%
%   See also OPF, QPS_MATPOWER.

%   MATPOWER
%   Copyright (c) 2000-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%----- initialization -----
%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

%% unpack data
mpc = get_mpc(om);
[baseMVA, bus, gen, branch, gencost] = ...
    deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch, mpc.gencost);
cp = get_cost_params(om);
% gen
[N, H, Cw] = deal(cp.N, cp.H, cp.Cw);
fparm = [cp.dd cp.rh cp.kk cp.mm];
Bf = userdata(om, 'Bf');
Pfinj = userdata(om, 'Pfinj');
[vv, ll] = get_idx(om);

%% problem dimensions
ipol = find(gencost(:, MODEL) == POLYNOMIAL); %% polynomial costs
ipwl = find(gencost(:, MODEL) == PW_LINEAR);  %% piece-wise linear costs
igen = find(bus(:, PD) > 0); %% generator bus
nb = size(bus, 1);          %% number of buses
nl = size(branch, 1);       %% number of branches
ng = size(gen, 1);          %% number of generator
nw = size(N, 1);            %% number of general cost vars, w
ny = getN(om, 'var', 'y');  %% number of piece-wise linear costs
nxyz = getN(om, 'var');     %% total number of control vars of all types
iload = find(bus(:,PD));

%% linear constraints & variable bounds
[A, l, u] = linear_constraints(om);
[x0, xmin, xmax] = getv(om);

%% set up objective function of the form: f = 1/2 * X'*HH*X + CC'*X
%% where X = [x;y;z]. First set up as quadratic function of w,
%% f = 1/2 * w'*HHw*w + CCw'*w, where w = diag(M) * (N*X - Rhat). We
%% will be building on the (optionally present) user supplied parameters.

%% piece-wise linear costs
any_pwl = (ny > 0);
if any_pwl
    Npwl = sparse(ones(ny,1), vv.i1.y:vv.iN.y, 1, 1, nxyz);     %% sum of y vars
    Hpwl = 0;
    Cpwl = 1;
    fparm_pwl = [1 0 0 1];
else
    Npwl = sparse(0, nxyz);
    Hpwl = [];
    Cpwl = [];
    fparm_pwl = [];
end

%% quadratic costs
npol = length(ipol);
if any(find(gencost(ipol, NCOST) > 3))
    error('DC opf cannot handle polynomial costs with higher than quadratic order.');
end
iqdr = find(gencost(ipol, NCOST) == 3);
ilin = find(gencost(ipol, NCOST) == 2);
polycf = zeros(npol, 3);                            %% quadratic coeffs for Pg
if ~isempty(iqdr)
  polycf(iqdr, :)   = gencost(ipol(iqdr), COST:COST+2);
end
polycf(ilin, 2:3) = gencost(ipol(ilin), COST:COST+1);
polycf = polycf * diag([ baseMVA^2 baseMVA 1]);     %% convert to p.u.
Npol = sparse(1:npol, vv.i1.Pg-1+ipol, 1, npol, nxyz);         %% Pg vars
Hpol = sparse(1:npol, 1:npol, 2*polycf(:, 1), npol, npol);
Cpol = polycf(:, 2);
fparm_pol = ones(npol,1) * [ 1 0 0 1 ];

%% combine with user costs
NN = [ Npwl; Npol; N ];
HHw = [ Hpwl, sparse(any_pwl, npol+nw);
        sparse(npol, any_pwl), Hpol, sparse(npol, nw);
        sparse(nw, any_pwl+npol), H   ];
CCw = [Cpwl; Cpol; Cw];
ffparm = [ fparm_pwl; fparm_pol; fparm ];

%% transform quadratic coefficients for w into coefficients for X
nnw = any_pwl+npol+nw;
M   = sparse(1:nnw, 1:nnw, ffparm(:, 4), nnw, nnw);
MR  = M * ffparm(:, 2);
HMR = HHw * MR;
MN  = M * NN;
HH = MN' * HHw * MN;
CC = full(MN' * (CCw - HMR));
C0 = 1/2 * MR' * HMR + sum(polycf(:, 3));   %% constant term of cost

%% options
if isempty(HH) || ~any(any(HH))
    model = 'LP';
else
    model = 'QP';
end
opt = mpopt2qpopt(mpopt, model);

%% try to select an interior initial point, unless requested not to
if mpopt.opf.init_from_mpc ~= 1 && ...
        (strcmp(opt.alg, 'MIPS') || strcmp(opt.alg, 'IPOPT'))
    Varefs = bus(bus(:, BUS_TYPE) == REF, VA) * (pi/180);

    lb = xmin; ub = xmax;
    lb(xmin == -Inf) = -1e10;   %% replace Inf with numerical proxies
    ub(xmax ==  Inf) =  1e10;
    x0 = (lb + ub) / 2;         %% set x0 mid-way between bounds
    k = find(xmin == -Inf & xmax < Inf);    %% if only bounded above
    x0(k) = xmax(k) - 1;                    %% set just below upper bound
    k = find(xmin > -Inf & xmax == Inf);    %% if only bounded below
    x0(k) = xmin(k) + 1;                    %% set just above lower bound
    x0(vv.i1.Va:vv.iN.Va) = Varefs(1);  %% angles set to first reference angle
    if ny > 0
        ipwl = find(gencost(:, MODEL) == PW_LINEAR);
        c = gencost(sub2ind(size(gencost), ipwl, NCOST+2*gencost(ipwl, NCOST)));    %% largest y-value in CCV data
        x0(vv.i1.y:vv.iN.y) = max(c) + 0.1 * abs(max(c));
    end
end

%% set input for multiple time OPF

HHH = HH;
CCC = CC;
AAA = A;
lll = l;
uuu = u;
ll = l;
uu = u;
xminnn = xmin;
xmaxxx = xmax;
x000 = x0;

%% Profil Beban
T = 24;
%  T = 2;
%  fkali = [1 1.019];
%kasus 2 - aman
fkali = [1 1.019 1 1.0028 1.0633 1.1013 1.0812 1.158 1.22 1.248 1.253 1.197 1.229 1.278 1.265 1.258 1.247 1.295 1.349 1.331 1.293 1.205 1.145 1.095];

%kasus 3 - beban puncak malam nambah
% fkali = [1 1.019 1 1.0028 1.0633 1.1013 1.0812 1.158 1.22 1.248 1.253 1.197 1.229 1.278 1.265 1.258 1.247 1.395 1.449 1.431 1.393 1.205 1.145 1.095];

%kasus 4 - shedding 5x pada pagi hari 
% fkali = [1 1.019 5 1.0028 1.0633 1.1013 1.0812 1.158 1.22 1.248 1.253 1.197 1.229 1.278 1.265 1.258 1.247 1.295 1.349 1.331 1.293 1.205 1.145 1.095];%

%Sumatera 1 januari
% fkali =[1.14 1.10 1.06 1.03	1.07 1.11 1.09 1.02	1.01 1.01 1.01 1.01	1.02 1.01 1.00 1.01	1.06 1.16 1.41 1.40 1.36 1.29 1.19 1.13]; 



%% weighted sum
w1 = 0.95*(1/243.57);       % 0.95
w2 = 0.05*(1/0.9646);        % 0.05
nbat = 4;
BattUnit = [3 5 8 11];

for i = 1:ng
    HH(nb+i,nb+i)=w1*(HH(nb+i,nb+i));
    CC(nb+i)=w1*(CC(nb+i));
end

for j = 1:nbat
    HH(nb+BattUnit(j),nb+BattUnit(j))=(1/w1)*w2*(HH(nb+BattUnit(j),nb+BattUnit(j)));
end

HHH = HH;
CCC = CC;

%% extended matrixes for Pload objective funct.
Pcolumn = bus(:,3); 
LoadBus = find(Pcolumn);        % bus loadbus yang nggak nol
nloadbus = length (LoadBus);    % jumlah loadbus yang gak nol

% matriks A
% A;
% assignin('base','xmin',xmin);
% assignin('base','xmax',xmax);
zeromat = zeros((nb+nl),(nloadbus));  % nambahin matriks nol buat nambah obj function load shedding
for i = 1:nloadbus;                   % nilai 1 di matriks tambahan obj function load shedding
    zeromat(LoadBus(i),i) = 1;
end

A = [A zeromat];
AAA = A;

% matriks H
% HH;
zeromatH = zeros((nb+ng),nloadbus); %nambahin matriks nol 
HH = [HH zeromatH];
zeromatHH = zeros(nloadbus,(nb+ng+nloadbus));
HH = [HH;zeromatHH];
HHH = HH;

% matriks G
% CC;
addCC = -10^8*ones(nloadbus,1) ;
% addCC = zeros(nloadbus,1) ;
CC = [CC;addCC] ;
CCC = CC ;

% xmin
% xmin;
addxmin = zeros(nloadbus,1) ;
xmin = [xmin;addxmin] ;
xminnn = xmin ;

% xmax
% xmax ;
y = T-1;
xmaxxx = xmax;

for i = 1:y
    addmax = Pcolumn(LoadBus) * fkali(i)/baseMVA;
    xmax = [xmax;addmax];
    xmax = [xmax;xmaxxx];
end
y = y+1;
addmaxx = Pcolumn (LoadBus)*fkali(y)/baseMVA;
xmax = [xmax;addmaxx];

% lb
% l;
l(LoadBus) = 0;
l (1:nb)   = 0;
lll = l;

% ub
% u;
u(LoadBus) = 0;
u (1:nb) = 0;
uuu = u;

%% expand matrixes due to additional considered time
for i = 1:T-1;
    HH = blkdiag(HH,HHH);
    CC = [CC; CCC];
%     assignin('base','sebelumAAA',A);
    A = blkdiag(A,AAA);
    ll(igen) = fkali(i+1)*lll(igen);
    uu(igen) = fkali(i+1)*uuu(igen);
    l = [l;ll];
    u = [u;uu];
    xmin = [xmin;xminnn];
    xmax = [xmax;xmaxxx];
    x0 = [x0;x000];
end

%% Expand A matrixes due to ramp rate constraint
add01 = zeros(ng*(T-1), length(A(1,:)));
ramprate = gen(:, RAMP_10)/(baseMVA)*1;

for j=1:T-1
for i=1:ng
    add01((j-1)*ng + i, ((j-1)*(nb+ng+nloadbus)) + nb + i) = -1;
    add01((j-1)*ng + i, (j*(nb+ng+nloadbus))+ nb + i) = 1;
end
    l = [l; -ramprate];
    u = [u; ramprate];
end
% assignin('base','sebelumramp',A);
A = [A; add01];

%% Inject PV to System
nPV = 3;
UnitPV = [4 7 12];
cuaca1 = [0;0;0;0;0;1.44000000000000;7.20000000000000;10.0800000000000;51.8400000000000;64.8000000000000;86.4000000000000;87.8400000000000;84.9600000000000;76.3200000000000;61.9200000000000;38.8800000000000;8.64000000000000;1.44000000000000;0;0;0;0;0;0];
cuaca2 = [0;0;0;0;0;1.20000000000000;6;8.40000000000000;43.2000000000000;54;72;73.2000000000000;70.8000000000000;63.6000000000000;51.6000000000000;32.4000000000000;7.20000000000000;1.20000000000000;0;0;0;0;0;0];
cuaca3 = [0;0;0;0;0;0.960000000000000;4.80000000000000;6.72000000000000;34.5600000000000;43.2000000000000;57.6000000000000;58.5600000000000;56.6400000000000;50.8800000000000;41.2800000000000;25.9200000000000;5.76000000000000;0.960000000000000;0;0;0;0;0;0];
    
PV1 = (cuaca1)/(baseMVA*1000); % 1 hari L28, 2 hari L52 
PV2 = (cuaca2)/(baseMVA*1000); % 1 hari L28, 2 hari L52
PV3 = (cuaca3)/(baseMVA*1000); % 1 hari L28, 2 hari L52

% PV1 = (xlsread('cuaca_sunnyboy',1,'H3:H26'))/(baseMVA*1000) % 1 hari L28, 2 hari L52 
% PV2 = (xlsread('cuaca_sunnyboy',1,'G3:G26'))/(baseMVA*1000); % 1 hari L28, 2 hari L52
% PV3 = (xlsread('cuaca_sunnyboy',1,'I3:I26'))/(baseMVA*1000); % 1 hari L28, 2 hari L52

DayaPV = [PV1 PV2 PV3];

for k = 1:nPV
    for i = 1:T
         xmax((nloadbus+ng-UnitPV(k))*(i-1)+(nb+UnitPV(k))*(i))= DayaPV(i,k);
    end
end

%% Energi baterai mba Dini
add03 = zeros(T,length(A(1,:)));

%% Battery Data validasi acak 
%      C(Ah) V   Charging Rate   Discharging Rate  Initial SOC       
% Batt = [2460 72 5 5 0.2;    %Derating 82% dari tesis Mba Dini
%         1640 72 5 5 0.4;
%         2050 72 5 5 0.6;
%         1845 72 5 5 0.8;
%         1640 72 5 5 1.0];
% Batt = [2460 72 5 5 0.2;    %Test Program
%         2050 72 5 5 0.6;
%         1845 72 5 5 0.8;
%         1640 72 5 5 1.0;
%         1640 72 5 5 0.4];
%% Capacity from DE 
Batt = [0 72 5 5 0.8;    %last known best
        0 72 5 5 0.8;
        0 72 5 5 0.8;
        0 72 5 5 0.8;
%         0 72 5 5 0.8
        ];
batt = evalin('base','currentbatt');
Batt(:,1)= batt(1,1:4)';

nbat = 4;
BattUnit = [3 5 8 11];

C = zeros(nbat);
v = zeros(nbat);
soc0 = zeros(nbat);
Ah0 = zeros(nbat);
Eo = zeros(nbat);
Emin = zeros(nbat);
Emax = zeros(nbat);

for i = 1:nbat
    C(i)    = Batt(i,1);     
    v(i)    = Batt(i,2);
    soc0(i) = Batt (i,5);
    Ah0(i)  = C(i)*soc0(i);
    Eo(i)   = C(i)*soc0(i)*v(i)*1e-6;   % MWh
    Emin(i) = C(i)*0.2*v(i)*1e-6;       % MWh
    Emax(i) = C(i)*1*v(i)*1e-6;         % MWh
end

for i = 1:nbat
    for j = 1:T
         add03(T*(i-1)+j:T*i,(nb+ng+nloadbus)*(j-1)+nb+ BattUnit(i))= -1;
         lowerbat (T*(i-1)+1:T*i,:) = (Emin(i)-Eo(i))/(baseMVA);
         upperbat (T*(i-1)+1:T*i,:) = (Emax(i)-Eo(i))/(baseMVA);
    end
end
% assignin('base','lowerbatt',lowerbat);
% assignin('base','upperbatt',upperbat);
% assignin('base','sebelumbatt',A);
A = [A;add03];
% assignin('base','Afinal',A);
% assignin('base','lsblm',l);
l = [l;lowerbat];
% assignin('base','l',l);
% assignin('base','usblm',u);
u = [u;upperbat];
% assignin('base','u',u);

%% -----  run opf  -----
[x, f, info, output, lambda] = qps_matpower(HH, CC, A, l, u, xmin, xmax, x0, opt);
success = (info == 1);
info
% f

assignin('base','konvergen',info);
% size(x);
% assignin('base','All',All); 
 
%% extract data
All = [];
Pdemand = [];
Allxmin = [];
Allxmax = [];

for i = 1:T
All = [All x(((i-1)*(nb+ng+nloadbus))+1:((i-1)*(nb+ng+nloadbus))+(nb+ng+nloadbus))]; %x per jam
Pdemand = [Pdemand bus(iload, PD)*fkali(i)]; %Ploadbus sebelum di shedding, aktual mw
end

Pdemand = Pdemand*10^3   ;                   % Pdemand, aktual kW
teta = All (1:nb,:);                         % sudut tegangan perjam (p.u)
Pgen = All(nb+1:nb+ng,:)*100*1000;           % Pgenerator per jam, aktual kW
Pgen_mw = All(nb+1:nb+ng,:)*100;             % Pgenerator per jam, aktual MW
Pload = All(nb+ng+1:nb+ng+nloadbus,:)*10^5 ; % P di bus load setelah shedding per jam, aktual kW
Pshed = Pdemand-Pload;                       % P yang di shedding, aktual kW

% matriksH = full(HH);
% matriksC = CC;
% A = full(A);

%% Kapasitas saluran
frombus = branch(:,1);
tobus = branch (:,2);
xline = branch (:,4);
Ptrans = zeros(41,24);

for i = 1:T
    for j = 1:nl
        Ptrans(j,i)= (teta(frombus(j),i)-teta(tobus(j),i))/xline(j);
    end
end
Ptrans = Ptrans*100*1000; %Ptrans aktual, kw


%% Generator Cost 
a_coeff = gencost(:,5); 
b_coeff = gencost(:,6);
c_coeff = gencost(:,7);
Cost = zeros(13,24);

for i = 1:T
    for j = 1:ng
        Cost(j,i)= a_coeff(j)*Pgen_mw(j,i)^2+b_coeff(j)*Pgen_mw(j,i)+c_coeff(j);
    end
end

%% Battery Cost
%Ilangin biaya baterai dari mba Dini
for i= 1:ng
    for j = 1:nbat
        if i==BattUnit(j)
            Cost(i,:)= 0;
        else
        end
    end
end


%% Total seluruh biaya
TotalCost = sum(Cost); %Total biaya operasi per jam
% Total24 = sum(TotalCost); % Without Shedding

% With Shedding
CostShedding = sum(Pshed*0.094); %Biaya shedding Rp1467/kWh
CostSheddingsum = sum(CostShedding);
TotalCostsum = sum(TotalCost);
Total24 = TotalCostsum+CostSheddingsum; %Grand total biaya operasi

%% Baterai 
Pbat_sq = zeros(4,24);
SOC = zeros(4,24);
Ebatt = zeros(4,24);
AhBat = zeros(4,24);
UsedAh = zeros(4,24);

Eo = Eo*1e3;                 % kW
for i = 1:nbat
    Pbat (i,:) = Pgen(BattUnit(i),:);
    for j = 1:T
        AhBat(i,j)= Pbat(i,j)/v(i)*1e3;    % itu dikali 1e3, biar jadi Ah
        if j == 1
            Ebatt(i,j) = Eo(i)- Pbat(i,j);
            UsedAh(i,j) = Ah0(i)-AhBat(i,j);
            SOC(i,j) = UsedAh(i,j)/C(i);
        else
            Ebatt(i,j) = Ebatt(i,j-1) - Pbat(i,j);
            UsedAh(i,j) = UsedAh(i,j-1)-AhBat(i,j);
            SOC(i,j) = UsedAh(i,j)/C(i);
        end
    end
end
% SOC;

for i = 1:nbat
    for j = 1:T
        Pbat_sq(i,j) = Pbat(i,j)^2;
    end
end
Pbat_sq=Pbat_sq';
BatteryTransaction= sqrt(sum(Pbat_sq));
TotalBat = sum(BatteryTransaction);

excel = [Total24 TotalCostsum CostSheddingsum sum(sum(Pshed)) TotalBat BatteryTransaction];

%% MAGIC PUT UR VAR TO WORKSPACE!!!
assignin('base','Cost_Op',Total24);
assignin('base','excel',excel); 
% assignin('base','costmatrix',Cost); 
% assignin('base','BattEnergy',BatteryTransaction);
% assignin('base','Cost_24',TotalCost);
% assignin('base','Pgen',Pgen); 
% assignin('base','Pgen_mw',Pgen_mw); 
% assignin('base','Pdemand',Pdemand); 
% assignin('base','Pload',Pload); 
% assignin('base','Pshed',Pshed); 
% assignin('base','SOC',SOC); 
% save('quadprog.mat','A','matriksH','matriksC','x','xmax','xmin','l','u')
% save('final.mat','teta','Pgen','Pshed')
% xlswrite('weighted.xlsx',excel,3,'C47')

%% graphic plot - Vito -----------------------------------------------------------------------
% xlswrite('final.xlsx',Pshed,1,'B3')
% xlswrite('final.xlsx',Pgen,2,'C3')
% xlswrite('final.xlsx',Pdemand,3,'B3')
% xlswrite('final.xlsx',Pload,4,'B3')
% xlswrite('final.xlsx',Ptrans,6,'D3')
% xlswrite('final.xlsx',SOC,7,'B3')
% xlswrite('final.xlsx',Ebatt,8,'B3')
% xlswrite('final.xlsx',excel,9,'B3')
% 
% % plot controlable source 
% figure;
% plot(1:T,Pgen(1,:),'k-.',1:T,Pgen(2,:),'k-d',1:T,Pgen(6,:),'k-p',1:T,Pgen(9,:),'k-*',1:T,Pgen(10,:),'k-',1:T,Pgen(13,:),'k:');xlabel('Time(Hour)');ylabel('Power Generation (KW)')
% legend('PLN','MT1','Die1','Die2','Die3','MT2')
% xlabel('Time (hour)')
% ylabel('Power (kW)')
% axis([0 24 -5 110])
% set(gcf,'color','white')
% %  
% % % plot pv & baterai
% figure;
% plot(1:T,Pgen(4,:),'k:*',1:T,Pgen(7,:),'k:p',1:T,Pgen(12,:),'k:x',1:T,Pgen(3,:),'g-',1:T,Pgen(5,:),'b:',1:T,Pgen(8,:),'r-.',1:T,Pgen(11,:),'m-*');xlabel('Time(Hour)');ylabel('Power Generation (KW)')
% legend('PV1','PV2','PV3','ES-1','ES-2','ES-3','ES-4')
% xlabel('Time (hour)')
% ylabel('Power (kW)')
% set(gcf,'color','white')
% 
% % SOC baterai;
% figure;
% plot(1:T,SOC(1,:),'g-',1:T,SOC(2,:),'b:',1:T,SOC(3,:),'r-.',1:T,SOC(4,:),'m-*')
% legend('ES-1','ES-2','ES-3','ES-4')
% xlabel('Time (hour)')
% ylabel('SOC') 
% axis([0 24 0.0 1.1])
% set(gcf,'color','white')
% --------------------------------------------------------------------------------------------
% figure (1);
% plot(1:T,Pgen(1,:),'k-d',1:T,Pgen(2,:),'k-p',1:T,Pgen(3,:),'k-.');xlabel('Time (Hour)');ylabel('Power Generation (MW)')
% legend('P1','P2','P3')
% 
% figure (2);
% plot(1:T,Pload(1,:),'r-d',1:T,Pload(2,:),'r-p',1:T,Pload(3,:),'r-.');xlabel('Time (Hour)');ylabel('Load (MW)')
% legend('Load1','Load2','Load3')
% 
% figure;
% plot(1:T,Pshed(1,:),'b-d',1:T,Pshed(2,:),'b-p',1:T,Pshed(3,:),'b-.');xlabel('Time (Hour)');ylabel('Shedding (MW)')
% legend('Shedding1','Shedding2','Shedding3')
% axis([0 24 0 200])

%% plot baterai
% energi baterai
% figure;
% plot(1:T,Ebatt(1,:),'g-',1:T,Ebatt(2,:),'b:',1:T,Ebatt(3,:),'r-.',1:T,Ebatt(4,:),'m-*')
% legend('E1','E2','E3','E4')
% xlabel('Time (hour)')
% ylabel('Energy (wh)') 
% set(gcf,'color','white')
% 


% figure;
% plot(1:T,Pgen(3,:),'k-',1:T,Pgen(5,:),'k:',1:T,Pgen(8,:),'k-.',1:T,Pgen(11,:),'k-*')
% legend('Bat1','Bat2','Bat3','Bat4')
% xlabel('Time (hour)')
% ylabel('Power (kW)')
%% plot 30 bus
% figure;
% plot(1:T,Pgen(1,:),'k-.',1:T,Pgen(2,:),'b-',1:T,Pgen(3,:),'g-',1:T,Pgen(4,:),'k:*',1:T,Pgen(5,:),'b:',1:T,Pgen(6,:),'r-' ....
%     ,1:T,Pgen(7,:),'k:p',1:T,Pgen(8,:),'r-p',1:T,Pgen(9,:),'r-.');xlabel('Time (Hour)');ylabel('Power Generation (KW)')
% legend('GRID','Mikr1','Pbat1','PV1','Pbat2','Dies1','PV2','Pbat3','Dies2')

%%-----  calculate return values  -----
if ~any(isnan(x))
    %% update solution data
    Va = x(vv.i1.Va:vv.iN.Va);
    Pg = x(vv.i1.Pg:vv.iN.Pg);
    f = f + C0;

    %% update voltages & generator outputs
    bus(:, VM) = ones(nb, 1);
    bus(:, VA) = Va * 180/pi;
    gen(:, PG) = Pg * baseMVA;

    %% compute branch flows
    branch(:, [QF, QT]) = zeros(nl, 2);
    branch(:, PF) = (Bf * Va + Pfinj) * baseMVA;
    branch(:, PT) = -branch(:, PF);
end

%% package up results
mu_l = lambda.mu_l;
mu_u = lambda.mu_u;
muLB = lambda.lower;
muUB = lambda.upper;

%% update Lagrange multipliers (ERROR)
% il = find(branch(:, RATE_A) ~= 0 & branch(:, RATE_A) < 1e10);
% bus(:, [LAM_P, LAM_Q, MU_VMIN, MU_VMAX]) = zeros(nb, 4);
% gen(:, [MU_PMIN, MU_PMAX, MU_QMIN, MU_QMAX]) = zeros(size(gen, 1), 4);
% branch(:, [MU_SF, MU_ST]) = zeros(nl, 2);
% bus(:, LAM_P)       = (mu_u(ll.i1.Pmis:ll.iN.Pmis) - mu_l(ll.i1.Pmis:ll.iN.Pmis)) / baseMVA;
% branch(il, MU_SF)   = mu_u(ll.i1.Pf:ll.iN.Pf) / baseMVA;
% branch(il, MU_ST)   = mu_l(ll.i1.Pf:ll.iN.Pf) / baseMVA;
% gen(:, MU_PMIN)     = muLB(vv.i1.Pg:vv.iN.Pg) / baseMVA;
% gen(:, MU_PMAX)     = muUB(vv.i1.Pg:vv.iN.Pg) / baseMVA;
pimul = [
  mu_l - mu_u;
 -ones(ny>0, 1);    %% dummy entry corresponding to linear cost row in A (in MINOS)
  muLB - muUB
];

mu = struct( ...
  'var', struct('l', muLB, 'u', muUB), ...
  'lin', struct('l', mu_l, 'u', mu_u) );

results = mpc;
[results.bus, results.branch, results.gen, ...
    results.om, results.x, results.mu, results.f] = ...
        deal(bus, branch, gen, om, x, mu, f);
% toc
raw = struct('xr', x, 'pimul', pimul, 'info', info, 'output', output);
