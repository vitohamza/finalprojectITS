function mpc = case30
%CASE30    Power flow data for 30 bus, 6 generator case.
%   Please see CASEFORMAT for details on the case file format.
%
%   Based on data from ...
%     Alsac, O. & Stott, B., "Optimal Load Flow with Steady State Security",
%     IEEE Transactions on Power Apparatus and Systems, Vol. PAS 93, No. 3,
%     1974, pp. 745-751.
%   ... with branch parameters rounded to nearest 0.01, shunt values divided
%   by 100 and shunt on bus 10 moved to bus 5, load at bus 5 zeroed out.
%   Generator locations, costs and limits and bus areas were taken from ...
%     Ferrero, R.W., Shahidehpour, S.M., Ramesh, V.C., "Transaction analysis
%     in deregulated power systems using game theory", IEEE Transactions on
%     Power Systems, Vol. 12, No. 3, Aug 1997, pp. 1340-1347.
%   Generator Q limits were derived from Alsac & Stott, using their Pmax
%   capacities. V limits and line |S| limits taken from Alsac & Stott.

%   MATPOWER

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	3	0         0       0	0       1	1	0	135	1	1.05	0.95;
	2	2	0.01205   12.7	  0	0       1	1	0	135	1	1.1     0.95;
	3	1	0.001371  1.2     0	0       1	1	0	135	1	1.05	0.95;
	4	1	0.004343  1.6     0	0       1	1	0	135	1	1.05	0.95;
	5	1	0         0       0	0.19	1	1	0	135	1	1.05	0.95;
	6	1	0         0       0	0       1	1	0	135	1	1.05	0.95;
	7	1	0.013028  10.9	  0	0       1	1	0	135	1	1.05	0.95;
	8	1	0.017143  30      0	0       1	1	0	135	1	1.05	0.95;
	9	1	0         0       0	0       1	1	0	135	1	1.05	0.95;
	10	1	0.003314  2       0	0       3	1	0	135	1	1.05	0.95;
	11	1	0         0       0	0       1	1	0	135	1	1.05	0.95;
	12	1	0.0064	  7.5     0	0       2	1	0	135	1	1.05	0.95;
	13	2	0         0       0	0       2	1	0	135	1	1.1     0.95;
	14	1	0.003543  1.6     0	0       2	1	0	135	1	1.05	0.95;
	15	1	0.004685  2.5     0	0       2	1	0	135	1	1.05	0.95;
	16	1	0.002     1.8     0	0       2	1	0	135	1	1.05	0.95;
	17	1	0.005143  5.8     0	0       2	1	0	135	1	1.05	0.95;
	18	1	0.001828  0.9     0	0       2	1	0	135	1	1.05	0.95;
	19	1	0.005429  3.4     0	0       2	1	0	135	1	1.05	0.95;
	20	1	0.001257  0.7 	  0	0       2	1	0	135	1	1.05	0.95;
	21	1	0.01	  11.2	  0	0       3	1	0	135	1	1.05	0.95;
	22	2	0         0       0	0       3	1	0	135	1	1.1     0.95;
	23	2	0.001828  1.6     0	0       2	1	0	135	1	1.1     0.95;
	24	1	0.00497   6.7     0	0.04	3	1	0	135	1	1.05	0.95;
	25	1	0         0       0	0       3	1	0	135	1	1.05	0.95;
	26	1	0.002     2.3     0	0       3	1	0	135	1	1.05	0.95;
	27	2	0         0       0	0       3	1	0	135	1	1.1     0.95;
	28	1	0         0       0	0       1	1	0	135	1	1.05	0.95;
	29	1	0.00137   0.9     0	0       3	1	0	135	1	1.05	0.95;
	30	1	0.006057  1.9     0	0       3	1	0	135	1	1.05	0.95;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
% mpc.gen = [
% 	1	23.54	0	150 	-20	1	100	1	80	0	0	0	0	0	0	0	0	80	0	0	0;  %--80--
% 	2	60.97	0	60      -20	1	100	1	80	0	0	0	0	0	0	0	0	80	0	0	0;  %--80--
% 	22	21.59	0	62.5	-15	1	100	1	50	0	0	0	0	0	0	0	0	80	0	0	0;  %--50-- bus 22
% 	27	26.91	0	48.7	-15	1	100	1	55	0	0	0	0	0	0	0	0	80	0	0	0;  %--55-- bus 27
% 	23	19.2	0	40      -10	1	100	1	30	0	0	0	0	0	0	0	0	80	0	0	0;  %--30-- bus 23
% 	13	37      0	44.7	-15	1	100	1	40	0	0	0	0	0	0	0	0	80	0	0	0;  %--40-- bus 13
% ];

mpc.gen = [
	1	0.02354	0	150 	-20	1	100	1	0.080	0       0	0	0	0	0	0	0	0.015	0	0	0;  % PLN
	2	0.06097	0	60      -20	1	100	1	0.045	0       0	0	0	0	0	0	0	1.050  	0	0	0;  % PV 1    
	13	0.037   0	44.7	-15	1	100	1	0.055   0       0	0	0	0	0	0	0	0.015	0	0	0;  % Diesel   
	22	0.02159	0	62.5	-15	1	100	1	0.0126	-0.0126	0	0	0	0	0	0	0	0.0252	0	0	0;  % Baterai
	23	0.0192	0	40      -10	1	100	1	0.0126	-0.0126	0	0	0	0	0	0	0	0.0252	0	0	0;  % Baterai    
	27	0.02691	0	48.7	-15	1	100	1	0.0550	0       0	0	0	0	0	0	0	1.050	0	0	0;  % PV 2
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1	2	0.02	0.06	0.03	0.130	0.130	0.130	0	0	1	-360	360;
	1	3	0.05	0.19	0.02	0.130	0.130	0.130	0	0	1	-360	360;
	2	4	0.06	0.17	0.02	0.065	0.065	0.065	0	0	1	-360	360;
	3	4	0.01	0.04	0       0.130	0.130	0.130	0	0	1	-360	360;
	2	5	0.05	0.2     0.02	0.130	0.130	0.130	0	0	1	-360	360;
	2	6	0.06	0.18	0.02	0.065	0.065	0.065	0	0	1	-360	360;
	4	6	0.01	0.04	0       0.090	0.090	0.090	0	0	1	-360	360;
	5	7	0.05	0.12	0.01	0.070	0.070	0.070	0	0	1	-360	360;
	6	7	0.03	0.08	0.01	0.130	0.130	0.130	0	0	1	-360	360;
	6	8	0.01	0.04	0       0.032	0.032	0.032	0	0	1	-360	360;
	6	9	0       0.21	0       0.065	0.065	0.065	0	0	1	-360	360;
	6	10	0       0.56	0       0.032	0.032	0.032	0	0	1	-360	360;
	9	11	0       0.21	0       0.065	0.065	0.065	0	0	1	-360	360;
	9	10	0       0.11	0       0.065	0.065	0.065	0	0	1	-360	360;
	4	12	0       0.26	0       0.065	0.065	0.065	0	0	1	-360	360;
	12	13	0       0.14	0       0.065	0.065	0.065	0	0	1	-360	360;
	12	14	0.12	0.26	0       0.032	0.032	0.032	0	0	1	-360	360;
	12	15	0.07	0.13	0       0.032	0.032	0.032	0	0	1	-360	360;
	12	16	0.09	0.2     0       0.032	0.032	0.032	0	0	1	-360	360;
	14	15	0.22	0.2     0       0.016	0.016	0.016	0	0	1	-360	360;
	16	17	0.08	0.19	0       0.016	0.016	0.016	0	0	1	-360	360;
	15	18	0.11	0.22	0       0.016	0.016	0.016	0	0	1	-360	360;
	18	19	0.06	0.13	0       0.016	0.016	0.016	0	0	1	-360	360;
	19	20	0.03	0.07	0       0.032	0.032	0.032	0	0	1	-360	360;
	10	20	0.09	0.21	0       0.032	0.032	0.032	0	0	1	-360	360;
	10	17	0.03	0.08	0       0.032	0.032	0.032	0	0	1	-360	360;
	10	21	0.03	0.07	0       0.032	0.032	0.032	0	0	1	-360	360;
	10	22	0.07	0.15	0       0.032	0.032	0.032	0	0	1	-360	360;
	21	22	0.01	0.02	0       0.032	0.032	0.032	0	0	1	-360	360;
	15	23	0.1     0.2     0       0.016	0.016	0.016	0	0	1	-360	360;
	22	24	0.12	0.18	0       0.016	0.016	0.016	0	0	1	-360	360;
	23	24	0.13	0.27	0       0.016	0.016	0.016	0	0	1	-360	360;
	24	25	0.19	0.33	0       0.016	0.016	0.016	0	0	1	-360	360;
	25	26	0.25	0.38	0       0.016	0.016	0.016	0	0	1	-360	360;
	25	27	0.11	0.21	0       0.016	0.016	0.016	0	0	1	-360	360;
	28	27	0       0.4     0       0.065	0.065	0.065	0	0	1	-360	360;
	27	29	0.22	0.42	0       0.016	0.016	0.016	0	0	1	-360	360;
	27	30	0.32	0.6     0       0.016	0.016	0.016	0	0	1	-360	360;
	29	30	0.24	0.45	0       0.016	0.016	0.016	0	0	1	-360	360;
	8	28	0.06	0.2     0.02	0.032	0.032	0.032	0	0	1	-360	360;
	6	28	0.02	0.06	0.01	0.032	0.032	0.032	0	0	1	-360	360;
];

%%-----  OPF Data  -----%%
%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
% mpc.gencost = [
% 	2	0	0	3	0.02	2       0; % gen 1
% 	2	0	0	3	0.0175	1.75	0; % gen 2
% 	2	0	0	3	0.0625	1       0; % gen 22
% 	2	0	0	3	0.00834	3.25	0; % gen 27
% 	2	0	0	3	0.025	3       0; % gen 23
% 	2	0	0	3	0.025	3       0; % gen 13
% ];

mpc.gencost = [
	2	0	0	3	0.1625	1       0; % PLN
	2	0	0	3	0 		1       0; % PV 1 
	2	0	0	3	0.0175  1.75    0; % diesel 
	2	0	0	3	1	    1     0; % Baterai 
	2	0	0	3	1   	1     0; % Baterai    
	2	0	0	3	0   	1	    0; % PV2
];