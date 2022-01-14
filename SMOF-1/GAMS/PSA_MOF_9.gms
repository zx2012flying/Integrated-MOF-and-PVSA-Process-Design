Set

comp                     'propane PA: pass, propene PE: adsorb'                          /PA, PE/

t                        'time discretization'                                           /0 * 60/
t_start(t)               'starting point'                                                /0/
t_pre(t)                 'time discretization for pressurization'                        /1 * 10/
t_ads(t)                 'time discretization for adsorption'                            /11 * 30/
t_rin(t)                 'time discretization for rinse'                                 /31 * 40/
t_des(t)                 'time discretization for desorption'                            /41 * 60/

z                        'column axis discretization'                                    /0 * 21/
z_start(z)               'axis starting point'                                           /1/
z_mid(z)                 'axis middle section'                                           /2 * 19/
z_end(z)                 'axis ending point'                                             /20/
z_left(z)                'left ghost'                                                    /0/
z_axis(z)                'axis point'                                                    /1 * 20/
z_right(z)               'right ghost'                                                   /21/
;

Parameters

Qsat_1_PA                'isotherm parameter: mol/kg'                                    /0.0572/
Qsat_2_PA                'isotherm parameter: mol/kg'                                    /1.0855/
b_1_PA                   'isotherm parameter'                                            /2.3024/
b_2_PA                   'isotherm parameter'                                            /0.1/

*Qsat_1_PE                'isotherm parameter: mol/kg'                                    /0.6694/
*Qsat_2_PE                'isotherm parameter: mol/kg'                                    /1.0256/
*b_1_PE                   'isotherm parameter'                                            /4.2127/
*b_2_PE                   'isotherm parameter'                                            /0.5714/
Qsat_1_PE                'isotherm parameter: mol/kg'                                    /0.3997/
Qsat_2_PE                'isotherm parameter: mol/kg'                                    /1.1623/
b_1_PE                   'isotherm parameter'                                            /1.4786/
b_2_PE                   'isotherm parameter'                                            /1.4786/
;

Parameters

u_ini(t, z)              'velocity initial value:                m/s'
P_ini(t, z)              'pressure initial value:                Pa'
yE_ini(t, z)             'PE composition initial value:          /'
yA_ini(t, z)             'PE composition initial value:          /'
;

$CALL 'del PSA_MOF_9.gdx'
$CALL GDXXRW I=PSA_MOF_9.xlsx O=PSA_MOF_9.gdx Index=In_index!A1

$GDXIN PSA_MOF_9.gdx
$LOAD u_ini, P_ini, yE_ini, yA_ini

Parameters

yE_feed                  'input PE composition'                                          /0.85/
yA_feed                  'input PA composition'
yE_rinse                 'rinse input PE composition'                                    /0.99/
yA_rinse                 'rinse input PA composition'
P_input                  'input pressure from FCC unit           Pa'                     /202650/
P_atm                    'atmosphere pressure                    Pa'                     /101325/

MW_PA                    'molecular weight of PA:                kg/mol'                 /0.0441/
MW_PE                    'molecular weight of PA:                kg/mol'                 /0.04208/
kA                       'mass transfer coefficient PA:          1/s'                    /0.8307/
kE                       'mass transfer coefficient PE:          1/s'                    /0.6118/

Temp                     'isothermal PSA temperature (Dai Tang): K'                      /300/
Length                   'length of column pilot plant           m'                      /2/
Diameter                 'diameter of column pilot plant         m'                      /0.4/
Porosity_bed             'porosity of column bed                 /'                      /0.43/
rho_s                    'particle density of solid adsorbent:   kg/m3'                  /703/
viscosity                'viscosity of gas:                      kg/m/s'                 /8e-6/
Radius_SP                'solid particle radius                  m'                      /5e-3/

pi                       'pi'                                                            /3.1416/
R                        'gas constant                           kg*m2/k/mol/s2'         /8.3145/
Eff                      '1 / efficiency of vacuum and compressor'                       /1.333/
IC                       'isentropic coeff for PE k = 1.15,      k/k-1'                  /7.67/

Nubz                     'number of z interval+1'
Dz                       'delta z: m'
DzS                      'square of delta z'
t_pre_des                'total duration                                 60'

CSA                      'cross section area                     m2      '
EM                       'parameter Ergun equation               kg/m3/s '
FM                       'parameter Ergun equation               1/m     '
M                        'parameter in mass balance                      '
ERTI                     'Eff * R * Temp * IC                             0.026'
CPRT                     'CSA * Porosity_bed / (R * Temp)'
;

yA_feed = 1 - yE_feed;
yA_rinse = 1 - yE_rinse;

Nubz = card(z_axis)+1;
Dz = Length / card(z_axis);
DzS = Dz * Dz;
t_pre_des = card(t_pre) + card(t_ads) + card(t_rin) + card(t_des);

CSA = pi * Diameter * Diameter / 4;
EM = 150 * viscosity * (1 - Porosity_bed) * (1 - Porosity_bed) / (4 * Radius_SP * Radius_SP * Porosity_bed * Porosity_bed * Porosity_bed);
FM = 1.75 * (1 - Porosity_bed) / (2 * Radius_SP * Porosity_bed * Porosity_bed * Porosity_bed);
M = (1 - Porosity_bed) * R * rho_s / Porosity_bed;
ERTI = Eff * IC * R * Temp / 1e6;
CPRT = CSA * Porosity_bed / (R * Temp);

Positive variables

Dt_pre                   'time interval at pressurization                second'
Dt_ads                   'time interval at adsorption                    second'
Dt_rin                   'time interval at resin                         second'
Dt_des                   'time interval at desorption                    second'
P_high                   'high pressure for adsorption                   Pa'
P_low                    'low pressure for desorption                    Pa'
P_end                    'adsorption rinse end pressure                  Pa'

yA(t, z)                 'gas molar fraction of component PA             /'
qA(t, z)                 'solid loading of PA                            mol/kg'
qA_Star(t, z)            'equilibrium loading of PA                      mol/kg'
yE(t, z)                 'gas molar fraction of component PE             /'
qE(t, z)                 'solid loading of PE                            mol/kg'
qE_Star(t, z)            'equilibrium loading of PE                      mol/kg'
P(t, z)                  'pressure in the column                         Pa'

PA_pre_feed              'PA feed at pre: / R / Temp                     mol'
PA_ads_feed              'PA feed at ads: / R / Temp                     mol'
PA_ads_out               'PA output at pre: / R / Temp                   mol'
PA_rin_feed              'PA feed at rin: / R / Temp                     mol'
PA_rin_out               'PA output at rin: / R / Temp                   mol'
PA_des_exhaust           'PA exhaust at des: / R / Temp                  mol'

PE_pre_feed              'PE feed at pre: / R / Temp                     mol'
PE_ads_feed              'PE feed at ads: / R / Temp                     mol'
PE_ads_out               'PE output at pre: / R / Temp                   mol'
PE_rin_feed              'PE feed at rin: / R / Temp                     mol'
PE_rin_out               'PE output at rin: / R / Temp                   mol'
PE_des_exhaust           'PE exhaust at des: / R / Temp                  mol'
;

Variables

PurityA                  'purity of PA                                   /'
PurityE                  'purity of PE                                   /'
RecoveryA                'recovery of PA                                 /'
RecoveryE                'recovery of PE                                 /'
PE_product               '1% PA + 99% PE product                         mol'
Mass_balanceA            'overall mass balances of PA: input - output    mol'
Mass_balanceE            'overall mass balances of PE: input - output    mol'

W_12_feed                'vacuum pump workload at pre+ads feed           MJ'
W_23_out                 'vacuum pump workload at ads+rin out            MJ'
W_3_feed                 'vacuum pump workload at rin feed               MJ'
W_4_out                  'vacuum pump workload at des out                MJ'

u(t, z)                  'gas velocity in the column                     m/s'
OBJ                      'objective function                             J/mol';

P_high.lo = P_atm;
P_high.up = P_input;
P_high.l = 2e5;

P_end.lo = P_atm;
P_end.up = P_input;
P_end.l = 1.8e5;

P_low.lo = 1e3;
P_low.up = P_atm;
P_low.l = 1e4;

P.lo(t, z) = 1e3;
*P.up(t, z) = P_input;
P.l(t, z) = P_ini(t, z);

Dt_pre.lo = 0.1;
Dt_pre.up = 3;
Dt_pre.l = 2;

Dt_ads.lo = 0.1;
Dt_ads.up = 10;
Dt_ads.l = 3;

Dt_rin.lo = 0.1;
Dt_rin.up = 20;
Dt_rin.l = 5;

Dt_des.lo = 0.1;
Dt_des.up = 30;
Dt_des.l = 30;

yA.lo(t, z) = 0;
yA.up(t, z) = 1;
yA.l(t, z) = yA_ini(t, z);

yE.lo(t, z) = 0;
yE.up(t, z) = 1;
yE.l(t, z) = yE_ini(t, z);

qA.lo(t, z) = 0;
qA.up(t, z) = 6;
qA.l(t, z) = 0.2;

qE.lo(t, z) = 0;
qE.up(t, z) = 6;
qE.l(t, z) = 3.8;

qA_Star.lo(t, z) = 0;
qA_Star.up(t, z) = 6;
qA_Star.l(t, z) = 0.6;

qE_Star.lo(t, z) = 0;
qE_Star.up(t, z) = 6;
qE_Star.l(t, z) = 8;

u.lo(t, z) = -40;
u.up(t, z) = 40;
u.l(t, z) = u_ini(t, z);

PA_pre_feed.lo = 1e-2;
PA_pre_feed.up = 2000;
PA_pre_feed.l = 10;

PA_ads_feed.lo = 1e-2;
PA_ads_feed.up = 2000;
PA_ads_feed.l = 20;

PA_ads_out.lo = 1e-2;
PA_ads_out.up = 2000;
PA_ads_out.l = 15;

PA_rin_feed.lo = 1e-3;
PA_rin_feed.up = 100;
PA_rin_feed.l = 1;

PA_rin_out.lo = 1e-2;
PA_rin_out.up = 2000;
PA_rin_out.l = 8;

PA_des_exhaust.lo = 1e-3;
PA_des_exhaust.up = 1000;
PA_des_exhaust.l = 1;

PE_pre_feed.lo = 1e-2;
PE_pre_feed.up = 5000;
PE_pre_feed.l = 50;

PE_ads_feed.lo = 1e-2;
PE_ads_feed.up = 5000;
PE_ads_feed.l = 100;

PE_ads_out.lo = 1e-2;
PE_ads_out.up = 5000;
PE_ads_out.l = 100;

PE_rin_feed.lo = 1e-2;
PE_rin_feed.up = 5000;
PE_rin_feed.l = 50;

PE_rin_out.lo = 1e-2;
PE_rin_out.up = 5000;
PE_rin_out.l = 50;

PE_des_exhaust.lo = 1e-2;
PE_des_exhaust.up = 5000;
PE_des_exhaust.l = 60;

PE_product.lo = 1e-6;
PE_product.l = 5;

Mass_balanceA.l = 0;
Mass_balanceE.l = 0;

W_12_feed.l = 0;
W_23_out.l = 0;
W_3_feed.l = 0.8;
W_4_out.l = 0;

OBJ.l = 1;

Equations

LGEq00, LGEq01, LGEq02, LGEq03, LGEq04
ISEq00, ISEq01
PREq00, PREq01, PREq02, PREq03, PREq04, PREq05, PREq06, PREq07, PREq08, PREq09, PREq10
ADEq00, ADEq01, ADEq02, ADEq03, ADEq04, ADEq05, ADEq06, ADEq07, ADEq08, ADEq09, ADEq10
REEq00, REEq01, REEq02, REEq03, REEq04, REEq05, REEq06, REEq07, REEq08, REEq09, REEq10
DEEq00, DEEq01, DEEq02, DEEq03, DEEq04, DEEq05, DEEq06, DEEq07, DEEq08, DEEq09, DEEq10
CSEq00, CSEq01, CSEq02, CSEq03
MBEq00, MBEq01, MBEq02, MBEq03, MBEq04, MBEq05, MBEq06, MBEq07, MBEq08, MBEq09, MBEq10, MBEq11
SSEq00, SSEq01, SSEq02, SSEq03, SSEq04, SSEq05, SSEq06, SSEq07, SSEq08
ECEq00, ECEq01, ECEq02, ECEq03, ECEq04
;
*****************************************
*          logical constraint           *
*****************************************
LGEq00(t, z)$(t_ads(t) and z_axis(z))    ..      u(t, z)                          =g= 0                                          ;
LGEq01(t, z)$(t_rin(t) and z_axis(z))    ..      u(t, z)                          =g= 0                                          ;
LGEq02(t, z)                             ..      yA(t, z)                         =e= 1 - yE(t, z)                               ;
LGEq03                                   ..      P_high                           =g= 1.05 * P_end                               ;
LGEq04                                   ..      P_end                            =g= 1.05 * P_low                               ;
*****************************************
*      Dual-site Langmuir isotherm      *
*****************************************
ISEq00(t, z)$z_axis(z)                   ..      qA_Star(t, z)                    =e= Qsat_1_PA * b_1_PA * P(t, z) * yA(t, z) / (1e5 + b_1_PA * P(t, z) * yA(t, z) + b_1_PE * P(t, z) * yE(t, z)) + Qsat_2_PA * b_2_PA * P(t, z) * yA(t, z) / (1e5 + b_2_PA * P(t, z) * yA(t, z) + b_2_PE * P(t, z) * yE(t, z));
ISEq01(t, z)$z_axis(z)                   ..      qE_Star(t, z)                    =e= Qsat_1_PE * b_1_PE * P(t, z) * yE(t, z) / (1e5 + b_1_PA * P(t, z) * yA(t, z) + b_1_PE * P(t, z) * yE(t, z)) + Qsat_2_PE * b_2_PE * P(t, z) * yE(t, z) / (1e5 + b_2_PA * P(t, z) * yA(t, z) + b_2_PE * P(t, z) * yE(t, z));
*****************************************
*         Pressurization step           *
*****************************************
PREq00(t, z)$(t_pre(t) and z_left(z))    ..      yE(t, z)                         =e= yE_feed                                    ;
PREq01(t, z)$(t_pre(t) and z_left(z))    ..      P(t, z)                          =e= P_high                                     ;
PREq02(t, z)$(t_pre(t) and z_right(z))   ..      yE(t, z)                         =e= yE(t, z - 1)                               ;
PREq03(t, z)$(t_pre(t) and z_right(z))   ..      P(t, z)                          =e= P(t, z - 1)                                ;
PREq04(t, z)$(t_pre(t) and z_right(z))   ..      u(t, z)                          =e= u(t, z - 1)                                ;
PREq05(t, z)$(t_pre(t) and z_right(z))   ..      u(t, z)                          =e= 0                                          ;

PREq06(t, z)$(t_pre(t) and z_axis(z))    ..      qA(t, z) - qA(t-1, z)            =e= Dt_pre * kA * (qA_Star(t, z) - qA(t, z))   ;
PREq07(t, z)$(t_pre(t) and z_axis(z))    ..      qE(t, z) - qE(t-1, z)            =e= Dt_pre * kE * (qE_Star(t, z) - qE(t, z))   ;

PREq08(t, z)$(t_pre(t) and z_axis(z))    ..      (yE(t, z) - yE(t-1, z)) / Dt_pre =e= - u(t, z) * (yE(t, z) - yE(t, z-1)) / Dz - M * Temp * ( (qE(t, z) - qE(t-1, z)) - yE(t, z) * (qE(t, z) - qE(t-1, z) + qA(t, z) - qA(t-1, z)) ) / (P(t, z) * Dt_pre)        ;
PREq09(t, z)$(t_pre(t) and z_axis(z))    ..      (P(t, z) - P(t-1, z)) / Dt_pre   =e= - (P(t, z) * u(t, z) - P(t, z-1) * u(t, z-1)) / Dz - M * Temp * ( qE(t,z) - qE(t-1, z) + qA(t,z) - qA(t-1, z) ) / Dt_pre                                                   ;
PREq10(t, z)$(t_pre(t) and z_axis(z))    ..      (P(t, z) - P(t, z-1)) / Dz       =e= - EM * u(t, z-1) - FM * u(t, z-1) * ((u(t, z-1)*u(t, z-1))**0.5) * P(t, z-1) * (MW_PA * yA(t, z-1) + MW_PE * yE(t, z-1)) / (R * Temp)                                      ;
*****************************************
*           Adsorption step             *
*****************************************
ADEq00(t, z)$(t_ads(t) and z_left(z))    ..      yE(t, z)                         =e= yE_feed                                    ;
ADEq01(t, z)$(t_ads(t) and z_left(z))    ..      P(t, z)                          =e= P_high                                     ;
ADEq02(t, z)$(t_ads(t) and z_right(z))   ..      yE(t, z)                         =e= yE(t, z - 1)                               ;
ADEq03(t, z)$(t_ads(t) and z_right(z))   ..      P(t, z)                          =e= P_end                                      ;
ADEq04(t, z)$(t_ads(t) and z_right(z))   ..      u(t, z)                          =e= u(t, z - 1)                                ;

ADEq05(t, z)$(t_ads(t) and z_axis(z))    ..      qA(t, z) - qA(t-1, z)            =e= Dt_ads * kA * (qA_Star(t, z) - qA(t, z))   ;
ADEq06(t, z)$(t_ads(t) and z_axis(z))    ..      qE(t, z) - qE(t-1, z)            =e= Dt_ads * kE * (qE_Star(t, z) - qE(t, z))   ;

ADEq07(t, z)$(t_ads(t) and z_axis(z))    ..      (yE(t, z) - yE(t-1, z)) / Dt_ads =e= - u(t, z) * (yE(t, z) - yE(t, z-1)) / Dz - M * Temp * ( (qE(t, z) - qE(t-1, z)) - yE(t, z) * (qE(t, z) - qE(t-1, z) + qA(t, z) - qA(t-1, z)) ) / (P(t, z) * Dt_ads)        ;
ADEq08(t, z)$(t_ads(t) and z_axis(z))    ..      (P(t, z) - P(t-1, z)) / Dt_ads   =e= - (P(t, z) * u(t, z) - P(t, z-1) * u(t, z-1)) / Dz - M * Temp * ( qE(t,z) - qE(t-1, z) + qA(t,z) - qA(t-1, z) ) / Dt_ads                                                   ;
ADEq09(t, z)$(t_ads(t) and z_axis(z))    ..      (P(t, z) - P(t, z-1)) / Dz       =e= - EM * u(t, z-1) - FM * u(t, z-1) * ((u(t, z-1)*u(t, z-1))**0.5) * P(t, z-1) * (MW_PA * yA(t, z-1) + MW_PE * yE(t, z-1)) / (R * Temp)                                      ;
ADEq10(t, z)$(t_ads(t) and z_right(z))   ..      (P(t, z) - P(t, z-1)) / Dz       =e= - EM * u(t, z-1) - FM * u(t, z-1) * ((u(t, z-1)*u(t, z-1))**0.5) * P(t, z-1) * (MW_PA * yA(t, z-1) + MW_PE * yE(t, z-1)) / (R * Temp)                                      ;
*****************************************
*              resin step               *
*****************************************
REEq00(t, z)$(t_rin(t) and z_left(z))    ..      yE(t, z)                         =e= yE_rinse                                   ;
REEq01(t, z)$(t_rin(t) and z_left(z))    ..      P(t, z)                          =e= P_high                                     ;
REEq02(t, z)$(t_rin(t) and z_right(z))   ..      yE(t, z)                         =e= yE(t, z - 1)                               ;
REEq03(t, z)$(t_rin(t) and z_right(z))   ..      P(t, z)                          =e= P_end                                      ;
REEq04(t, z)$(t_rin(t) and z_right(z))   ..      u(t, z)                          =e= u(t, z - 1)                                ;

REEq05(t, z)$(t_rin(t) and z_axis(z))    ..      qA(t, z) - qA(t-1, z)            =e= Dt_rin * kA * (qA_Star(t, z) - qA(t, z))   ;
REEq06(t, z)$(t_rin(t) and z_axis(z))    ..      qE(t, z) - qE(t-1, z)            =e= Dt_rin * kE * (qE_Star(t, z) - qE(t, z))   ;

REEq07(t, z)$(t_rin(t) and z_axis(z))    ..      (yE(t, z) - yE(t-1, z)) / Dt_rin =e= - u(t, z) * (yE(t, z) - yE(t, z-1)) / Dz - M * Temp * ( (qE(t, z) - qE(t-1, z)) - yE(t, z) * (qE(t, z) - qE(t-1, z) + qA(t, z) - qA(t-1, z)) ) / (P(t, z) * Dt_rin)        ;
REEq08(t, z)$(t_rin(t) and z_axis(z))    ..      (P(t, z) - P(t-1, z)) / Dt_rin   =e= - (P(t, z) * u(t, z) - P(t, z-1) * u(t, z-1)) / Dz - M * Temp * ( qE(t,z) - qE(t-1, z) + qA(t,z) - qA(t-1, z) ) / Dt_rin                                                   ;
REEq09(t, z)$(t_rin(t) and z_axis(z))    ..      (P(t, z) - P(t, z-1)) / Dz       =e= - EM * u(t, z-1) - FM * u(t, z-1) * ((u(t, z-1)*u(t, z-1))**0.5) * P(t, z-1) * (MW_PA * yA(t, z-1) + MW_PE * yE(t, z-1)) / (R * Temp)                                      ;
REEq10(t, z)$(t_rin(t) and z_right(z))   ..      (P(t, z) - P(t, z-1)) / Dz       =e= - EM * u(t, z-1) - FM * u(t, z-1) * ((u(t, z-1)*u(t, z-1))**0.5) * P(t, z-1) * (MW_PA * yA(t, z-1) + MW_PE * yE(t, z-1)) / (R * Temp)                                      ;
*****************************************
*           Desorption step             *
*****************************************
DEEq00(t, z)$(t_des(t) and z_left(z))    ..      yE(t, z)                         =e= yE(t, z + 1)                               ;
DEEq01(t, z)$(t_des(t) and z_left(z))    ..      P(t, z)                          =e= P_low                                      ;
DEEq02(t, z)$(t_des(t) and z_left(z))    ..      u(t, z)                          =e= u(t, z+1)                                  ;
DEEq03(t, z)$(t_des(t) and z_right(z))   ..      yE(t, z)                         =e= yE(t, z - 1)                               ;
DEEq04(t, z)$(t_des(t) and z_right(z))   ..      P(t, z)                          =e= P(t, z - 1)                                ;
DEEq05(t, z)$(t_des(t) and z_right(z))   ..      u(t, z)                          =e= 0                                          ;

DEEq06(t, z)$(t_des(t) and z_axis(z))    ..      qA(t, z) - qA(t-1, z)            =e= Dt_des * kA * (qA_Star(t, z) - qA(t, z))   ;
DEEq07(t, z)$(t_des(t) and z_axis(z))    ..      qE(t, z) - qE(t-1, z)            =e= Dt_des * kE * (qE_Star(t, z) - qE(t, z))   ;

DEEq08(t, z)$(t_des(t) and z_axis(z))    ..      (yE(t, z) - yE(t-1, z)) / Dt_des =e= - u(t, z) * (yE(t, z+1) - yE(t, z)) / Dz - M * Temp * ( (qE(t, z) - qE(t-1, z)) - yE(t, z) * (qE(t, z) - qE(t-1, z) + qA(t, z) - qA(t-1, z)) ) / (P(t, z) * Dt_des)        ;
DEEq09(t, z)$(t_des(t) and z_axis(z))    ..      (P(t, z) - P(t-1, z)) / Dt_des   =e= - (P(t, z+1) * u(t, z+1) - P(t, z) * u(t, z)) / Dz - M * Temp * ( qE(t,z) - qE(t-1, z) + qA(t,z) - qA(t-1, z) ) / Dt_des                                                   ;
DEEq10(t, z)$(t_des(t) and z_axis(z))    ..      (P(t, z) - P(t, z-1)) / Dz       =e= - EM * u(t, z) - FM * u(t, z) * ((u(t, z)*u(t, z))**0.5) * P(t, z) * (MW_PA * yA(t, z) + MW_PE * yE(t, z)) / (R * Temp)                                                    ;
*****************************************
*       cyclic steady state: CSS        *
*****************************************
CSEq00(t, z)$(t_start(t) and z_axis(z))  ..      yE(t, z) - yE(t+t_pre_des, z)    =e= 0                                          ;
CSEq01(t, z)$(t_start(t) and z_axis(z))  ..      P(t, z)  - P(t+t_pre_des, z)     =e= 0                                          ;
CSEq02(t, z)$(t_start(t) and z_axis(z))  ..      qE(t, z) - qE(t+t_pre_des, z)    =e= 0                                          ;
CSEq03(t, z)$(t_start(t) and z_axis(z))  ..      qA(t, z) - qA(t+t_pre_des, z)    =e= 0                                          ;
*****************************************
*              mass balances            *
*****************************************
MBEq00(z)$z_left(z)                      ..      PA_pre_feed     =e= P_high * sum(t_pre, u(t_pre, z)) * yA_feed * Dt_pre * CPRT                                          ;
MBEq01(z)$z_left(z)                      ..      PA_ads_feed     =e= P_high * sum(t_ads, u(t_ads, z)) * yA_feed * Dt_ads * CPRT                                          ;
MBEq02(z)$z_end(z)                       ..      PA_ads_out      =e= sum(t_ads, P(t_ads, z) * u(t_ads, z) * yA(t_ads, z)) * Dt_ads * CPRT                                ;
MBEq03(z)$z_left(z)                      ..      PA_rin_feed     =e= P_high * sum(t_rin, u(t_rin, z)) * yA_rinse * Dt_rin * CPRT                                         ;
MBEq04(z)$z_end(z)                       ..      PA_rin_out      =e= sum(t_rin, P(t_rin, z) * u(t_rin, z) * yA(t_rin, z)) * Dt_rin * CPRT                                ;
MBEq05(z)$z_start(z)                     ..      PA_des_exhaust  =e= - sum(t_des, P(t_des, z) * u(t_des, z) * yA(t_des, z)) * Dt_des * CPRT                              ;
MBEq06(z)$z_left(z)                      ..      PE_pre_feed     =e= P_high * sum(t_pre, u(t_pre, z)) * yE_feed * Dt_pre * CPRT                                          ;
MBEq07(z)$z_left(z)                      ..      PE_ads_feed     =e= P_high * sum(t_ads, u(t_ads, z)) * yE_feed * Dt_ads * CPRT                                          ;
MBEq08(z)$z_end(z)                       ..      PE_ads_out      =e= sum(t_ads, P(t_ads, z) * u(t_ads, z) * yE(t_ads, z)) * Dt_ads * CPRT                                ;
MBEq09(z)$z_left(z)                      ..      PE_rin_feed     =e= P_high * sum(t_rin, u(t_rin, z)) * yE_rinse * Dt_rin * CPRT                                         ;
MBEq10(z)$z_end(z)                       ..      PE_rin_out      =e= sum(t_rin, P(t_rin, z) * u(t_rin, z) * yE(t_rin, z)) * Dt_rin * CPRT                                ;
MBEq11(z)$z_start(z)                     ..      PE_des_exhaust  =e= - sum(t_des, P(t_des, z) * u(t_des, z) * yE(t_des, z)) * Dt_des * CPRT                              ;
*****************************************
*       separation specifications       *
*****************************************
SSEq00                                   ..      PurityE         =e= PE_des_exhaust / (PA_des_exhaust + PE_des_exhaust + 1E-10)                                          ;
SSEq01                                   ..      PurityA         =e= (PA_ads_out + PA_rin_out) / (PA_ads_out + PA_rin_out + PE_ads_out + PE_rin_out + 1E-10)             ;
SSEq02                                   ..      RecoveryE       =e= (PE_des_exhaust - PE_rin_feed) / (PE_pre_feed + PE_ads_feed + 1E-10)                                ;
SSEq03                                   ..      RecoveryA       =e= (PA_ads_out + PA_rin_out) / (PA_pre_feed + PA_ads_feed + 1E-10)                                     ;
SSEq04                                   ..      PE_product      =e= PE_des_exhaust - PE_rin_feed                                                                        ;
SSEq05                                   ..      Mass_balanceA   =e= (PA_pre_feed + PA_ads_feed + PA_rin_feed) - (PA_ads_out + PA_rin_out + PA_des_exhaust)              ;
SSEq06                                   ..      Mass_balanceE   =e= (PE_pre_feed + PE_ads_feed + PE_rin_feed) - (PE_ads_out + PE_rin_out + PE_des_exhaust)              ;
SSEq07                                   ..      PurityE         =g= 0.99                                                                                                ;
SSEq08                                   ..      RecoveryE       =g= 0.3                                                                                                  ;
*****************************************
*          energy consumption           *
*****************************************
*ECEq00                                   ..      W_12_feed       =e= ERTI * (PA_pre_feed + PA_ads_feed + PE_pre_feed + PE_ads_feed) * ((P_high / P_input)**(1/IC) - 1)     ;
*ECEq00                                   ..      W_12_feed       =e= ERTI * (PA_pre_feed + PA_ads_feed + PE_pre_feed + PE_ads_feed) * ((P_atm / P_high)**(1/IC) - 1)     ;
ECEq00                                   ..      W_12_feed       =e= 0;

*ECEq01                                   ..      W_23_out        =e= ERTI * (PA_ads_out + PA_rin_out + PE_ads_out + PE_rin_out) * ((P_atm / P_end)**(1/IC) - 1)          ;
ECEq01                                   ..      W_23_out        =e= 0                                                                                                   ;

ECEq02                                   ..      W_3_feed        =e= ERTI * (PA_rin_feed + PE_rin_feed) * ((P_high / P_low)**(1/IC) - 1)                                 ;
ECEq03                                   ..      W_4_out         =e= ERTI * (PA_des_exhaust + PE_des_exhaust) * ((P_atm / P_low)**(1/IC) - 1)                            ;
ECEq04                                   ..      OBJ             =e= (W_12_feed + W_23_out + W_3_feed + W_4_out) / (1e-10 + PE_product)                                  ;

***********************************************************************************

Model PSA_BenchMark  /all/;

option NLP = CONOPT4;
option reslim = 36000;
option domlim = 500000;

*solve PSA_BenchMark using NLP maximizing PurityE;
*solve PSA_BenchMark using NLP maximizing RecoveryE;
solve PSA_BenchMark using NLP minimizing OBJ;

option OBJ:6;
display OBJ.l;






























