import numpy as np 
import matplotlib.pyplot as plt
from sympy.solvers import solve
from sympy import Symbol

U_e=[]; U_r=[]; U_f=[]; U_i=[]; sU_e=[]; sU_r=[]; sU_f=[]; sU_i=[]
Ti=[0];Pi=[0];hi=[0];nu=[0];si=[0];
######################Heat Transfer Coefficient##########################################################
u_e=0.46                     #External Wall [W/m2ꞏK]
u_r=0.437                    #Roof          [W/m2ꞏK]
u_f=0.437                    #Floor         [W/m2ꞏK]
u_i=0.75                     #Internal Wall [W/m2ꞏK]
x_e = Symbol('x_e')
x_r = Symbol('x_r')
x_f = Symbol('x_f')
x_i = Symbol('x_i')
xe=solve(1/(0.02/0.8   + 0.38/0.81     + 0.004/209.7     + x_e/0.025 + 1/8 + 1/23)-u_e,x_e)
xr=solve(1/(0.2/1.19   + 0.004/209.7   + x_r/0.025  +                  1/9 + 1/23)-u_r,x_r)
xf=solve(1/(0.04/0.32  + 0.004/209.7   + x_f/0.025  +                  1/9 + 1/23)-u_f,x_f)
xi=solve(1/(0.004/209.7+ x_i/0.025                  +                  1/9 + 1/23)-u_i,x_i)

ß_wet=1.075
ß_m=1.09
ß=0.003
λ_0=0.025
t_in1=12.5+273.15
λ_e1=ß_wet*ß_m*λ_0*(1+ß*t_in1)
λ_r1=ß_wet*ß_m*λ_0*(1+ß*t_in1)
λ_f1=ß_wet*ß_m*λ_0*(1+ß*t_in1)
λ_i=ß_wet*ß_m*λ_0*(1+ß*t_in1)
x_e1 = Symbol('x_e1')
x_r1 = Symbol('x_r1')
x_f1 = Symbol('x_f1')
x_i = Symbol('x_i')
xe1=solve(1/(0.02/0.8   + 0.38/0.81     + 0.004/209.7     + x_e1/λ_e1 + 1/8 + 1/23)-u_e,x_e1)
xr1=solve(1/(0.2/1.19   + 0.004/209.7   + x_r1/λ_r1  +                  1/9 + 1/23)-u_r,x_r1)
xf1=solve(1/(0.04/0.32  + 0.004/209.7   + x_f1/λ_f1  +                  1/9 + 1/23)-u_f,x_f1)
xi =solve(1/(0.004/209.7+ x_i/λ_i                    +                  1/9 + 1/23)-u_i,x_i)

print('\nIsulation is dry')
print('Thick of Ex Wall= ',round(xe[0]*1000,1),' [mm]')
print('Thick of Roof   = ',round(xr[0]*1000,1),' [mm]')
print('Thick of Floor  = ',round(xf[0]*1000,1),' [mm]')
print('Thick of In Wall= ',round(xi[0]*1000,1),' [mm]')

print('\nMoisture is within insulation')
print('Thick of Ex Wall 1      = ',round(xe1[0]*1000,2),'[mm]')
print('Thick of Roof cooling   = ',round(xr1[0]*1000,2),'[mm]')
print('Thick of Floor cooling  = ',round(xf1[0]*1000,2),'[mm]')
print('Thick of In Wall        = ',round(xi[0]*1000,2),'[mm]\n')
###############################Temparature_Conditions####################################################
T_a=37                       #Max ambient Temperature        [C]
T_c=12.5                     #Temperature of cooling room    [C]
T_A=25                       #Temperature of Ante room       [C]
T_g=20                       #Ground temperature assumed     [C]
T_r=T_a+5                    #Temperature of roof            [C]   Table 4.1
T_ew=T_a+3                   #Temperature of Walls East&West [C]   Table 4.1
T_s=T_a+2                    #Temperature of Walls South     [C]   Table 4.1
###############################Areas_of_Walls############################################################
A_N=43*3.49                  #Area of North Walls                           [Square meter]
A_EWr=3.49*2*12                 #Area of East and West walls of cooling room   [Square meter]
A_EWa=3.49*4.7*2                #Area of East and West walls of Ante room      [Square meter]
A_Sr=3.49*4.6*2                 #Area of South walls of Cooling room           [Square meter]
A_Sa=3.49*33.8                  #Area of South walls of Ante room              [Square meter]
A_Rr=43*12                   #Area of roof of Rooms                         [Square meter]
A_Ra=4.7*33.8                #Area of roof of Ante room                     [Square meter]
A_Fr=A_Rr                    #Area of roof of Rooms                         [Square meter]
A_Fa=A_Ra                    #Area of roof of Ante Room                     [Square meter]

################################Heat transfer through enclosure##########################################
Q_n =u_e* (T_a-T_c) *A_N                     #Heat lossing from north Walls
Q_s =u_e*((T_s-T_c) *A_Sr +(T_s-T_A) *A_Sa)  #Heat lossing from south Walls
Q_ew=u_e*((T_ew-T_c)*A_EWr+(T_ew-T_A)*A_EWa) #Heat lossing from E&W Walls
Q_r =u_r*((T_r-T_c) *A_Rr +(T_r-T_A) *A_Ra)  #Heat lossing from Roofs
Q_g =u_f* (T_g-T_c) *A_Fr                    #Heat lossing from ground
Q_c =u_e* (T_A-T_c) *A_Sa                    #Heat lossing between Ante & Rooms
Q_transfer=Q_n+Q_s+Q_ew+Q_r+Q_g+Q_c          #Total Heat lossing of Plant
print ('Q Transfer=',round(Q_transfer/1000,3), '[kW]')##Q_transfer=43*3*0.46*(T_a-T_c)+12*3*2*0.46*(T_ew-T_c)+4.7*3*2*0.46*(T_s-T_c)+43*12*0.46*0.95*(T_r-T_c)+43*12*0.46*0.95*(T_g-T_c)+34*3*0.46*(T_A-T_c)

################################Goods dηils############################################################
m_b=6*9*16*18   *3           #Mass of goods (Bananas) per room
h_resp=0.0597e-3             #heat of respiration [kW/kg]

################################Comperssor operation time per day########################################
h_day=20                     #Huors a day

################################Safe Factor##############################################################
S_F=0.3                      #[%]

################################Loads####################################################################
Q_respiration=m_b*h_resp*3600*24           #[kJ/Day]
Q_transfer=Q_transfer*3600*24/1000         #[kJ/Day]
Q_overall=Q_transfer+Q_respiration         #[kJ/Day]
Q_Refrigerator=Q_overall*(1+S_F)/(3600*h_day)/3
print('\nPlant dηils')
print('Load of Respiration  per Day              = '+str('%-8.3f'%float(Q_respiration))+'[kJ/Day]')
print('Load of Heat Transfer Lossing per Day     = '+str('%-8.f'%float(Q_transfer))+'[kJ/Day]')
print('Comperssor operation time per day:          '+str('%-8.f'%float(h_day))+'[h/day]')
print('Design factor                               '+str('%-8.f'%float(S_F*100))+'[%]')
print('Overall Cooling Capacity                  = '+str('%-8.f'%float(Q_overall))+'[kJ/Day]')
print('Final Overall Cooling Capacity            = '+str('%-8.3f'%float(Q_Refrigerator))+'[kW/unit]')

#########################################################################################################
#################Refrigerator Reheat Cycle Thermodenamic Properties (R134a)##############################
prprtz=[[6,3.622,253.909,0,0.9283],
[11,3.622,258.511,0.0579,0.9446],
[59.03,13.19,286.108,0,0.9446],
[50,13.19,275.272,0,0.9115],
[50,13.19,123.478,0,0.4418],
[45,13.19,115.751,0,0.4177],
[6,3.622,115.751,0,0.4334],
[6,3.622,253.909,0,0.9283]]
for i in prprtz:
    for j in i:
        if j==i[0]:
            Ti.append(j)
        if j==i[1]:
            Pi.append(j)
        if j==i[2]:
            hi.append(j)
        if j==i[3]:
            nu.append(j)
        if j==i[4]:
            si.append(j)

#########################################################################################################
#############################Thermodynamic Calculation###################################################
c=0.03
m=1
η_e_motor=0.98 
b=0.0025
p_i_fr=40 
T_ev=Ti[1]
T_cd=Ti[4]
T_1= Ti[2]+273.15
thη =T_1-(T_ev+273.15)
alpha=1.12
bη=0.5
q_dot_cold=round(Q_Refrigerator,5) ##[kW]
q_cold=hi[1]-hi[7]
q_v=q_cold/nu[2]
q_cd=hi[3]-hi[5]
w_cyc=hi[3]-hi[2]
m_dot_wf=q_dot_cold/q_cold
V_dot_real=m_dot_wf*nu[2]
lambda_c=1-c*((Pi[3]/Pi[1])**(1/m)-1)
lambda_w=((T_ev+273.15)+thη)/(alpha*(T_cd+273.15)+0.5*thη )
Lambda=lambda_c*lambda_w
V_dot_h=V_dot_real/Lambda
η_i=lambda_w+b*Ti[1]
P_dot_th=m_dot_wf*w_cyc
P_dot_i=P_dot_th/η_i
P_dot_fr=V_dot_h*p_i_fr
P_dot_e=P_dot_i+P_dot_fr
P_dot_el=P_dot_e/η_e_motor
COP_star_carnot_RM=(T_ev+273.15)/(T_cd-T_ev)
COP_th_RM=q_dot_cold/P_dot_th 
COP_real_RM=q_dot_cold/P_dot_e 
η_th_RM=COP_th_RM/COP_star_carnot_RM
η_real_RM=COP_real_RM/COP_star_carnot_RM

print('\nRefrigerant is R134a')
print('Refrigeration Capacity              [Q_cold] ='+str('%-10.3f' %float(q_dot_cold))+'[kW]')
print('Condenser Temperature                 [T_cd] ='+str('%-10.f' %float(T_cd))+'[C]')
print('Ambient average Temperature            [T_0] ='+str('%-10.f' %float(T_a))+'[C]')
print('Cooling rooms Temperature              [T_c] ='+str('%-10.1f' %float(T_c))+'[C]')
print('Evaporator Temperature                [T_ev] ='+str('%-10.f' %float(T_ev))+'[C]')
print('Specific refrigerating capacity     [q_cold] ='+str('%-8.3f' %float(q_cold))+'[kJ/kg]')
print('Specific refrigerating capacity        [q_v] ='+str('%-8.3f' %float(q_v))+'[kJ/m3]')
print('Specific work of the Compressur       [w_cm] ='+str('%-8.3f' %float(w_cyc))+'[kJ/kg]')
print('Specific Specific condenser heat      [q_cd] ='+str('%-8.3f' %float(q_cd))+'[kJ/kg]')
print('Mass flow rate of the working fluid      [m] ='+str('%-8.3f' %float(m_dot_wf))+'[kg/s]')
print('Real displacement of the compressor  [Vreal] ='+str('%-8.3f' %float(V_dot_real))+'[kg/s]')
print('Volumetric efficiency of the compr.      [λ] ='+str('%-8.3f' %float(Lambda))+'[kg/s]')
print('The volumetric efficiency              [λ_c] ='+str('%-8.3f' %float(lambda_c))+'[kg/s]')
print('The volumetric efficiency              [λ_w] ='+str('%-8.3f' %float(lambda_w))+'[kg/s]')
print('Theoretical displacement of the compr. [V_h] ='+str('%-8.3f' %float(V_dot_h))+'[kg/s]')
print('Theoretical power of the compressor   [P_th] ='+str('%-8.3f' %float(P_dot_th))+'[kW]')
print('Indicator power of the compressor      [P_i] ='+str('%-8.3f' %float(P_dot_i))+'[kW]')
print('Indicator efficiency                   [η_i] ='+str('%-8.3f' %float(η_i))+'[kW]')
print('Power of friction of the compressor   [P_fr] ='+str('%-8.3f' %float(P_dot_fr))+'[kW]')
print('Effective power of the compressor      [P_e] ='+str('%-8.3f' %float(P_dot_e))+'[kW]')
print('Required electrical power of Compr.   [P_el] ='+str('%-8.3f' %float(P_dot_el))+'[kW]')

print('COP of the Carnot Cycle  =',round(COP_star_carnot_RM,3))
print('Theoretical value of COP =',round(COP_th_RM,3))
print('Real value of COP        =',round(COP_real_RM,3))
print('Theoretical effectiveness of the refrigeration cycle  =',round(η_th_RM,3))
print('Real effectiveness of the refrigeration cycle         =',round(η_real_RM,3))

plt.subplot(2, 1, 1)
plt.plot(si[1:], np.array(Ti[1:]),'-r',si[1:], np.array(Ti[1:]),'.k',linewidth = 1)
plt.title('R134a')
plt.xlabel('Entropt [kJ/kg/K]').set_color('red')
plt.ylabel('Temperature [K]')
plt.grid('on')

plt.subplot(2, 1, 2)
plt.plot(hi[1:], np.array(Pi[1:]),'-b',hi[1:], np.array(Pi[1:]),'.k')
plt.xlabel('Enthalpt [kJ/kg]').set_color('blue')
plt.ylabel('Pressure [Bar]')
plt.grid()
plt.show()
##input('Press Enter')
