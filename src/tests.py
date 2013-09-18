import ACHP
import matplotlib.pyplot as plt
import numpy as np
from CoolProp.Plots import Ph

iT = ACHP.iT
iP = ACHP.iP
iQ = ACHP.iQ

B = ACHP.BrazedPlateHeatExchanger()

fig = plt.figure(figsize=(6,3))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
ax1.set_yscale('log')

for mdot_ratio_h_c in np.logspace(-2,2):
    # Cold stream inlet
    B.State_c_inlet = ACHP.CoolPropStateClassSI("Propane")
    B.State_c_inlet.update(iT,290,iQ,1)
    B.State_c_inlet.update(iT,270,iP,B.State_c_inlet.p())
    B.mdot_c = 0.01
    
    # Hot stream inlet
    B.State_h_inlet = ACHP.CoolPropStateClassSI("Propane")
    B.State_h_inlet.update(iT,320,iQ,1)
    B.State_h_inlet.update(iT,360,iP,B.State_h_inlet.p())
    B.mdot_h = B.mdot_c*mdot_ratio_h_c
    
    B.SaturationStates()
    B.DetermineQmax()
    
    ph = [B.State_h_inlet.p()/1000 for i in range(len(B.EnthalpyList_h))]
    pc = [B.State_c_inlet.p()/1000 for i in range(len(B.EnthalpyList_h))]
    hh = [h/1000 for h in B.EnthalpyList_h]
    hc = [h/1000 for h in B.EnthalpyList_c]
    hhnorm = (np.array(hh)-hh[0])/(hh[-1]-hh[0])
    hcnorm = (np.array(hc)-hc[0])/(hc[-1]-hc[0])
    Th = [T for T in B.TemperatureList_h]
    Tc = [T for T in B.TemperatureList_c]
    
    Ph('Propane',axis = ax1)
    ax1.plot(hh, ph, 'ro-')
    ax1.plot(hc, pc, 'bo-')
    ax1.set_title('')
    ax1.set_ylim(ymin = 100)
    
    ax2.plot(hcnorm,Tc,'bo-')
    ax2.plot(hcnorm,Th,'ro-')
    
    plt.tight_layout()    
    plt.savefig('{s:08d}'.format(s=int(mdot_ratio_h_c*1000))+'.png')    
    ax1.cla()
    ax2.cla()