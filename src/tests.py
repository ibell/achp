import ACHP
import matplotlib.pyplot as plt
import numpy as np
from CoolProp.Plots import Ph

iT = ACHP.iT
iP = ACHP.iP
iQ = ACHP.iQ

B = ACHP.BrazedPlateHeatExchanger()

fig = plt.figure(figsize=(6,3))

for mdot_ratio_h_c in np.logspace(np.log10(0.1),np.log10(5),200):
    # Cold stream inlet
    B.State_c_inlet = ACHP.CoolPropStateClassSI("Propane")
    B.State_c_inlet.update(iT,290,iQ,0.3)
    B.State_c_inlet.update(iT,270,iP,B.State_c_inlet.p())
    B.mdot_c = 0.01
    
    # Hot stream inlet
    B.State_h_inlet = ACHP.CoolPropStateClassSI("Propane")
    B.State_h_inlet.update(iT,320,iQ,1)
    B.State_h_inlet.update(iT,360,iP,B.State_h_inlet.p())
    B.mdot_h = B.mdot_c*mdot_ratio_h_c
    
    B.SaturationStates()
    Qmax = B.DetermineQmax()
    
    ph = [B.State_h_inlet.p()/1000 for i in range(len(B.EnthalpyList_h))]
    pc = [B.State_c_inlet.p()/1000 for i in range(len(B.EnthalpyList_h))]
    hh = np.array([h/1000 for h in B.EnthalpyList_h])
    hc = np.array([h/1000 for h in B.EnthalpyList_c])
    hhnorm = B.mdot_h*(np.array(hh)-hh[0])/Qmax*1000
    hcnorm = B.mdot_c*(np.array(hc)-hc[0])/Qmax*1000
    Th = [T for T in B.TemperatureList_h]
    Tc = [T for T in B.TemperatureList_c]
    
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    ax1.set_yscale('log')

    Ph('Propane',axis = ax1)
    ax1.plot(hh, ph, 'ro-')
    ax1.plot(hc, pc, 'bo-')
    ax1.set_title('')
    ax1.set_ylim(ymin = 100)
    
    ax2.plot(hcnorm, Tc, 'bo-')
    ax2.plot(hhnorm, Th, 'ro-')
    ax2.set_xlim(0,1)
    
    plt.suptitle(str(mdot_ratio_h_c))
    
    plt.tight_layout()    
    plt.savefig('{s:08d}'.format(s=int(mdot_ratio_h_c*1000))+'.png')    
    fig.clf()