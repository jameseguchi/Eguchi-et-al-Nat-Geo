# Python Script to run C-O cycle model from Eguchi et al. 2019 Nature Geoscience. Model is designed to track 
#C fluxes, reservoirs, and isotopes, as well as atmospheric O levels. Model was
#designed to investigate the relationship between large atmospheric
#oxygenation events and large C isotope excursions. In the model, these events are
#driven by changing the C emissions at Mid-ocean ridges, which in turn will
#change the flux of C leaving the atmosphere as carbonates and organic C.
#Changes in organic C reservoirs will drive changes in atmospheric O levels.
#C isotope excursion is driven by the relatively quick release of carbonates at
#arcs and the delayed release of organic C at ocean island volcanoes. See
#Eguchi et al. 2019 Nature Geoscience for more model details.
#

import numpy as np
import matplotlib.pyplot as plt

# time domain
# model time domain
t0 = 0       # Model start time [Myr]   
tf = 5000    # Model end time [Myr]
t  = np.linspace(t0,tf,(tf-t0+1))
tnum=len(t)
onset=1000 #time of first tectonic transtion [Myr]
tchange=2700 #time of second tectonic change [Myr]
crb_tau=30 #delay time for release of carbonates at arcs
org_tau=350#delay time for release of organic C at OIBs

# Constants
k=0.1 # weathering constant
forg=0.2 #fraction of C from atmosphere as organics
xcrb=1-forg # fraction of carbonates
chi_org=0.6 # fraction of organics that are subducted from the surface
alpha_org=0.0 # fraction of subducted organics that come out at arcs
orgc4=0.0 # fraction that remains in the mantle
chi_carb=chi_org # fraction of carbonates that are subducted from teh surface
alpha_crb=1.0 # fraction of subducted carbonates that come out at arcs
crbc4=0.0 # fraction that remains in the mantle
arcc1=0.0 # scalar of how much primitive c outgasses at arc

# initial conditions
c_atmi=0 # mass of atm-ocean C reservoir [g]
c_crbi=0 # mass of crustal carbonate C reservoir [g]
c_orgi=0 # mass of crustal organic C reservoir [g]
c_mcrbi=0 # mass of mantle carbonate C reservoir [g]
c_morgi=0 # mass of mantle organic C reservoir [g]
c_mntli=1e23 # mass of mantle primordial C reservoir [g]
d13C_atmi=-5 # d13c of C in atm-ocean [permill]
d13C_crbi=0  # d13c of C in carbonate [permill]
d13C_orgi=-25 # d13c of C in organic C [permill]
d13C_prim=-5 # d13c of C in primitive mantle C [permill]


F_mori=1e13#intial MORB C flux [g/Myr]
F_oibi=1e13#initial OIB flux [g/Myr]
F_arci=1e13#initial arc flux [g/Myr]

#Changes to MORB C fluxes to drive transitions

morb_change1=1e16 # MORB flux after Initial transition [g/Myr]
morb_change2=1e19 #MORB flux after 2nd tectonic transition [g/Myr]



#-----------No Need to change anything after this line-------------------------

F_orgi=forg*k*c_atmi
F_crbi=xcrb*k*c_atmi
F_sorgi=0
F_scrbi=0




# reservoirs
c_atm=np.zeros(tnum)
c_atm[0]=c_atmi
o_atm=np.zeros(tnum)
c_crb=np.zeros(tnum)
c_crb[0]=c_crbi
c_org=np.zeros(tnum)
c_org[0]=c_orgi
c_mcrb=np.zeros(tnum)
for index, item in enumerate(c_mcrb):
        c_mcrb[index] = c_mcrbi
c_morg=np.zeros(tnum)
for index, item in enumerate(c_morg):
        c_morg[index] = c_morgi
c_mntl=np.zeros(tnum)
c_mntl[0]=c_mntli

# isotopes
d13C_atm=np.zeros(tnum)
d13C_atm[0]=d13C_atmi
d13C_crb=np.zeros(tnum)
d13C_crb[0]=d13C_crbi
d13C_org=np.zeros(tnum)
d13C_org[0]=d13C_orgi
d13C_prm=np.ones(tnum)
d13C_prm=d13C_prm*d13C_prim
d13C_arc=np.zeros(tnum)
d13C_oib=np.zeros(tnum)
d13C_oib[0]=d13C_prim
d13C_mor=np.zeros(tnum)
for index, item in enumerate(d13C_mor):
        d13C_mor[index] = -5

# fluxes
F_org=np.zeros(tnum)
F_org[0]=F_orgi
F_crb=np.zeros(tnum)
F_crb[0]=F_crbi
F_sorg=np.zeros(tnum)
F_sorg[0]=F_sorgi
F_scrb=np.zeros(tnum)
F_scrb[0]=F_scrbi
F_mor=np.zeros(tnum)
for index, item in enumerate(F_mor):
        F_mor[index] =F_mori
F_mor[0]=F_mori
F_oib=np.zeros(tnum)
F_oib[0]=F_oibi
F_arc=np.zeros(tnum)
F_arc[0]=F_arci
for index, item in enumerate(F_arc):
        F_arc[index]=F_arci
F_tot=np.zeros(tnum)
F_tot[0]=F_oib[0]+F_arc[0]+F_mor[0]


# iterate through model
# preconvection
for time in range(1,tnum):

    if (t0+time)<(onset):

        c_atm[time]=c_atm[time-1]+(F_tot[time-1]-(F_org[time-1]+F_crb[time-1]))
        c_mntl[time]=c_mntl[time-1]-F_oib[time-1]
        c_crb[time]=c_crb[time-1]+(F_crb[time-1]-F_scrb[time-1])
        c_org[time]=c_org[time-1]+(F_org[time-1]-F_sorg[time-1])

        F_oib[time]=F_oibi
        F_arc[time]=F_arci
        F_tot[time]=F_oib[time]
        F_org[time]=forg*k*c_atm[time]
        F_crb[time]=xcrb*k*c_atm[time]
        F_tot[time]=F_oib[time]+F_arc[time]+F_mor[time]

        d13C_oib[time]=d13C_prim
        d13C_atm[time]=d13C_oib[time]
        d13C_crb[time]=d13C_atm[time]+5
        d13C_org[time]=d13C_atm[time]-20

    elif (t0+time)<(onset+crb_tau):
        atm_Fout=F_org[time-1]+F_crb[time-1]

        c_atm[time]=c_atm[time-1]+(F_tot[time-1]-(atm_Fout))
        c_crb[time]=c_crb[time-1]+(F_crb[time-1]-chi_carb*F_scrb[time-1])
        c_org[time]=c_org[time-1]+(F_org[time-1]-chi_org*F_sorg[time-1])
        c_mcrb[time]=c_mcrb[time-1]+F_scrb[time-1]
        c_morg[time]=c_morg[time-1]+F_sorg[time-1]
        c_mntl[time]=c_mntl[time-1]-F_oib[time-1]-F_mor[time-1]

        F_oib[time]=F_oibi
        F_mor[time]=morb_change1
        F_org[time]=forg*k*c_atm[time]
        F_crb[time]=xcrb*k*c_atm[time]
        F_scrb[time]=chi_carb*F_crb[time]
        F_sorg[time]=chi_org*F_org[time]
        F_tot[time]=F_oib[time]+F_arc[time]+F_mor[time]

        d13C_oib[time]=d13C_prim
        d13C_atm[time]=d13C_oib[time]
        d13C_crb[time]=d13C_atm[time]+5
        d13C_org[time]=d13C_atm[time]-20

    # outgassing only at arcs
    elif (t0+time)<(onset+org_tau):

        atm_Fout=F_org[time-1]+F_crb[time-1]

        c_atm[time]=c_atm[time-1]+(F_tot[time-1]-(atm_Fout))
        c_crb[time]=c_crb[time-1]+(F_crb[time-1]-F_scrb[time-1])
        c_org[time]=c_org[time-1]+(F_org[time-1]-F_sorg[time-1])
        c_mcrb[time]=c_mcrb[time-1]+F_scrb[time-1]-(1-crbc4)*alpha_crb*F_scrb[time-crb_tau]
        c_morg[time]=c_morg[time-1]+F_sorg[time-1]-(1-orgc4)*alpha_org*F_sorg[time-crb_tau]
        c_mntl[time]=c_mntl[time-1]-F_oib[time-1]-F_mor[time-1]-arcc1*(alpha_crb*F_scrb[time-crb_tau]+alpha_org*F_sorg[time-crb_tau])

        F_oib[time]=F_oibi
        F_mor[time]=morb_change1
        F_org[time]=forg*k*c_atm[time]
        F_crb[time]=xcrb*k*c_atm[time]
        F_scrb[time]=chi_carb*F_crb[time]
        F_sorg[time]=chi_org*F_org[time]
        Farccrb=alpha_crb*F_scrb[time-crb_tau]
        Farcorg=alpha_org*F_sorg[time-crb_tau]
        Farcmntl=(arcc1)*(Farccrb+Farcorg)
        F_arc[time]=Farcmntl+Farccrb+Farcorg
   
        F_tot[time]=F_oib[time]+F_arc[time]+F_mor[time]

        d13C_oib[time]=d13C_prim
        d13C_arcorg=Farcorg/F_arc[time]*d13C_org[time-crb_tau]
        d13C_arccrb=Farccrb/F_arc[time]*d13C_crb[time-crb_tau]
        d13C_arcmntl=Farcmntl/F_arc[time]*d13C_prm[time]
        d13C_arc[time]=d13C_arcorg+d13C_arccrb+d13C_arcmntl
        d13C_atm[time]=F_oib[time]/F_tot[time]*d13C_oib[time]+F_arc[time]/F_tot[time]*d13C_arc[time]+F_mor[time]/F_tot[time]*d13C_mor[time]
        d13C_crb[time]=d13C_atm[time]+5
        d13C_org[time]=d13C_atm[time]-20

    # have all systems going
    elif ((t0+time)>=(onset+org_tau)) and ((t0+time)<tchange):
    #else:

        atm_Fout=F_org[time-1]+F_crb[time-1]

        c_atm[time]=c_atm[time-1]+(F_tot[time-1]-(atm_Fout))
        c_crb[time]=c_crb[time-1]+(F_crb[time-1]-F_scrb[time-1])
        c_org[time]=c_org[time-1]+(F_org[time-1]-F_sorg[time-1])
        c_mcrb[time]=c_mcrb[time-1]+F_scrb[time-1]-(1-crbc4)*alpha_crb*F_scrb[time-crb_tau]-(1-crbc4)*(1-alpha_crb)*F_scrb[time-org_tau]
        c_morg[time]=c_morg[time-1]+F_sorg[time-1]-(1-orgc4)*alpha_org*F_sorg[time-crb_tau]-(1-orgc4)*(1-alpha_org)*F_sorg[time-org_tau]
        c_mntl[time]=c_mntl[time-1]-F_oib[time-1]*0-F_mor[time-1]-arcc1*(alpha_crb*F_scrb[time-crb_tau]+alpha_org*F_sorg[time-crb_tau])-arcc1*(crbc4*F_scrb[time-org_tau]+orgc4*F_sorg[time-org_tau])

        F_mor[time]=morb_change1
        F_org[time]=forg*k*c_atm[time]
        F_crb[time]=xcrb*k*c_atm[time]
        F_scrb[time]=chi_carb*F_crb[time]
        F_sorg[time]=chi_org*F_org[time]
        Farccrb=alpha_crb*F_scrb[time-crb_tau]
        Farcorg=alpha_org*F_sorg[time-crb_tau]
        Farcmntl=(arcc1)*(Farccrb+Farcorg)
        F_arc[time]=Farcmntl+Farccrb+Farcorg
        Foibcrb=(1-alpha_crb)*F_scrb[time-org_tau]
        Foiborg=(1-alpha_org)*F_sorg[time-org_tau]
        Foibmntl=F_oibi
        F_oib[time]=Foibcrb+Foiborg+Foibmntl
        
        F_tot[time]=F_oib[time]+F_arc[time]+F_mor[time]

        d13C_oiborg=Foiborg/F_oib[time]*d13C_org[time-org_tau]
        d13C_oibcrb=Foibcrb/F_oib[time]*d13C_crb[time-org_tau]
        d13C_oibmntl=Foibmntl/F_oib[time]*d13C_prm[time]
        d13C_oib[time]=d13C_oiborg+d13C_oibcrb+d13C_oibmntl
        d13C_arcorg=Farcorg/F_arc[time]*d13C_org[time-crb_tau]
        d13C_arccrb=Farccrb/F_arc[time]*d13C_crb[time-crb_tau]
        d13C_arcmntl=Farcmntl/F_arc[time]*d13C_prm[time]
        d13C_arc[time]=d13C_arcorg+d13C_arccrb+d13C_arcmntl
        d13C_atm[time]=F_oib[time]/F_tot[time]*d13C_oib[time]+F_arc[time]/F_tot[time]*d13C_arc[time]+F_mor[time]/F_tot[time]*d13C_mor[time]
        d13C_crb[time]=d13C_atm[time]+5
        d13C_org[time]=d13C_atm[time]-20
    
    elif (t0+time)>=tchange:
    #else:

        atm_Fout=F_org[time-1]+F_crb[time-1]
                
        F_mor[time]=morb_change2
        c_atm[time]=c_atm[time-1]+(F_tot[time-1]-(atm_Fout))
            
        c_crb[time]=c_crb[time-1]+(F_crb[time-1]-F_scrb[time-1])
        c_org[time]=c_org[time-1]+(F_org[time-1]-F_sorg[time-1])
        c_mcrb[time]=c_mcrb[time-1]+F_scrb[time-1]-(1-crbc4)*alpha_crb*F_scrb[time-crb_tau]-(1-crbc4)*(1-alpha_crb)*F_scrb[time-org_tau]
        c_morg[time]=c_morg[time-1]+F_sorg[time-1]-(1-orgc4)*alpha_org*F_sorg[time-crb_tau]-(1-orgc4)*(1-alpha_org)*F_sorg[time-org_tau]
        c_mntl[time]=c_mntl[time-1]-F_oib[time-1]*0-F_mor[time-1]-arcc1*(alpha_crb*F_scrb[time-crb_tau]+alpha_org*F_sorg[time-crb_tau])-arcc1*(crbc4*F_scrb[time-org_tau]+orgc4*F_sorg[time-org_tau])
         
        F_org[time]=forg*k*c_atm[time]
        F_crb[time]=xcrb*k*c_atm[time]
        F_scrb[time]=chi_carb*F_crb[time]
        F_sorg[time]=chi_org*F_org[time]
        Farccrb=alpha_crb*F_scrb[time-crb_tau]
        Farcorg=alpha_org*F_sorg[time-crb_tau]
        Farcmntl=(arcc1)*(Farccrb+Farcorg)
        F_arc[time]=Farcmntl+Farccrb+Farcorg
        Foibcrb=(1-alpha_crb)*F_scrb[time-org_tau]
        Foiborg=(1-alpha_org)*F_sorg[time-org_tau]
        Foibmntl=F_oibi
        F_oib[time]=(Foibcrb+Foiborg+Foibmntl)
        
        F_tot[time]=F_oib[time]+F_arc[time]+F_mor[time]

        d13C_oiborg=Foiborg/F_oib[time]*d13C_org[time-org_tau]
        d13C_oibcrb=Foibcrb/F_oib[time]*d13C_crb[time-org_tau]
        d13C_oibmntl=Foibmntl/F_oib[time]*d13C_prm[time]
        d13C_oib[time]=d13C_oiborg+d13C_oibcrb+d13C_oibmntl
        d13C_arcorg=Farcorg/F_arc[time]*d13C_org[time-crb_tau]
        d13C_arccrb=Farccrb/F_arc[time]*d13C_crb[time-crb_tau]
        d13C_arcmntl=Farcmntl/F_arc[time]*d13C_prm[time]
        d13C_arc[time]=d13C_arcorg+d13C_arccrb+d13C_arcmntl
        d13C_atm[time]=F_oib[time]/F_tot[time]*d13C_oib[time]+F_arc[time]/F_tot[time]*d13C_arc[time]+F_mor[time]/F_tot[time]*d13C_mor[time]
        d13C_crb[time]=d13C_atm[time]+5
        d13C_org[time]=d13C_atm[time]-20

# visualize data
fig=plt.figure(1,[4,12])
plt.subplot(4,1,1)
plt.plot(t,d13C_crb)
plt.xlim([t0,tf])
plt.ylim([-10,15])
plt.title('Isotopes')
plt.legend(frameon=False)


plt.subplot(4,1,2)
plt.semilogy(t,(c_morg+c_org)/1e21)
plt.xlim([t0,tf])
plt.title(r'Atmospheric $O_2$')
plt.legend(frameon=False)


plt.subplot(4,1,3)
plt.semilogy(t,F_oib,'k',label='oib')
plt.semilogy(t,F_arc,'b',label='arc')
plt.semilogy(t,F_mor,'y',label='mor')
plt.semilogy(t,F_crb+F_org,'r',label='weathering')
plt.semilogy(t,F_tot,'g',label='tot volc')
plt.xlim([t0,tf])
plt.title('Fluxes')
plt.legend(frameon=False)


plt.subplot(4,1,4)
plt.semilogy(t,c_atm,'b',label='c atm ocean')
plt.semilogy(t,c_org,'r',label='c crustal org')
plt.semilogy(t,c_morg,'g',label='c mantle org')
plt.semilogy(t,c_mntl,'y',label='Prim mantle C')
plt.semilogy(t,c_crb,'m',label='crustal carb')
plt.semilogy(t,c_mcrb,'c',label='mantle carb')
plt.legend(frameon=False)
plt.xlim([t0,tf])
plt.show()
