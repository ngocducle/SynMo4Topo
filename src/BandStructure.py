import matplotlib.pyplot as plt

### ============================================================================== ###
###                                                                                ###
### This module contains all the functions to plot the band structure for          ###
### different cases: in the whole Brillouin zone, in the vicinity of the M-point   ###
###                                                                                ###
### ============================================================================== ###

##### FUNCTION: plot the band structure in the whole Brillouin zone
def PlotBand_BrillouinZone(number,freqs,Nk,namesave,show_fig):
    fig, ax = plt.subplots()
    ax.plot(number, freqs)
    plt.vlines(Nk+1,0.0,1.0,linestyle='dashed',color='black')
    plt.vlines(2*(Nk+1),0.0,1.0,linestyle='dashed',color='black')
    plt.xlim(0,3*(Nk+1))
    plt.ylim(0,0.5)
    tick_locs = [i*(Nk+1) for i in range(4)]
    tick_labs = [r'$\Gamma$','X','M',r'$\Gamma$']
    ax.set_xticks(tick_locs)
    ax.set_xticklabels(tick_labs,size=16)
    ax.set_ylabel(r'$\omega a / (2 \pi c)$', fontsize = 14)
    plt.title(namesave,fontsize=14)
    plt.savefig(namesave+'.png')
    if show_fig == 'Yes':
        plt.show()

##### FUNCTION: plot the band structure in the vicinity of the M-point
def PlotBand_M(number,freqs,Nk,namesave,show_fig):
    fig, ax = plt.subplots()
    ax.plot(number, freqs)
    plt.vlines(Nk+1,0.0,1.0,linestyle='dashed',color='black')
    plt.vlines(2*(Nk+1),0.0,1.0,linestyle='dashed',color='black')
    plt.xlim(0,3*(Nk+1))
    plt.ylim(0,0.5)
    tick_locs = [i*(Nk+1) for i in range(4)]
    tick_labs = [r'$\Gamma \prime$','X $\prime$','M $\prime$',r'$\Gamma \prime$']
    ax.set_xticks(tick_locs)
    ax.set_xticklabels(tick_labs,size=16)
    ax.set_ylabel(r'$\omega a / (2 \pi c)$', fontsize = 14)
    plt.title(namesave,fontsize=14)
    plt.savefig(namesave+'.png')
    if show_fig == 'Yes':
        plt.show()


##### FUNCTION: plot the band structure along the Gamma-M point
def PlotBand_GammaM(number,freqs,Nk,namesave,show_fig):
    fig, ax = plt.subplots()
    ax.plot(number, freqs)
    tick_locs = [0,Nk+1]
    tick_labs = [r'$\omega a / (2 \pi c)$']
    ax.set_xticks(tick_locs)
    ax.set_xticklabels(tick_labs,size=16)
    ax.set_ylabel(r'$\omega a / (2 \pi c)$', fontsize = 14)
    plt.title(namesave,fontsize=14)
    plt.savefig(namesave+'.png')
    if show_fig == 'Yes':
      plt.show()

##### FUNCTION: plot the band structure in the whole Brillouin zone 
###   for square cell with rhombus hole, 
###   breaking the C4 symmetry of the square lattice
def PlotBand_BrillouinZone_Scell_Rhole(number,freqs,Nk,namesave,show_fig):
    fig, ax = plt.subplots()
    ax.plot(number, freqs)
    plt.vlines(Nk+1,0.0,1.0,linestyle='dashed',color='black')       # X 
    plt.vlines(2*(Nk+1),0.0,1.0,linestyle='dashed',color='black')   # M 
    plt.vlines(3*(Nk+1),0.0,1.0,linestyle='dashed',color='black')   # Gamma 
    plt.vlines(4*(Nk+1),0.0,1.0,linestyle='dashed',color='black')   # Y  
    plt.xlim(0,5*(Nk+1))
    plt.ylim(0,0.5)
    tick_locs = [i*(Nk+1) for i in range(6)]
    tick_labs = [r'$\Gamma$','X','M',r'$\Gamma$','Y','M']
    ax.set_xticks(tick_locs)
    ax.set_xticklabels(tick_labs,size=16)
    ax.set_ylabel(r'$\omega a / (2 \pi c)$', fontsize = 14)
    plt.title(namesave,fontsize=14)
    plt.savefig(namesave+'.png')
    if show_fig == 'Yes':
        plt.show()