import matplotlib.pyplot as plt


##### FUNCTION: plot the band structure
def PlotBand_BrillouinZone(number,freqs,Nk,namesave):
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
    plt.show()