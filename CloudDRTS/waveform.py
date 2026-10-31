from common import np, plt, figure, OUT

def main():
    fig,ax=plt.subplots(figsize=(7.7,2.5))
    for name,label in [('interface_reference','Unsplit'),('split_k20','Full split, +100 µs'),('thevenin_k0','Static equivalent, +0 µs')]:
        a=np.genfromtxt(OUT/(name+'.csv'),delimiter=',',names=True)
        ax.plot(a['t_s']*1000,a['vr_V'],label=label)
    ax.axvspan(5,9,alpha=.09,color='#37877a')
    ax.set(xlim=(3,11),xlabel='Time (ms)',ylabel='Receiving voltage (V)')
    ax.legend(ncol=3,fontsize=8); figure(fig,'waveform')
if __name__=='__main__': main()
