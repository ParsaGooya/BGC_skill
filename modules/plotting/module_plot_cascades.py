import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator,
                               FormatStrFormatter,
                               AutoMinorLocator)
from matplotlib.collections import LineCollection
from pathlib import Path


class set_plot_arguments():
    def __init__(self, stype, xx_max, xx_i, xx_step, units=''):
        if stype == 'cascade':
            self.xlabel = ''
            self.ylabel = f'{units}' #r'Pg C yr$^{-1}$'
            self.label_ts1 = 'Observed'
            self.label_ts2 = 'Uninitialized'
            self.label_ts3 = 'Initialized'
            self.color_ts1 = 'k'
            self.color_ts2 = 'gray'
            self.color_ts3 = 'r'
            self.linestyle_ts1 = 'solid' 
            self.linestyle_ts2 = 'solid'
            self.linestyle_ts3 = 'None'
            self.linewidth_ts1 = 5
            self.linewidth_ts2 = 5
            self.linewidth_ts3 = 0
            self.marker_ts1 = 'o'                 
            self.marker_ts2 = 'o'                 
            self.marker_ts3 = 'o'                 
            self.markersize_ts1 = 0
            self.markersize_ts2 = 0
            self.markersize_ts3 = 10
            self.show_all_leads = True                 
            self.xx_ticks_loctn = np.arange(0,xx_max,step=xx_step)
            self.xx_ticks_label = np.arange(0,xx_max,step=xx_step) + xx_i
            self.nx_ticks_minor = 0
            self.show_grid = False
            self.ymin = 9. 
            self.ymax = 13. 
            self.figsizex = 24
            self.figsizey = 8
        else:
            print(stype+' --type not defined')
            quit()



def colored_line(x,y,line_width,
                  cmap='jet_r'):
    from matplotlib.collections import LineCollection
    N = 501
    xx = np.linspace(x.min(),
                     x.max(), N)
    yy = np.interp(xx, x, y)
    points = np.array([xx, yy]).T.reshape(-1,1,2)
    segments = np.concatenate([points[:-1],
                               points[1:]],
                              axis=1)
    lc = LineCollection(segments,
                        cmap=cmap,
                        linewidth=line_width)
    lc.set_array(xx)
    return lc



def get_time_series_to_plot_cascade(obs,
                                    ds_sim=None,
                                    ds_hnd=None):
    '''
     Data must be 2D of the form: initial years x lead time in month
    '''
    
    nmth = 12
    n1, n2 = obs.shape
    nld = max([1,int(n2/nmth)])
    
    ts1 = np.zeros(n1+nld+1) - np.nan
    ts1[1:n1+1]  = obs.reshape(n1,
                               nld,
                               nmth).mean(axis=2)[:,0]
    
    ts2 = np.zeros(n1+nld+1) - np.nan
    ts3 = np.zeros((n1+nld,nld)) - np.nan
    
    if ds_sim is not None:
        n1_sim, n2_sim = ds_sim.shape
        ts2[1:n1+1]  = ds_sim.reshape(n1,
                                      nld,
                                      nmth).mean(axis=2)[:,0]
    if ds_hnd is not None:
        n1_hnd, n2_hnd = ds_hnd.shape
        ts3[1:n1+1,:] = ds_hnd.reshape(n1,
                                       nld,
                                       nmth).mean(axis=2)     
    
    return [ts1,
            ts2,
            ts3]


def plot_cascade(ts1=None,
                 ts2=None,
                 ts3=None,
                 title='',
                 uncertainty=False,
                 units='',
                 ymin=None,
                 ymax=None,
                 xx_max=40,
                 xx_i=1982,
                 xx_step=5,
                 text1='', 
                 xtext1=None,                 
                 xtext2=None,                 
                 ytext1=None,                 
                 ytext2=None,                 
                 text2='',
                 shadesize_left=None,
                 shadesize_right=None,
                 show_samples=False,
                 fig_dir='./Figures',
                 fig_name=None,
                 save=False):
    
    fnt_size = 24
    
    set_arg = set_plot_arguments('cascade',
                                 xx_max,
                                 xx_i,
                                 xx_step,
                                 units)        

    fig, ax = plt.subplots(figsize=(set_arg.figsizex,
                                    set_arg.figsizey))
    
    if ymin is None and ymax is None:
        ymin = set_arg.ymin
        ymax = set_arg.ymax    
    ax.set_ylim([ymin,ymax])
    ax.tick_params(which='major',
                   width=3,
                   length=10)
    if set_arg.nx_ticks_minor > 0:
        ax.xaxis.set_minor_locator(AutoMinorLocator(set_arg.nx_ticks_minor))
        ax.tick_params(which='minor',
                       width=3,
                       length=5)
    plt.setp(ax.spines.values(),
             linewidth=4)
    if title is not None:
        plt.title(title,
                  size=fnt_size+15)
    plt.xlabel(set_arg.xlabel,
               size=fnt_size+8)
    plt.ylabel(set_arg.ylabel,
               size=fnt_size+8)
    plt.xticks(size=fnt_size+4)
    plt.yticks(size=fnt_size+4)
    if ~np.isnan(np.max(set_arg.xx_ticks_label)) == True:
        plt.xticks(set_arg.xx_ticks_loctn,
                   set_arg.xx_ticks_label,
                   size=fnt_size+4)
        
    if show_samples:
#        ax.axvspan(1,shadesize_left,
#                   alpha=0.1,color='blue')        
        ax.axvspan(shadesize_left,shadesize_left+shadesize_right,
                   alpha=0.1,color='red')        
        if xtext1 is None:
            xtext1 = 2
        if xtext2 is None:
            xtext2 = shadesize_left + 1
        if ytext1 is None:
            ytext1 = ymin - .5 
        if ytext2 is None:
            ytext2 = ymax - .5
        # text1 = '' #'training set'
        # text2 = 'test years'
        plt.text(xtext1,
                 ytext1,
                 text1,
                 weight='bold',
                 fontsize=fnt_size+3)
        plt.text(xtext2,
                 ytext2,
                 text2,
                 fontsize=fnt_size+3)

### his
    if ts2 is not None: 
        plt.plot(ts2,
                 color=set_arg.color_ts2,
                 linestyle=set_arg.linestyle_ts2,
                 linewidth=set_arg.linewidth_ts2,
                 marker=set_arg.marker_ts2,                 
                 markersize=set_arg.markersize_ts2,
                 label=set_arg.label_ts2)
#### hnd
    if ts3 is not None: 
        nld=len(ts3[1])
        marker_ts = ["" for idot in range(nld)]
        marker_ts[0] = "o"
        plt.plot(ts3[:,0],
                 color=set_arg.color_ts3,
                 linestyle=set_arg.linestyle_ts3,
                 linewidth=set_arg.linewidth_ts3,
                 marker=set_arg.marker_ts3,                 
                 markersize=set_arg.markersize_ts3,
                 label=set_arg.label_ts3)
        if set_arg.show_all_leads:
            line_width = 3
            for iinic in range(len(ts3[:,0])):
                if np.isnan(ts3[iinic,0]) == False:
                    if iinic+1 < len(ts3[:,0]) and np.isnan(ts3[iinic+1,0]) == True:
                        line_width = 3 #8 #6
                    xx_inic   = np.array(np.arange(iinic,nld+iinic))
                    lc = colored_line(xx_inic,
                                      ts3[iinic,:],
                                      line_width)
                    plt.gca().add_collection(lc)
#### obs
    if ts1 is not None: #len(ts1) > 1:
        plt.plot(ts1,
                 color=set_arg.color_ts1,
                 linestyle=set_arg.linestyle_ts1,
                 linewidth=set_arg.linewidth_ts1,
                 marker=set_arg.marker_ts1,                 
                 markersize=set_arg.markersize_ts1,
                 label=set_arg.label_ts1)

    if ts1 is None:
        print('ts1 not provided --will not appear in plot')
    if ts2 is None:
        print('ts2 not provided --will not appear in plot')
    if ts3 is None:
        print('ts3 not provided --will not appear in plot')
        
    handles, labels = plt.gca().get_legend_handles_labels()
    
    order = [0]
    if ts1 is not None and ts3 is not None:
        order = [1,0]
    if ts1 is not None and ts2 is not None and ts3 is not None:
        order = [2,1,0]
        
    if ts1 is not None or ts2 is not None or ts3 is not None:
        plt.legend([handles[idx] for idx in order],
                   [labels[idx]  for idx in order],
                   fontsize=fnt_size-1,#-2,
                   frameon=False)
    
    if uncertainty:

        plt.plot(lr,
                 color=set_arg.color_ts1,
                 linestyle='dashed',
                    alpha = 0.5)
        
        plt.plot(ur,
                 color=set_arg.color_ts1,
                 linestyle='dashed',
                    alpha = 0.5)

    if save:
        Path(fig_dir).mkdir(parents=True, exist_ok=True)
        plt.savefig(f'{fig_dir}/{fig_name}',
                    bbox_inches='tight',
                    dpi=300)    
    plt.show()
    plt.close()                   


