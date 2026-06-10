import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from pathlib import Path
import pandas as pd

def plot_bars(dict_data,
            kind = 'bar',
           title='',
           xlabel='',
           ylabel='',
           xmin=-5,
           xmax=25,
           bar_width=.7,
           color_dict=None,
           rot_xlabel=45,
           ylim=(-.15,1),
           leg=True,
           ncol=1, 
           bbox=(.25,.5,.5,.5),
           show_hline=True,
           y_hline=0,
           fontsize=18,
           figsize=(15, 8),
           dir_name=None,
           file_name=None,
           show=False,
           save=False):
    
        df = pd.DataFrame(dict_data)    
    
        ax = df.plot(kind=kind,
                     width=bar_width,
                     color=color_dict, 
                     rot=rot_xlabel,
                     ylim=ylim,
                     fontsize=fontsize,
                     figsize=figsize)
    
        ax.set_title(f'{title}',
                     pad=8,
                     fontdict={'fontsize':fontsize})
        
        ax.set_xlabel(xlabel,
                      fontdict={'fontsize':fontsize})
        
        ax.set_ylabel(ylabel,
                      fontdict={'fontsize':fontsize})

        if show_hline:
            ax.hlines(y_hline,
                      xmin,
                      xmax,
                      color='k',
                      linestyles='dashed')
            
        ax.get_legend().remove()
        if leg:
            ax.legend(loc='best',
                      bbox_to_anchor=bbox,
                      frameon=False,
                      ncol=ncol,
                      fontsize=fontsize-2)   
            
        if save:
            Path(dir_name).mkdir(parents=True,
                                exist_ok=True)
            plt.savefig(f'{dir_name}/{file_name}.png',
                        bbox_inches='tight',
                        dpi=300)          
        if show is False:
            plt.close()
