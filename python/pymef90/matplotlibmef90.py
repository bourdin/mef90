def matplotlibdefaults(palette='medium',useTex=False):
    from matplotlib import rcParams
    lightgrey = '#CBCBCB'
    grey = '#8C8C8C'
    darkgrey  = '#4D4D4D'
    if palette == 'pastel':
        cs = ['#fbb4ae', '#b3cde3', '#ccebc5', '#decbe4', '#fed9a6', '#ffffcc', '#e5d8bd', '#fddaec']
    elif palette == 'light':
        cs =  ['#8abde6', 'fbb258', '90cd97', 'f6aac9' , 'bfa554' , 'bc99c7' , 'eddd46' , 'f07e6e']
    elif palette == 'medium':
        cs = ['#5da5da', 'faa43a' , '60bd68' , 'f17cb0' , 'b2912f' , 'b276b2' , 'decf3f' , 'f15854']
    elif palette == 'dark':
        cs =   ['#265dab', '#df5dab', '#059748' , '#e5120b' , '#9d722a' , '#7b3a96' , '#c7b42e' , '#cb201e']
    else:
        print('Unknown palette: Using medium')
        cs = ['#5da5da', 'faa43a' , '60bd68' , 'f17cb0' , 'b2912f' , 'b276b2' , 'decf3f' , 'f15854']

    rcParams['axes.labelsize']   = 12
    rcParams['axes.facecolor']   = 'none'   # axes background color
    rcParams['axes.edgecolor']   = grey  # axes edge color
    rcParams['axes.labelcolor']  = darkgrey
    rcParams['axes.color_cycle'] = cs

    rcParams['xtick.labelsize']   = 12
    rcParams['ytick.labelsize']   = 12
    rcParams['xtick.color']       = grey
    rcParams['ytick.color']       = grey
    rcParams['xtick.direction']   = 'out'
    rcParams['ytick.direction']   = 'out'
    rcParams['xtick.major.width'] = 2
    rcParams['ytick.major.width'] = 2
    rcParams['xtick.major.size']  = 8
    rcParams['ytick.major.size']  = 8

    rcParams['legend.fontsize'] = 12
    rcParams['font.family'] = 'serif'
    rcParams['text.usetex'] = useTex
    if useTex:
        rcParams['font.serif'] = 'Computer Modern Roman'
    else:
        rcParams['font.serif'] = 'Times'

    rcParams['lines.linewidth'] = 3.0
    rcParams['lines.markersize'] = 8     
    rcParams['lines.markeredgewidth'] = 0
    rcParams['lines.solid_joinstyle'] = 'round'

    rcParams['figure.facecolor'] = '#FFFFFF'    # figure facecolor; 0.75 is scalar gray

        

    rcParams['axes.linewidth'] = 2.0
    rcParams['axes.titlesize'] = 12
    rcParams['text.color'] = darkgrey

    rcParams['grid.color'] = lightgrey
    rcParams['grid.linestyle'] = '-'
    rcParams['grid.linewidth'] = 0.25     # in points
    rcParams['grid.alpha'] = .5     # transparency, between 0.0 and 1.0

    rcParams['legend.frameon'] = False


def setspines():
    import matplotlib.pylab
    for i in matplotlib.pylab.get_fignums():
        for j in matplotlib.pylab.figure(i).get_axes():
            j.spines['top'].set_color('none')
            j.spines['right'].set_color('none')
            j.tick_params(axis='both',top='off',right='off',which='both',colors='#8C8C8C')
            #j.spines['left'].set_position(('outward',10))
            #j.spines['bottom'].set_position(('outward',10))
            #j.spines['left'].set_position(('axes', -0.05))
            #j.spines['bottom'].set_position(('axes', -0.05))
    return 0
