"""General plotting settings for matplotlib
Example in:
https://matplotlib.org/3.3.2/tutorials/introductory/customizing.html
"""

import matplotlib
import matplotlib.pyplot as plt

RCPARAM = {
    # Default text/line sizes
    'lines.linewidth'   :   2,
    'lines.markersize'  :   10,
    'font.size'         :   12,

    'xtick.labelsize'   :   12,
    'xtick.direction'   :   'in',
    'ytick.labelsize'   :   12,
    'ytick.direction'   :   'in',
    'axes.labelsize'    :   12,
    'axes.titlesize'    :   12,

    # Latex
    'text.usetex'       : True,

    # Default figure size (important for 'plt.savefig'
    'figure.figsize': (4.7, 4.7),

    # Save figure options
    'savefig.dpi': 250,
    }
plt.style.use(RCPARAM)

