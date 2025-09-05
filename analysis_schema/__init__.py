from .schema import *
from .compute import SpikeDetectionCompute
import matplotlib as mpl

mpl.rcParams.update({
    # PDF backend settings
    "pdf.fonttype": 42,       # Embed fonts as Type 42 (TrueType)
    "ps.fonttype": 42,        # Same for PS
    "svg.fonttype": 'none',   # Keep text as text in SVG
    
    # Font settings
    #"font.family": "sans-serif",
    #"font.sans-serif": ["Arial"],  # Change to your preferred font
    "font.size": 12,          # Base font size
    
    # Axes and ticks
    "axes.labelsize": 12,
    "axes.titlesize": 16,
    "xtick.labelsize": 11,
    "ytick.labelsize": 12,
    
    # Lines and markers
    "lines.linewidth": 1,
    "lines.markersize": 5,
    
    # Figure size (in inches)
    #"figure.figsize": (3.5, 2.5),  # Good for single-column in papers
    #"figure.dpi": 300,
    
    # Savefig settings
    "savefig.dpi": 600,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02
})