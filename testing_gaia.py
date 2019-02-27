import numpy as np
import matplotlib as mpl
from astropy.table import Table, Column, join 
#mpl.use('Agg')
import matplotlib.pyplot as plt
from astropy import wcs
from astropy.wcs import WCS
from astropy.io import fits
import sys
import math
import os
import glob
import sys
from sortedcontainers import SortedDict
import datetime as dt
import imageio
import os
from PIL import Image
from matplotlib.colors import LogNorm
from astropy.nddata.utils import Cutout2D
from astropy import units as u
import datetime as dt 
import astropy.units as u
from astroML.crossmatch import crossmatch_angular 
from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats, sigma_clip
from astroquery.vizier import Vizier
from astropy.coordinates import Angle
import glob



######################### INPUT field values ###########################################################

field_RA_DEC = '62.5000000 -55.0000000' 
year = '2015'
month = '01' 
field =  '4hr'
 

########################################################################################################

gaiaDR2 = 'I/345'			
result = Vizier.query_region(field_RA_DEC, radius=Angle(1.5, "deg"), catalog=gaiaDR2)

#gaiaDR2_output = 'USNO_B1.ascii'
#result[0].write(output_cats + gaiaDR2_output, format='ascii', overwrite =True) 

#print(result[0])

gaiaDR2_table = Table()
gaiaDR2_table['RA'] = result[0]['RA_ICRS']
gaiaDR2_table['DEC'] = result[0]['DE_ICRS']
gaiaDR2_table['Gmag_Vega'] = result[0]['Gmag']
gaiaDR2_table['Gmag_Vega_err'] = result[0]['e_Gmag']
gaiaDR2_table['Gmag_AB'] = gaiaDR2_table['Gmag_Vega'] -0.08

print(gaiaDR2_table)
gaiaDR2_output = 'gaiaDR2.ascii'
gaiaDR2_table.write(gaiaDR2_output, format='ascii', overwrite= True)


