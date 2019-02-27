
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

input_path = '/mnt/dwf/archive_NOAO_data/data_outputs/2015/01/4hr/g_band/single/c4d_150116_053736_ooi_g_v1/SE_cats/'
#output_path = '/mnt/dwf/archive_NOAO_data/data_outputs/'
#output_cats = '/mnt/dwf/archive_NOAO_data/data_outputs/2015/01/4hr/vizier_cats/'

path_list = glob.glob(input_path)

'''
if not os.path.exists(output_cats):
	os.makedirs(output_cats)
else:
	pass 
'''
field_RA_DEC = '62.5000000 -55.0000000' 

AAVOS = 'II/336'			
#v = Vizier(columns = ['mag', '_RAJ2000', '_DEJ2000'], catalog = AAVOS)
#result = v.query_region(field_RA_DEC, radius=Angle(1.5, "deg"))
result = Vizier.query_region(field_RA_DEC, radius=Angle(1.5, "deg"), catalog=AAVOS)

AAVOS_table = Table()
AAVOS_table['RA'] = result[0]['RAJ2000']
AAVOS_table['DEC'] = result[0]['DEJ2000']
AAVOS_table['g_mag'] = result[0]['g_mag']
AAVOS_table['g_mag_err'] = result[0]['e_g_mag']
AAVOS_table['r_mag'] = result[0]['r_mag']
AAVOS_table['r_mag_err'] = result[0]['e_r_mag']
AAVOS_table['i_mag'] = result[0]['i_mag']
AAVOS_table['i_mag_err'] = result[0]['e_i_mag']


AAVOS_output = 'AAVOS.ascii'
AAVOS_table.write(output_cats + AAVOS_output, format='ascii', overwrite =True) 


USNO_B1 = 'I/284'			
result = Vizier.query_region(field_RA_DEC, radius=Angle(1.5, "deg"), catalog=USNO_B1)

USNO_B1_output = 'USNO_B1.ascii'
result[0].write( output_cats + USNO_B1_output, format='ascii', overwrite =True) 
print(result[0])
