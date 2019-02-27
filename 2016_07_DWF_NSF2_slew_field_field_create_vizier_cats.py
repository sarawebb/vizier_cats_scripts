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

field_RA_DEC = '322.0000000   -66.8000000' 
year = '2016'
month = '07' 
field =  'DWF_NSF2_slew'
 

########################################################################################################

output_path = '/mnt/dwf/archive_NOAO_data/data_outputs/'
output_cats = '/mnt/dwf/archive_NOAO_data/data_outputs/'+year+'/'+month+'/'+field+'/vizier_cats/'



if not os.path.exists(output_cats):
	os.makedirs(output_cats)
else:
	pass 

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

########################################################################################################


USNO_B1 = 'I/284'			
result = Vizier.query_region(field_RA_DEC, radius=Angle(1.5, "deg"), catalog=USNO_B1)


USNO_table = Table()
USNO_table['RA'] = result[0]['RAJ2000']
USNO_table['DEC'] = result[0]['DEJ2000']
USNO_table['R1mag'] = result[0]['R1mag']
USNO_table['R2mag'] = result[0]['R2mag']
USNO_table['Imag'] = result[0]['Imag']

USNO_output = 'USNO.ascii'
USNO_table.write(USNO_output, format='ascii', overwrite= True)

USNO_B1 = 'I/284'			
result = Vizier.query_region(field_RA_DEC, radius=Angle(1.5, "deg"), catalog=USNO_B1)

USNO_B1_output = 'USNO_B1.ascii'
result[0].write(output_cats + USNO_B1_output, format='ascii', overwrite =True) 
########################################################################################################

gaiaDR2 = 'I/345'			
result = Vizier.query_region(field_RA_DEC, radius=Angle(1.5, "deg"), catalog=gaiaDR2)


gaiaDR2_table = Table()
gaiaDR2_table['RA'] = result[0]['RA_ICRS']
gaiaDR2_table['DEC'] = result[0]['DE_ICRS']
gaiaDR2_table['Gmag_Vega'] = result[0]['Gmag']
gaiaDR2_table['Gmag_Vega_err'] = result[0]['e_Gmag']
gaiaDR2_table['Gmag_AB'] = gaiaDR2_table['Gmag_Vega'] -0.08

print(gaiaDR2_table)
gaiaDR2_output = 'gaiaDR2.ascii'
gaiaDR2_table.write(output_cats + gaiaDR2_output, format='ascii', overwrite= True)


