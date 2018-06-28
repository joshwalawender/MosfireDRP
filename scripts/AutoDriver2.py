import os

def nfiles(list_file):
    if os.path.exists(list_file) is False:
        return 0
    else:
        contents = IO.list_file_to_strings(file)
        return len(contents)

# Main Program

nf = True
maskname = 'long2pos_specphot (align)'
band = 'H'
longexp = True
obsfiles = ['Offset_7.0.txt','Offset_-7.0.txt']

if longexp is True:
    wavfiles = obsfiles
elif nfiles('Ne.txt') > 0:
    wavfiles = ['Ne.txt']
elif nfiles('Ar.txt') > 0:
    wavfiles = ['Ar.txt']
else:
    wavfiles = obsfiles

redfiles = [f"eps_{file}.fits" for file in obsfiles]

template = f'''
# If you have questions, please submit a ticket on the github issue page:
# https://github.com/Keck-DataReductionPipelines/MosfireDRP/issues

# User Options:
noninteractiveflag = {str(nf):5s} # Set this flag to autofit wavelength
                           # solution instead of manually verifying the fit.

# Imports
import matplotlib
matplotlib.use('TkAgg') # Force TkAgg backend for interactivity. This is
                        # critical to bypass a bug in the MacOSX backend.
from MOSFIRE import Options, Flats, Wavelength, Background, Rectify, Extract

maskname = '{maskname}'
band = '{band}'
obsfiles = {obsfiles}
wavfiles = {wavfiles}
redfiles = {redfiles}
Wavelength_file = 'lambda_solution_wave_stack_H_m170910_0401-0418.fits'

Flats.handle_flats('Flat.txt', maskname, band, Options.flat)

Wavelength.imcombine(wavfiles, maskname, band, Options.wavelength)
Wavelength.fit_lambda_interactively(maskname, band, obsfiles,
           Options.wavelength, noninteractive=noninteractiveflag)
Wavelength.fit_lambda(maskname, band, obsfiles, obsfiles, Options.wavelength)
Wavelength.apply_lambda_simple(maskname, band, obsfiles, Options.wavelength)

Background.handle_background(obsfiles, Wavelength_file, maskname, band,
           Options.wavelength)

Rectify.handle_rectification(maskname, redfiles, Wavelength_file, band,
        obsfiles, Options.wavelength)

Extract.extract_spectra(maskname, band, width=10,
        interactive=(not noninteractiveflag))
'''
