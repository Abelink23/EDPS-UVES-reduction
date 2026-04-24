
import os
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.constants import c

from scipy.signal import convolve

class ob():
    def __init__(self, file, orig='pipeline', bar_corr=False, to_vac=False,
                cut_edges=False, clean=False, tare=False):
        '''
        Function to extract the wavelength and flux from a UVES spectrum FITS file.
        It creates an object and optionally applies:
            - Barycentric correction
            - Air to vacuum conversion (Morton, Griesen, or IAU methods)
            - Cutting of the edges of the spectrum
            - Cleaning of spikes in the spectrum by different approaches.
            - Division of the flux by its median value.
            - Export of the spectrum in an ascii or FITS file.

            Parameters
            ----------
            file : str
                The path to the file to be extracted.
            orig : str, optional
                The origin of the file to be extracted. It can be:
                - 'pipeline' for files from the pipeline,
                - 'reduced' for master reduced files,
                - 'ascii' for ascii files.
        '''
        self.file = file

        if orig == 'reduced' or orig == 'master' or \
            file.split('/')[-1].upper().startswith('MASTER_'):
            self.extract_from_master()
        elif orig == 'pipeline':
            self.extract_from_pipeline()
        elif orig == 'ascii':
            self.extract_from_ascii()
        else:
            print(f"Origin {orig} not recognized. No extraction done.")
            return

        self.bar_corr_done = False
        if bar_corr:
            self.barycentric_correction()

        self.air_to_vac_done = False
        if to_vac in ['Morton', 'Griesen', 'IAU']:
            self.air_to_vacuum(method=to_vac)

        self.cut_edges_done = False
        if cut_edges:
            self.cut_edges()

        self.clean_spikes_done = False
        if clean:
            self.clean_spikes(method=clean)

        self.tare_done = False
        if tare:
            self.tare()

        self.clean_negative_flux()

    def extract_from_pipeline(self):
        ''' Extract wavelength and flux of spectra from the pipeline output files.
            #! CSYER1 = 0,000115 -> is this the offset I have to apply to match with the EDPS plot?
        '''
        sp = fits.open(self.file)
        self.header = sp[0].header

        if "ADP." in self.header['ARCFILE'].upper():
            self.wave = sp[1].data[0][0]
            self.flux = sp[1].data[0][1]
        else:
            self.dlam = self.header['CDELT1']   # Step of increase in wavelength
            lam0 = self.header['CRVAL1']        # Get the wavelength of the first pixel
            pix0 = self.header['CRPIX1']        # Reference pixel (generally 1)
            spec_length = self.header['NAXIS1'] # Length of the spectrum
            self.wave = lam0 + self.dlam*(np.arange(spec_length) - pix0 + 1)
            self.flux = sp[0].data

        # Get Dichroic info
        self.dich = self.header['HIERARCH ESO TPL NAME'].split(' ')[0].upper()

        # Get central wavelength, CCD, and slit width info
        self.ccd = self.header['HIERARCH ESO PRO CATG'].split('_')[-1]
        if self.ccd == 'REDL' or self.ccd == 'REDU' or self.ccd == 'RED':
            self.cwlen = int(self.header['HIERARCH ESO INS GRAT2 WLEN'])
            self.slit_width = self.header['HIERARCH ESO INS SLIT3 WID']
        elif self.ccd == 'BLUE' or self.ccd == 'BLU':
            self.cwlen = int(self.header['HIERARCH ESO INS GRAT1 WLEN'])
            self.slit_width = self.header['HIERARCH ESO INS SLIT2 WID']

        print(f'{self.header["ARCFILE"].replace(".fits","")} | {self.ccd} | {self.dich} | {self.cwlen} | {self.slit_width}"')
        sp.close()

    def extract_from_master(self):
        ''' Extract wavelength and flux of spectra from the master reduced files.
        '''
        sp = fits.open(self.file)
        self.header = sp[0].header
        self.wave = sp[1].data
        self.flux = sp[2].data
        self.flux_error = sp[3].data
        self.flux_cal = sp[4].data
        self.flux_cal_error = sp[5].data
        try:
            self.norm_flux = sp[6].data
        except:
            None
        sp.close()

        # Assume that the corrections have already been applied in the master reduced files
        self.bar_corr_done = True
        self.air_to_vac_done = True
        self.cut_edges_done = True
        self.clean_spikes_done = True

        # Get Dichroic info
        self.dich = self.header['HIERARCH ESO TPL NAME'].split(' ')[0].upper()

        # Get central wavelength, CCD, and slit width info
        self.ccd = self.header['HIERARCH ESO PRO CATG'].split('_')[-1]
        if self.ccd == 'REDL' or self.ccd == 'REDU':
            self.cwlen = int(self.header['HIERARCH ESO INS GRAT2 WLEN'])
            self.slit_width = self.header['HIERARCH ESO INS SLIT3 WID']
        elif self.ccd == 'BLUE':
            self.cwlen = int(self.header['HIERARCH ESO INS GRAT1 WLEN'])
            self.slit_width = self.header['HIERARCH ESO INS SLIT2 WID']

        # Calculate dlam as the median of the difference between consecutive wavelengths
        self.dlam = np.median(np.diff(self.wave))

        # Print summary of the spectrum
        print(f'{self.header["ARCFILE"].replace(".fits","")} | {self.ccd} | {self.dich} | {self.cwlen} | {self.slit_width}"')

    def extract_from_ascii(self):
        ''' Extract wavelength and flux from an ascii file.
        '''
        self.header = None

        data = np.loadtxt(self.file, skiprows=1)
        self.wave = data[:, 0]
        self.flux = data[:, 1]

        self.dlam = np.median(np.diff(self.wave))

        self.dich = None
        self.cwlen = None
        self.ccd = None

    def clean_negative_flux(self):
        ''' Remove negative flux values from the spectrum.
        '''
        if np.any(self.flux < 0):
            print(f"\033[93mWarning: {np.sum(self.flux < 0)} negative flux values found. They will be set to zero.\033[0m")
        self.flux[self.flux < 0] = 0.0

    def barycentric_correction(self):
        ''' Correct from barycentric velocity
        '''
        if 'BARYCORR' in self.header:
            vbar = self.header['BARYCORR']
            self.wave = self.wave*(1 + 1000*vbar/c.value)
            self.bar_corr_done = True
            print(f"Barycentric correction applied: velocity = {vbar:.2f} m/s.")
        else:
            print("Warning: BARYCORR keyword not found in header. No correction applied.")

    def air_to_vacuum(self, method='Morton'):
        ''' Air to vacuum correction applied to the wavelength using different formulas.
        '''
        if self.air_to_vac_done in ['Morton', 'Griesen', 'IAU']:
            print(f"Air to vacuum conversion already applied. No correction done.")
            return

        # Use Morton (1991, ApJS, 77, 119) - IAU standard
        if method == 'Morton':
            self.wave = self.wave * (1 + 2.735182e-4 + 131.4182/self.wave**2 + 2.76249e8/self.wave**4)
        # Use Griesen (2006, A&A 446, 747) - used by specutil
        if method == 'Griesen':
            self.wave = self.wave * (1 + 1e-6 * (287.6155 + 1.62887e5/self.wave**2 + 0.01360/self.wave**4))
        # Use as actually implemented by wcslib
        if method == 'IAU':
            s = (1.0 / (self.wave/1e10))**2 # change wave to meters
            n = 2.554e8 / (0.41e14 - s)
            n += 294.981e8 / (1.46e14 - s)
            n += 1.000064328
            self.wave = self.wave * n
        self.air_to_vac_done = method
        print(f"Air to vacuum conversion applied using {method} method.")

    def cut_edges(self):
        ''' Make cuts at the edges of the spectra to remove unwanted regions
        '''
        # Make cuts at the edges of the spectra to remove unwanted regions
        if self.ccd == 'BLUE' and self.cwlen == 346 and self.dich == 'DIC1':
            self.flux = self.flux[(self.wave > 3060) & (self.wave < 3860.5)]
            self.wave = self.wave[(self.wave > 3060) & (self.wave < 3860.5)]
        if self.ccd == 'BLUE' and self.cwlen == 437 and self.dich == 'DIC2':
            self.flux = self.flux[(self.wave > 3761.5) & (self.wave < 4979)]
            self.wave = self.wave[(self.wave > 3761.5) & (self.wave < 4979)]
        if self.ccd == 'REDL' and self.cwlen == 580 and self.dich == 'DIC1':
            self.flux = self.flux[(self.wave > 4792) & (self.wave < 5751)]
            self.wave = self.wave[(self.wave > 4792) & (self.wave < 5751)]
        if self.ccd == 'REDU' and self.cwlen == 580 and self.dich == 'DIC1':
            self.flux = self.flux[(self.wave > 5853) & (self.wave < 6791)]
            self.wave = self.wave[(self.wave > 5853) & (self.wave < 6791)]
        if self.ccd == 'REDL' and self.cwlen == 860 and self.dich == 'DIC2':
            self.flux = self.flux[(self.wave > 6712) & (self.wave < 8510)]
            self.wave = self.wave[(self.wave > 6712) & (self.wave < 8510)]
        if self.ccd == 'REDU' and self.cwlen == 860 and self.dich == 'DIC2':
            self.flux = self.flux[(self.wave > 8674) & (self.wave < 10411)]
            self.wave = self.wave[(self.wave > 8674) & (self.wave < 10411)]
        self.cut_edges_done = True
        print("Edges cut applied.")

    def clean_spikes(self, method, wl_split=100, dmin=0.05, protect_em_lines=False,
                    zs_cut=5, ker_sig=2, ker_iter=3, sig_g=None):
        '''
        Function to remove spikes from the spectra by different approaches.

        Parameters
        ----------
        method : str, optional
            Method for the removal strategy: zscore (default) or kernel.
        dmin : int/float, optional
            Minimum distance between flux and cleaned flux to consider for replacement.
            Default is 0.05.
        zs_cut : int/float, optional
            In 'zscore', threshold value for finding spikes. Default is 4.
            Tip: For noisy spectra, increase this value up to 7-9.
        wl_split : int/float, optional
            In 'kernel' method, wavelength size used to split the spectrum before
            applying the removal. Default is 100.
        ker_sig : float, optional
            In 'kernel' method, sigma clipping value used to remove spikes. Default is 2.
        ker_iter : int, optional
            In 'kernel' method, number of iterations used to remove spikes. Default is 3.
        ker_sig_g : float, optional
            In 'kernel' method, sigma of the gaussian function used to construct the kernel.
            Default is the theoretical sigma based on wavelength and resolution.
        protect_em_lines : boolean, optional
            If True, some emission lines will be masked from spike removal.
            Default is False.
        '''

        if method == 'zscore':
            # www.towardsdatascience.com/removing-spikes-from-raman-spectra-8a9fdda0ac22

            # First we calculate (nabla)x(i):
            delta_flux = np.diff(self.flux)
            median_int = np.median(delta_flux)
            mad_int = np.median([np.abs(delta_flux - median_int)])
            modified_z_scores = 0.6745 * (delta_flux - median_int) / mad_int
            # The multiplier 0.6745 is the 0.75th quartile of the standard normal
            # distribution, to which the median absolute deviation converges to.
            modified_z_scores =  np.concatenate(([0], np.abs(modified_z_scores)))

            self.flux_clean = np.where(modified_z_scores > zs_cut, np.nan, self.flux)

        elif method == 'kernel':
            # Range in angstroms in which the spectrum will be initially splitted
            wl_split = 100

            wl_range = self.wave[-1]-self.wave[0]
            if wl_range < wl_split:
                wl_split = wl_range

            if 1 < wl_range/wl_split <= 2:
                wl_split = wl_range/2 + 1
            elif wl_range/wl_split > 2:
                wl_split += (wl_range % wl_split) / int(wl_range/wl_split)

            n = int(wl_range/wl_split)

            flux_clean = []
            for i in range(n):

                mask = (self.wave >= self.wave[0]+i*wl_split) & (self.wave < self.wave[0]+(i+1)*wl_split)
                if i == range(n)[-1]:
                    mask[-1] = True # To catch the last flux value from not using <= above

                resolution = 100000 #! average resolution of UVES

                if ker_sig_g is None:
                    lambda0 = np.mean(self.wave[mask])
                    ker_sig_g = lambda0/(2.35482*resolution)
                else:
                    ker_sig_g = float(ker_sig_g)

                x = np.arange(-5*ker_sig_g, 5*ker_sig_g+self.dlam, self.dlam)
                gauss = np.exp(-0.5*(x/ker_sig_g)**2)
                kernel = gauss/np.trapz(gauss)

                convoluted = 1 + convolve(self.flux[mask] - 1, kernel, mode='same')

                flux_norm = self.flux[mask]/convoluted

                for i in range(ker_iter):
                    std = np.nanstd(flux_norm)
                    flux_norm = np.where(abs(flux_norm - 1) > ker_sig*std, np.nan, flux_norm)

                self.flux_clean=np.concatenate([flux_clean,
                                np.where(np.isnan(flux_norm), np.nan, self.flux[mask])])

        else:
            print(f"Method {method} not recognized. No cleaning applied.")
            return

        nans = np.isnan(self.flux_clean)
        x = lambda z: z.nonzero()[0]
        self.flux_clean[nans] = np.interp(x(nans), x(~nans), self.flux_clean[~nans])

        # Recover the original flux in the regions of the spectrum where telluric lines are
        self.flux_clean = np.where(
              ((self.wave > 3302.0) & (self.wave < 3304.0))
            | ((self.wave > 3933.0) & (self.wave < 3935.0)) # Is this telluric?
            | ((self.wave > 3967.0) & (self.wave < 3969.0)) # Is this telluric?
            | ((self.wave > 5885.0) & (self.wave < 6000.0))
            | ((self.wave > 6275.0) & (self.wave < 6330.0))
            | ((self.wave > 6440.0) & (self.wave < 6590.0))
            | ((self.wave > 6865.0) & (self.wave < 7400.0))
            | ((self.wave > 7585.0) & (self.wave < 7715.0))
            | ((self.wave > 7865.0) & (self.wave < 8380.0))
            | ((self.wave > 8915.0) & (self.wave < 9900.0)),
            self.flux, self.flux_clean)

        # Recover the original flux in the regions of the spectrum where problematic lines are
        self.flux_clean = np.where(
                ((self.wave > 3383.0) & (self.wave < 3384.0))
            |   ((self.wave > 3968.0) & (self.wave < 3970.0)),
                self.flux, self.flux_clean)

        # Recover the original flux in the regions of the spectrum where emission lines are
        if protect_em_lines == True:
            for wl_em in [3967.79,4958.911,5006.843,6300.304,6548.04,6583.46,6716.44,6730.82]:
                mask = (self.wave > wl_em-0.8) & (self.wave < wl_em+0.8)
                self.flux_clean[mask] = self.flux[mask]

        # Recover the original flux if the difference is lower than dmin value
        self.flux_clean = np.where(np.abs(self.flux - self.flux_clean) > dmin, self.flux_clean, self.flux)
        self.flux = self.flux_clean
        print(f"Spikes cleaned using {method} method.")

    def normalize(self, step=None, append_master=False):
        '''Normalize a complete spectrum by regions.

        It adds to the class the following attributes:
            self.norm_flux - Normalized flux array of the spectrum.
            self.continuum - Continuum fit of the spectrum.
            self.norm_mask - Mask of the points used for the normalization (1 good, 0 bad).
            self.lam - Wavelengths of the strong lines used for normalization.

        Parameters
        ----------
        step : int/float, optional
            Step in angstroms to select the regions for normalization.
            If None, it will be selected depending on the setup.
            Default is None.
        '''

        lam0 = np.min(self.wave)
        w = np.array(self.wave.astype(float))
        f = np.array(self.flux.astype(float))

        # Selection of the individual regions for normalization depending on the setup
        if step is not None:
            lam = np.arange(w.min() + step, w.max(), step)
        elif self.ccd == 'BLUE' and self.dich == 'DIC2' and self.cwlen == 437:
            print("Normalization regions selected for BLUE DIC2 437 setup.")
            lam = np.array([3815,3920,4060,4150,4290,4605,4830,4950], float)
        else:
            lam = np.array([3815,3920,4060,4150,4290,4605,4830,4950,5070,5320,5500,5620,5750,5840,
                            5980,6200,6460,6700,7100,7620,8000,8400,8710,8930], float)

        # Take regions within the range of the spectrum
        mask = (lam >= w.min()) & (lam <= w.max())
        lam = np.concatenate([lam[mask], [w.max()]]).astype(float)

        ws_list, fs_list, fn_list, yn_list, fk_list, snr_list = [], [], [], [], [], []
        # First region
        ws0, fs0, fn0, yn0, fk0, snr0, wlim1, wlim2, lam1 = normalize_slice(w, f, lam0, lam[0])
        ws0, fs0, fn0, yn0, fk0, snr0, wlim1, wlim2, lam1 = normalize_slice(
            w, f, lam0, lam[0], wlim1=wlim1, wlim2=wlim2, lam1=lam1, iter=1)
        ws_list.append(ws0); fs_list.append(fs0); fn_list.append(fn0); yn_list.append(yn0); fk_list.append(fk0)
        snr_list.append(snr0)

        # Rest of the regions
        for i in range(1, len(lam)):
            ws0, fs0, fn0, yn0, fk0, snr0, wlim1, wlim2, lam1 = normalize_slice(
                w, f, lam[i-1], lam[i], wlim1=wlim1, wlim2=wlim2, lam1=lam1)
            ws0, fs0, fn0, yn0, fk0, snr0, wlim1, wlim2, lam1 = normalize_slice(
                w, f, lam[i-1], lam[i], wlim1=wlim1, wlim2=wlim2, lam1=lam1, iter=1)
            ws_list.append(ws0); fs_list.append(fs0); fn_list.append(fn0); yn_list.append(yn0); fk_list.append(fk0)
            snr_list.append(snr0)

        ws = np.concatenate(ws_list)
        fs = np.concatenate(fs_list)
        fn = np.concatenate(fn_list)
        yn = np.concatenate(yn_list)
        fk = np.concatenate(fk_list)
        snr = np.array(snr_list)

        if ws is None or fn is None or len(ws) == 0 or len(fn) == 0:
            ft = np.zeros_like(w)
            ynt = np.zeros_like(w)
        else:
            ft = np.interp(w, ws, fn)
            ynt = np.interp(w, ws, yn)
        # global continuum interpolated
        # global mask passing from ws to w, trusting ws with w
        fkt = np.zeros_like(w, float)
        fkt[np.searchsorted(w, ws)] = fk

        self.norm_flux = ft
        self.continuum = ynt
        self.norm_mask = fkt
        self.snr = np.median(snr)

        # Add the normalization after flux_cal_error as extra extension to the master reduced FITS file if append_master is True. In this way the FITS will have the same structure as the master reduced files, and the normalized flux will be available for future use without needing to re-apply the normalization.
        if append_master:
            filename = self.file.replace('.fits', '_norm.fits')
            hdul = fits.open(self.file)
            hdul.append(fits.ImageHDU(data=self.norm_flux, header=self.header))
            if os.path.exists(filename):
                overwrite = input(f"File {filename.split('/')[-1]} already exists, overwrite? (y/n) ")
                if overwrite.lower() != 'y':
                    print("File not overwritten.")
                    return
            hdul.writeto(filename, overwrite=True)
            print(f"Normalized spectrum added as a new extension to {filename.split('/')[-1]}.")


    def tare(self):
        ''' Divide the flux by its median value
        '''
        if not self.tare_done:
            self.flux_median = np.median(self.flux)
            self.flux = self.flux / self.flux_median
            self.tare_done = True
            print("Flux/median applied: flux divided by its median value.")
        else:
            self.flux = self.flux * self.flux_median
            self.tare_done = False
            print("Flux/median reverted: flux multiplied by its median value.")

    def plot(self, y='flux', delta_y=0, alpha=1.0, label=None):
        '''
        Function to plot the spectrum.

        Parameters
        ----------
        y : str, optional
            'flux' for the original flux.
            'flux_cal' for the flux calibrated flux.
            'norm_flux' for the normalized flux.
            Default is 'flux'.
        delta_y : float, optional
            Value added to the flux for visualization purposes. Default is 0.
        alpha : float, optional
            Transparency of the plot. Default is 1.0.
        label : str, optional
            Label of the plot. Default is None.
        '''
        if not hasattr(self, y):
            raise ValueError(f"Attribute '{y}' not found in the object.")
            return

        if y == 'flux':
            flux = self.flux
        elif y == 'norm_flux':
            flux = self.norm_flux
        elif y == 'flux_cal':
            flux = self.flux_cal

        plt.plot(self.wave, flux + delta_y, alpha=alpha, lw=0.7,
                label=label or self.file.split('/')[-1])
        plt.xlabel('Wavelength (Angstrom)')
        if 'fluxcal_' in self.file.lower() or y == 'flux_cal':
            plt.ylabel('Flux (erg/s/cm**2)')
        elif 'spectrum_' in self.file.lower() or y == 'flux':
            plt.ylabel('Flux (ADU)')
        elif 'error_red' in self.file.lower():
            plt.ylabel('Flux error (ADU)')
        elif y == 'norm_flux':
            plt.ylabel('Normalized flux')
        else:
            plt.ylabel('Flux')

        plt.gcf().set_size_inches(12, 6)
        plt.tight_layout()
        plt.show(block=False)

    def plot_orders(self):
        '''
        Function to plot the orders of a resampled_ff_science_*.fits file.
        '''
        if len(range(self.flux.shape[0])) == 0:
            print("No orders found in the spectrum. Is this a resampled_ff_science_*.fits file?")
            return

        for i in range(self.flux.shape[0]):
            flux_orders = self.flux[i, :]
            # filter flux with zeroes
            mask = (flux_orders > 0) & (self.wave > 1)
            wave_i = self.wave[mask]
            flux_i = flux_orders[mask]
            plt.plot(wave_i, flux_i + i*1e3, alpha=0.7, lw=0.7)
            plt.text(max(wave_i), flux_i[-1] + i*1e3, f'Order {i+1}', color='black')

        # remove labels of y-axis
        plt.yticks([])
        plt.xlabel('Angstroms')
        plt.ylabel('Flux')
        plt.tight_layout()
        plt.show(block=False)

    def export_fits(self, tail=''):
        '''Function to export the current status of the file to a FITS file.

        Parameters
        ----------
        tail : str, optional
            Tail of the file added before the extension for its identification.
        '''
        if tail and not tail.startswith('_'):
            tail = '_' + tail

        filename = self.file.replace('.fits', tail + '.fits')
        hdu = fits.PrimaryHDU(data=self.flux, header=self.header)
        if os.path.exists(filename):
            overwrite = input(f"File {filename} already exists. Do you want to overwrite it? (y/n) ")
            if overwrite.lower() != 'y':
                print("File not overwritten.")
                return
        hdu.writeto(filename, overwrite=True)
        plt.close('all')

    def export_ascii(self, tail=''):
        '''
        Function to export the current wavelength and flux of the spectrum in the class
        into an ascii file.

            Parameters
            ----------
            tail : str, optional
                Tail of the file added before the extension for its identification.
                Default is ''.
        '''

        filename = self.file.replace('.fits', '')
        new_name = filename + tail + '.ascii'
        if os.path.exists(new_name):
            overwrite = input(f"File {new_name} already exists. Do you want to overwrite it? (y/n) ")
            if overwrite.lower() != 'y':
                print("File not overwritten.")
                return

        if hasattr(self, 'flux_error') and self.flux_error is not None:
            np.savetxt(new_name, np.c_[self.wave,self.flux,self.flux_error], fmt=['%.4f','%.6e','%.6e'],
            header='wave      flux   flux_error', comments='')
        else:
            np.savetxt(new_name, np.c_[self.wave,self.flux], fmt=['%.4f','%.6e'],
            header='wave      flux', comments='')

def normalize_slice(w, f, w0, w1, wlim1=None, wlim2=None, lam1=None, iter=None):

    '''
    Function to normalize the spectrum by fitting a polynomial to the continuum
    regions. It removes the strong lines from the spectrum and iteratively fits a
    polynomial to the remaining points, applying sigma-clipping to remove points
    affected by weak lines or noise.

    Parameters
    ----------
    w : array-like
        Wavelength array of the spectrum.
    f : array-like
        Flux array of the spectrum.
    w0 : float
        Lower limit of the wavelength region to be normalized.
    w1 : float
        Upper limit of the wavelength region to be normalized.
    wlim1 : array-like, optional
        Lower limits of the wavelength regions to be normalized.
    wlim2 : array-like, optional
        Upper limits of the wavelength regions to be normalized.
    lam1 : array-like, optional
        Wavelengths of the strong lines to be removed.
    iter : int, optional
        If not None, it indicates that this is an iterative call to the function, and the strong lines should not be added to lam1 again.

    Returns
    -------
    ws0 : array-like
        Wavelength array of the normalized spectrum in the given region.
    fs0 : array-like
        Flux array of the normalized spectrum in the given region.
    fn : array-like
        Normalized flux array of the spectrum in the given region.
    yn : array-like
        Continuum fit of the spectrum in the given region.
    fk : array-like
        Mask of the points used for the normalization (1 for good points, 0 for bad points) in the given region.
    wlim1 : array-like
        Updated lower limits of the wavelength regions to be normalized.
    wlim2 : array-like
        Updated upper limits of the wavelength regions to be normalized.
    lam1 : array-like
        Updated wavelengths of the strong lines to be removed.
    '''

    # Strong lines
    lamH = np.array([3712,3722,3735,3750,3771,3797,3835,3890,3970,
                     4102,4340,4860,6563,8438,8465,8500,8545,8600,8665,8748,8860], dtype=float)
    lamHeI  = np.array([3820,3926,4009.,4026.,4143.,4387.,4471.5,4713,4922,5875,6678])
    lamHeII = np.array([4200,4542,4686,5411])
    lamISM  = np.array([4430])

    lam0  = np.concatenate([lamH, lamHeI, lamHeII, lamISM])
    dlam0 = np.concatenate([lamH*0+5., lamHeI*0+2., lamHeII*0+2., lamISM*0+5.])
    tt = np.argsort(lam0)
    lam0 = lam0[tt]
    dlam0 = dlam0[tt]
    snr = np.nan

    # Cut region
    mask = (w > w0) & (w <= w1) & (f > 0)
    if np.sum(mask) < 50:
        print(f"Warning: Region {round(w0, 2)}-{round(w1, 2)} has too few points to be normalized.")
        # Return neutral arrays
        ws0 = w[mask] if mask.any() else np.array([])
        fs0 = f[mask] if mask.any() else np.array([])
        return ws0, fs0, np.ones_like(ws0), np.zeros_like(ws0), np.zeros_like(ws0), snr, wlim1, wlim2, lam1

    ws0 = w[mask]
    fs0 = f[mask]
    ws00 = ws0.copy()
    fs00 = fs0.copy()

    # Keywords initialization
    lam1 = [] if lam1 is None else list(lam1)
    if wlim1 is None or not isinstance(wlim1, np.ndarray) or len(wlim1) != len(lam0):
        wlim1 = np.zeros(len(lam0))
    if wlim2 is None or not isinstance(wlim2, np.ndarray) or len(wlim2) != len(lam0):
        wlim2 = np.zeros(len(lam0))

    # Calculates the FWHM of the line as average size of the line wings
    # The * 0.9 factor reduces that width by 10%.
    if iter is not None:
        dlam0 = 0.5*(np.abs(wlim2)+np.abs(wlim1))*0.9

    # Remove strong lines
    for j in range(len(lam0)):
        if ws0.size == 0:
            break
        if (ws0.min()-10 <= lam0[j]) and (ws0.max()+10 >= lam0[j]):
            good = (ws0 <= lam0[j]-dlam0[j]) | (ws0 >= lam0[j]+dlam0[j])
            ws0 = ws0[good]
            fs0 = fs0[good]
            if iter is None:
                lam1.append(lam0[j])

    ws = ws0
    fs = fs0

    # Initial adjustment of the continuum. If there are not enough points, return the original spectrum.
    ord = 1
    if len(ws) < ord+1:
        # Not enough points to fit a polynomial of the given order.
        print(f"Warning: not enough points to fit a polynomial of order {ord}")
        fn = np.ones_like(ws00)
        yn = np.zeros_like(ws00)
        fk = np.zeros_like(ws00)
        if np.shape(fn)!=np.shape(yn):
            print('fn return de normaliza region sin puntos: ', np.shape(fn))
            print('yn return de normaliza region sin puntos: ', np.shape(yn))
        return ws00, fs00, fn, yn, fk, snr, wlim1, wlim2, np.array(lam1)

    # Note: Simón-Díaz used to add an extra smoothing here (not implemented)
    dat = np.polyfit(ws, fs, ord)
    ys  = np.polyval(dat, ws)

    fn = fs / ys
    fk = np.ones_like(fn)
    mfn = np.mean(fn)
    sfn = np.std(fn)

    # Iterative sigma-clipping
    error0, error = 40., 20.
    value, step = 2., 0
    while abs(error-error0) >= 0.05 and ws.size > 0:
        for i in range(len(ws)):
            if fn[i] < mfn-sfn or fn[i] > mfn+value*sfn:
                bad = (fn < mfn - sfn) | (fn > mfn + value * sfn)
                fn[bad] = mfn + np.random.randn(np.sum(bad)) * sfn
                fk[bad] = 0.

        nn = np.where(fk > 0)[0]
        error0 = error
        if len(nn) < 2:
            break
        dat2 = np.polyfit(ws[nn], fn[nn], 2)
        ys2  = np.polyval(dat2, ws)
        error2 = np.std(fn - ys2)

        dat1 = np.polyfit(ws[nn], fn[nn], 1)
        ys1  = np.polyval(dat1, ws)
        error1 = np.std(fn - ys1)

        if error2 <= error1:
            ys = ys2
            error = abs(error2)
        else:
            ys = ys1
            error = abs(error1)

        if w0 >= 8300:
            dat = np.polyfit(ws[nn], fn[nn], 1)
            ys  = np.polyval(dat, ws)
            error = np.std(fn - ys)

        fn = fn / ys
        mfn = np.median(fn)
        sfn = np.std(fn)
        step += 1

    # Final normalization
    nn = np.where(fk > 0)[0]
    wp = ws[nn] if nn.size > 0 else np.array([])
    fp = fs[nn] if nn.size > 0 else np.array([])

    ord = 2 if w0 < 8300 else 1
    for _ in range(3):
        if wp.size < ord+1:
            break
        dat = np.polyfit(wp, fp, ord)
        yfit = np.polyval(dat, wp) / fp
        keep = np.abs(yfit - np.median(yfit)) <= 2*np.std(yfit)
        wp = wp[keep]
        fp = fp[keep]

    dat = np.polyfit(wp, fp, ord) if wp.size >= ord+1 else np.array([1.])
    yn     = np.polyval(dat, ws00) if wp.size >= ord+1 else np.ones_like(ws00)
    ynoise = np.polyval(dat, ws0)   if wp.size >= ord+1 else np.ones_like(ws0)

    fn = fs00 / yn # This is the key to why fk has fewer points. fn is redefined here, but fk is not.

    # Calculate the SNR - NOT ENABLED
    fnoise = fs0 / ynoise
    snr = calculate_snr(fnoise) if ws0.size > 0 else 0

    # Save the limits of the§ strong lines for the next iteration.
    if iter is None and len(lam1) > 0:
        lam1 = np.array(lam1)
        for line in lam1:
            tt = np.where(lam0 == line)[0]
            if wp.size > 0:
                left = np.where(wp <= line)[0]
                wlim1[tt] = wp[left.max()] - line if left.size > 0 else wp.min() - line
                right = np.where(wp >= line)[0]
                wlim2[tt] = wp[right.min()] - line if right.size > 0 else wp.max() - line
            else:
                wlim1[tt] = 0
                wlim2[tt] = 0

    # Fix small mismatch as fk was defined in ws, and everything else in ws00. It was coming out smaller
    fkt = np.zeros_like(ws00)       #fkt = np.zeros_like(ws00, dtype=float)
    idx = np.searchsorted(ws00, ws) #idx = np.searchsorted(ws00, ws)
    # ensure that we only assign valid indices: valid = (idx >= 0) & (idx < len(ws00)) [possible alternative]
    fkt[idx] = fk                   #fkt[idx[valid]] = fk[valid]

    return ws00, fs00, fn, yn, fkt, snr, wlim1, wlim2, lam1

def calculate_snr(flux):
    # Calculate the SNR of the spectrum with an iterative sigma clipping method
    me0 = np.nanmean(flux)
    st0 = np.nanstd(flux)
    snr0 = me0 / st0
    snr = snr0 + 100
    while abs(snr0 - snr) >= 1.0:
        snr = snr0
        tt = np.where(np.abs(flux - me0) <= 2.0 * st0)
        me0 = np.nanmean(flux[tt])
        st0 = np.nanstd(flux[tt])
        snr0 = me0 / st0

    return int(snr0)

def make_master(folder):
    ''' Function to create a master FITS and ASCII files from:
        - red_science_*.fits (contains the spectrum)
        - error_red_science_*.fits (contains the error of the spectrum)
        - fluxcal_science_*.fits (contains the flux-calibrated spectrum)
        - fluxcal_error_science_*.fits (contains the error of the flux-calibrated spectrum)
        - parameters.rc file (contains the non-default parameters used in the reduction)

    Note: the header of the FITS is taken from fluxcal_science file.
    '''

    # find each file in the folder
    ccd = None
    for root, dirs, files in os.walk(folder):
        if 'red_science_blue.fits' in files:
            root_dir = root
            ccd = 'BLUE'
            break
        elif 'red_science_redl.fits' in files:
            ccd = 'RED'
            root_dir = root
            break

    parameters = os.path.join(root_dir, 'parameters.rc')
    parameters = np.loadtxt(parameters, dtype=str, delimiter='=', comments='#')
    parameters_dict = {param[0].strip(): param[1].strip() for param in parameters}

    # load the components based on the CCD
    parts = []
    if ccd == 'BLUE':
        arms = [('blue', '')]
    elif ccd == 'RED':
        arms = [('redl', '_L'), ('redu', '_U')]

    for arm, ext_suffix in arms:
        r = ob(os.path.join(root_dir, f'red_science_{arm}.fits'),
                orig='pipeline', bar_corr=True, to_vac='Morton', cut_edges=True)
        e_r = ob(os.path.join(root_dir, f'error_red_science_{arm}.fits'),
                orig='pipeline', bar_corr=True, to_vac='Morton', cut_edges=True)
        fcal = ob(os.path.join(root_dir, f'fluxcal_science_{arm}.fits'),
                orig='pipeline', bar_corr=True, to_vac='Morton', cut_edges=True)
        fcal_e = ob(os.path.join(root_dir, f'fluxcal_error_science_{arm}.fits'),
                orig='pipeline', bar_corr=True, to_vac='Morton', cut_edges=True)
        parts.append((arm, ext_suffix, r, e_r, fcal, fcal_e))

    for filename_suffix, ext_suffix, r, e_r, fcal, fcal_e in parts:
        # Create the master FITS file
        # add the parameters.rc file as a header extension to the master FITS file
        for key, value in parameters_dict.items():
            fcal.header[key] = value

        hdu = fits.PrimaryHDU(header=fcal.header)
        hdu_wave = fits.ImageHDU(data=r.wave, name=f'WAVE{ext_suffix}')
        hdu_flux = fits.ImageHDU(data=r.flux, name=f'FLUX{ext_suffix}')
        hdu_error = fits.ImageHDU(data=e_r.flux, name=f'ERROR{ext_suffix}')
        hdu_fluxcal = fits.ImageHDU(data=fcal.flux, name=f'FLUXCAL{ext_suffix}')
        hdu_fluxcal_error = fits.ImageHDU(data=fcal_e.flux, name=f'FLUXCAL_ERROR{ext_suffix}')

        hdul = fits.HDUList([hdu, hdu_wave, hdu_flux, hdu_error, hdu_fluxcal, hdu_fluxcal_error])
        arcfile = r.header['ARCFILE'].replace('.fits', '')
        hdul.writeto(os.path.join(folder, f'MASTER_{filename_suffix.upper()}_{arcfile}.fits'), overwrite=True)

        # Also create an ASCII file
        np.savetxt(os.path.join(folder, f'MASTER_{filename_suffix.upper()}_{arcfile}.ascii'),
                                np.c_[r.wave, r.flux, e_r.flux, fcal.flux, fcal_e.flux],
                                fmt=['%.4f', '%.6e', '%.6e', '%.6e', '%.6e'],
        header='wave      flux         flux_error   fluxcal      fluxcal_error', comments='')


def plt_all_spec(file_type, tare=False, alpha=1.0, diff=False):
    folder = '/Users/adeburgo/Documents/pipelines/EDPS_data/UVES/object'
    all_files = []
    for root, dirs, files in os.walk(folder):
        for file in files:
            if file == file_type:
                path = os.path.join(root, file)
                all_files.append(path)
                date = os.path.getctime(path)
                date = datetime.fromtimestamp(date).strftime('%Y-%m-%d %H:%M:%S')
                print(f"Plotting {file} from {date}")
                spec = ob(path, orig='pipeline', bar_corr=True, to_vac='Morton', cut_edges=True, tare=tare)
                spec.plot(alpha=alpha, label=f"{date}")

    if diff and len(all_files) == 2:
        plt.figure(figsize=(12, 6))
        plt_diff(all_files[0], all_files[1], tare=tare, alpha=alpha)

    plt.legend()

def plt_diff(fits_file1, fits_file2, tare=False, alpha=1.0):
    spec1 = ob(fits_file1, orig='pipeline', bar_corr=True, tare=tare)
    wave1, flux1 = spec1.wave, spec1.flux
    spec2 = ob(fits_file2, orig='pipeline', bar_corr=True, tare=tare)
    wave2, flux2 = spec2.wave, spec2.flux

    # limit to common wavelength range
    wave_min = max(wave1.min(), wave2.min())
    wave_max = min(wave1.max(), wave2.max())
    mask1 = (wave1 >= wave_min) & (wave1 <= wave_max)
    mask2 = (wave2 >= wave_min) & (wave2 <= wave_max)
    wave = wave1[mask1]

    # interpolate flux2 to wave1
    flux2_interp = np.interp(wave, wave2[mask2], flux2[mask2])
    flux_diff = flux1[mask1]/flux2_interp

    plt.plot(wave, flux_diff, alpha=alpha, lw=0.7)
    plt.title('Spectrum from FITS file')
    plt.xlabel('Wavelength (Angstrom)')
    plt.ylabel('Flux (ADU)')

    plt.gcf().set_size_inches(12, 6)
    plt.tight_layout()
    plt.show(block=False)

def plt_diag_lines(master_folder, flux_cal=False):
    '''Function to plot the diagnostic lines of the spectrum
        for all the master FITS files in a folder.

        Parameters
        ----------
        master_folder : str
            Path to the folder containing the master FITS files.
        flux_cal : bool, optional
            If True, plot the flux-calibrated spectra. Default is False.
    '''

    setups = {
            'DIC1+346' : (3050, 3860),
            'DIC1+580' : (4790, 6790),
            'DIC2+437' : (3760, 4980),
            'DIC2+860' : (6700, 10400),
    }
    # missing OH, C2, C3 lines?
    lines = {
            3230: 'Ti II',
            3243: 'Ti II',
            3303.5: 'Na I',
            3385: 'Ti II',
            3721: 'Fe I',
            3874: 'CN',
            3934.8: 'Ca II K',
            3969.6: 'Ca II H',
            4234: 'CH+',
            4301.5: 'CH',
            4964: 'DIB',
            5894.6: 'Na I D',
            6196: 'DIB',
            6614: 'DIB',
            7667: 'K I',
            7701: 'K I'
            }

    master_and_lines = []
    for file in os.listdir(master_folder):
        if file.startswith('MASTER_') and file.endswith('.fits'):
            path = os.path.join(master_folder, file)
            m = ob(path, orig='reduced')
            setup = m.dich + '+' + str(m.cwlen)
            if not setup in setups:
                print(f"Setup {setup} not recognized. Skipping file {file}.")
                continue

            lines_i = [line for line, name in lines.items() if setups[setup][0] <= line <= setups[setup][1]]
            master_and_lines.append((m, lines_i))
            print(lines_i)

    # make a unique list of lines to plot
    lines_f = sorted(set([line for m, lines_i in master_and_lines for line in lines_i]))

    n_lines = len(lines_f)
    n_cols = 4
    n_rows = int(np.ceil(n_lines / n_cols))

    fig, axs = plt.subplots(n_rows, n_cols, figsize=(15, 2.5*n_rows))
    axs = axs.flatten()

    i = 0
    for line in lines_f:
        for m, lines_i in master_and_lines:
            if line in lines_i:
                delta = [5 if line in [5894.6,3874] else 2][0]
                mask = (m.wave > line-delta) & (m.wave < line+delta)
                if flux_cal:
                    flux = m.flux_cal[mask]
                else:
                    flux = m.flux[mask]
                axs[i].plot(m.wave[mask], flux, alpha=0.7, lw=0.7)
                axs[i].set_title(f'{lines[line]} ({line} Å)')
                axs[i].set_xlabel('Wavelength (Å)')
                axs[i].set_ylabel('Flux')
        i += 1

    # remove empty subplots
    for j in range(i, n_rows*n_cols):
        fig.delaxes(axs[j])

    plt.tight_layout()
    plt.show(block=False)

def info(folder='/Users/adeburgo/Documents/pipelines/EDPS_data/UVES/object/'):
    '''Function to print the information of the files in a folder.
    '''
    uves_files = ['red_science_blue.fits', 'red_science_redl.fits', 'red_science_redu.fits']
    for root, dirs, files in os.walk(folder):
        for file in files:
            filename = file.split('/')[-1]
            if filename in uves_files or \
                (filename.startswith('MASTER_') and filename.endswith('.fits')):
                path = os.path.join(root, file)
                print(f'\033[94m {path.split("/")[-1]} \033[0m')
                if 'MASTER_' in filename:
                    spec = ob(path, orig='reduced')
                else:
                    spec = ob(path)
