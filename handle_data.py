import os
import shutil
from astropy.io import fits

def split_red_blue_data(root_dir):
    ''' Split FITS files in the given directory into two subdirectories based on
        the CCD type (RED or BLUE) indicated in the "HIERARCH ESO PRO CATG" header keyword.

        Parameters
        ----------
        root_dir : str
            The root directory containing FITS files to be organized.
    '''

    # Walk through the directory tree
    for root, dirs, files in os.walk(root_dir):
        for file in files:
            # Check if it is a FITS file (case insensitive)
            if not file.lower().endswith('.fits') or file.startswith('._'):
                continue
            file_path = os.path.join(root, file)

            # Open the FITS file to read the header
            with fits.open(file_path) as hdul:
                header = hdul[0].header
                if len(hdul) == 3:
                    header.update(hdul[1].header)

                # Get CCD type from "HIERARCH ESO PRO CATG"
                if 'HIERARCH ESO PRO CATG' in header:
                    ccd = header['HIERARCH ESO PRO CATG'].split('_')[-1]
                else:
                    ccd = header['HIERARCH ESO INS PATH'].strip()                
                if 'RED' in ccd:
                    ccd_type = 'RED'
                elif 'BLUE' in ccd or 'BLU' in ccd:
                    ccd_type = 'BLUE'
                else:
                    print(f"Warning: Could not determine CCD type for {file}. Skipping.")
                    continue

            # Create the specific CCD folder inside the root directory
            target_folder = os.path.join(root_dir, ccd_type)
            os.makedirs(target_folder, exist_ok=True)

            # Copy the file
            destination_path = os.path.join(target_folder, file)
            shutil.move(file_path, destination_path)

            print(f"Moved: {file} -> {ccd_type}/")


def organize_EDPSdata(destination_folder):
    '''
    Organize FITS files from the EDPS pipeline by matching them with the corresponding
    files in the final_output folder. The function looks for files in the final_output
    folder, extracts the ARCFILE keyword from their headers, and then searches for
    matching files in the EDPS_data folder based on the same ARCFILE value. If a match
    is found, both the final_output file and the corresponding EDPS file are copied to
    the specified destination folder. In the case of the RED CCD the function takes into
    account that there is a _REDL_ and _REDU_ version of the same file by looking at the
    "ESO PRO CATG" keyword and matching it with the corresponding one in the EDPS files.

        Parameters
        ----------
        destination_folder : str
            The directory where the matched files will be copied.
    '''

    edps_path = '/Users/adeburgo/Documents/pipelines/'
    if not destination_folder.endswith('/'):
        destination_folder += '/'

    for root, dirs, files in os.walk(edps_path + 'final_output/'):
        for file in files:
            if file.startswith('._') or not file.endswith('.fits'):
                continue
            final_output_path = os.path.join(root, file)
            with fits.open(final_output_path) as hdul:
                header = hdul[0].header
                arcfile = header['ARCFILE']
                catg = header.get('HIERARCH ESO PRO CATG', '').split('_')[-1]
                print(f"Found {file} with ARCFILE: {arcfile} and CATG: {catg}")

            for root2, dirs2, files2 in os.walk(edps_path + 'EDPS_data/UVES/object/'):
                for file2 in files2:
                    if file2.startswith('._') or not file2.endswith('.fits'):
                        continue
                    if file2.startswith('fluxcal_science'):
                        edps_path2 = os.path.join(root2, file2)
                        with fits.open(edps_path2) as hdul2:
                            header2 = hdul2[0].header
                            arcfile2 = header2['ARCFILE']
                            catg2 = header2.get('HIERARCH ESO PRO CATG', '').split('_')[-1]
                        if arcfile == arcfile2 and catg == catg2:
                            shutil.copy2(final_output_path, destination_folder)
                            shutil.copy2(edps_path2, destination_folder + \
                                        final_output_path.split('/')[-1].replace('SPECTRUM', 'FLUXCAL'))
                            print(f"Copied {final_output_path.split('/')[-1]} and {edps_path2.split('/')[-1]} to {destination_folder}")
                            break
                else:
                    continue
                break

def organize_fits_in_EDPSfolder_by_object(root_dir):
    ''' Organize FITS files in the given directory into
        subdirectories based on the OBJECT header keyword.

        Parameters
        ----------
        root_dir : str
            The root directory containing FITS files to be organized.
        '''

    # Name of the new directory where sorted folders will live
    output_base_name = "final_output"
    output_dir = os.path.join(root_dir, output_base_name)

    # Create the main output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    print(f"Scanning directory: {root_dir+'EDPS_data/'}")
    print(f"Output directory: {output_dir}")
    print("-" * 30)

    files_processed = 0
    errors = 0

    # Walk through the directory tree
    for root, dirs, files in os.walk(root_dir+'/EDPS_data/'):

        # Skip the output directory itself so we don't recursively process copies
        if output_base_name in root:
            continue

        for file in files:
            # Check if it is a FITS file (case insensitive)
            if not file.lower().endswith('.fits') or file.startswith('._'):
                continue
            if file.startswith('fluxcal_science') or file.startswith('merged_science'):
                file_path = os.path.join(root, file)

                try:
                    # Open the FITS file to read the header
                    # 'memmap=True' helps with large files by not loading the data segment
                    with fits.open(file_path, memmap=True) as hdul:
                        # Usually the OBJECT keyword is in the Primary Header (index 0)
                        header = hdul[0].header

                        # Get OBJECT value, default to "UNKNOWN" if missing
                        obj_name = header['OBJECT'].strip().replace(' ', '')

                        # Get PIPEFILE and ARCFILE
                        pipefile = header['PIPEFILE'].strip().replace(' ', '')
                        arcfile = header['ARCFILE'].strip().replace(' ', '')

                    # SANITIZATION:
                    # 1. Remove spaces
                    folder_name = obj_name.replace(" ", "")
                    # 2. (Optional) Remove slashes to prevent path errors
                    folder_name = folder_name.replace("/", "_").replace("\\", "_")

                    # Create the specific Object folder inside the output directory
                    target_folder = os.path.join(output_dir, folder_name)
                    os.makedirs(target_folder, exist_ok=True)

                    # Rename the file depending on its type
                    file = pipefile.replace('.fits', '_').upper() + arcfile

                    # Copy the file
                    # shutil.copy2 preserves file metadata (timestamps, etc.)
                    destination_path = os.path.join(target_folder, file)
                    shutil.copy2(file_path, destination_path)

                    print(f"Copied: {file} -> {folder_name}/")
                    files_processed += 1

                except Exception as e:
                    print(f"Error processing {file}: {e}")
                    errors += 1

    print("-" * 40)
    print(f"Done. Processed {files_processed} files with {errors} errors.")

def organize_fits_by_object_in_folders(root_dir):
    ''' Organize FITS files in the given directory into
        subdirectories based on the OBJECT header keyword.

        Parameters
        ----------
        root_dir : str
            The root directory containing FITS files to be organized.
    '''

    # Walk through the directory tree
    for root, dirs, files in os.walk(root_dir):
        for file in files:
            # Check if it is a FITS file (case insensitive)
            if not file.lower().endswith('.fits') or file.startswith('._'):
                continue
            file_path = os.path.join(root, file)

            try:
                # Open the FITS file to read the header
                # 'memmap=True' helps with large files by not loading the data segment
                with fits.open(file_path, memmap=True) as hdul:
                    # Usually the OBJECT keyword is in the Primary Header (index 0)
                    header = hdul[0].header

                    # Get OBJECT value, default to "UNKNOWN" if missing
                    obj_name = header['OBJECT'].strip().replace(' ', '')

                # Create the specific Object folder inside the root directory
                target_folder = os.path.join(root_dir, obj_name)
                os.makedirs(target_folder, exist_ok=True)

                # Copy the file
                # shutil.copy2 preserves file metadata (timestamps, etc.)
                destination_path = os.path.join(target_folder, file)
                shutil.copy2(file_path, destination_path)

                print(f"Copied: {file} -> {obj_name}/")

            except Exception as e:
                print(f"Error processing {file}: {e}")


def fake_master_response(file):
    #! FOR SOME REASON I HAD TO MAKE A COPY OF THE ORIGINAL FILE AND EDIT THE COPY
    shutil.copyfile(file, file.replace('.fits', '_fix.fits'))

    # Manually fix headers of the STD response file:
    hdu = fits.open(file.replace('.fits', '_fix.fits'))

    hdu0 = hdu[0]
    hdu1 = hdu[1]
    hdu2 = hdu[2]

    header0 = hdu0.header

    color = header0['HIERARCH ESO PRO CATG'].split('_')[-1]

    header0['HIERARCH ESO PRO CATG'] = 'MASTER_RESPONSE_'+color
    header1 = hdu1.header
    header1['TTYPE1'] = 'LAMBDA'
    header1['TTYPE2'] = 'FLUX_CONV'

    data1 = hdu1.data
    data2 = hdu2.data

    data1['LAMBDA'] = data2['wavelength']
    data1['FLUX_CONV'] = data2['response_raw']

    hdu.writeto(file.replace('.fits', '_fix.fits'), overwrite=True)
    hdu.close()

def fake_date_flat(file, new_date):
    # check that the new_date is in the format YYYY-MM-DD
    if len(new_date) != 10 or new_date[4] != '-' or new_date[7] != '-':
        print("Error: new_date must be in the format YYYY-MM-DD")
        return

    old_date = file.split('UVES.')[1].split('T')[0]
    shutil.copyfile(file, file.replace(old_date, new_date))

    # Manually fix headers of the flat file:
    hdu = fits.open(file.replace(old_date, new_date), mode='update')

    hdu0 = hdu[0]

    header0 = hdu0.header

    old_date = header0['DATE'].split('T')[0]
    header0['DATE'] = header0['DATE'].replace(old_date, new_date)
    old_date = header0['DATE-OBS'].split('T')[0]
    header0['DATE-OBS'] = header0['DATE-OBS'].replace(old_date, new_date)
    old_date = header0['HIERARCH ESO OBS START'].split('T')[0]
    header0['HIERARCH ESO OBS START'] = header0['HIERARCH ESO OBS START'].replace(old_date, new_date)
    old_date = header0['HIERARCH ESO TPL START'].split('T')[0]
    header0['HIERARCH ESO TPL START'] = header0['HIERARCH ESO TPL START'].replace(old_date, new_date)
    old_date = header0['ARCFILE'].split('UVES.')[1].split('T')[0]
    header0['ARCFILE'] = header0['ARCFILE'].replace(old_date, new_date)

    header0['LST'] = 30069.702
    header0['MJD-OBS'] = 60901.645622128
    header0['UTC'] = 55773.000

    hdu.flush()
    hdu.close()

