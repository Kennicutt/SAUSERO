from reduction.reduction_osirisplus import *
from aligning.aligning_osirisplus import *
from astroMetry.astrometry_osirisplus import *
from photometry.photometry_osirisplus import *

from astropy import units as u

import argparse, time
import os, json

parser = argparse.ArgumentParser(
                     prog = 'OsirisDRP',
                     description = 'This software reduces the observations carried out by OSIRIS \
                         in BBI mode. It can reduce any filter configuration and it can be used \
                         with observation afected by fringing (Sloan_z).')

parser.add_argument('-pr','--program', help='Select the GTC program.',
                     required=True, type=str)

parser.add_argument('-bl','--block', help='Select the block of the program.',
                     required=True, type=str)

args = parser.parse_args()

def readJSON():
    return json.load(open("/home/fabricio.perez/Core/Proyectos/SAUSERO/configuration.json"))

def Results(PATH, ZP, eZP, MASK, filt, ext_info = extinction_dict):
    """ Esta función añade información relevante a la imagen de ciencia
    que se va a entregar al PI. Esta información consiste en la magnitud
    instrumental para el filtro empleado y su error, las unidades empleadas,
    el valor de la extinción utilizado y su error, dividir la imagen por
    el tiempo de exposición y aplicar la máscara.

    Args:
        PATH (string): Directorio que contiene los resultados preliminares.
        ZP (float): La magnitud instrumental estimada para el filtro.
        eZP (float): El error de la magnitud instrumental.
        MASK (bool): Máscara de píxeles malos.
        filt (string): Filtro empleado para la adquisición de las imágenes.
    """
    ic = ccdp.ImageFileCollection(PATH, keywords='*', glob_include='ast*')
    for sky in ['SKY', 'NOSKY']:
        fname = ic.files_filtered(include_path=True, filtro=filt, ssky=sky)[0]
        frame = CCDData.read(fname, unit='adu')
        print(f'Load science image: {os.path.basename(fname)}')
        hd = frame.header
        hd['ZP'] = (ZP, 'ZeroPoint estimation')
        hd['eZP'] = (eZP, 'Error ZeroPoint estimation')
        print(f'ZeroPoint information added to header')
        extinction, e_extinction = ext_info[filt]
        hd['EXT']=(extinction, 'Filter extinction')
        hd['eEXT']=(e_extinction, 'Error filter extinction')
        print(f'Extinction information added to header')
        frame.header = hd
        frame.unit = u.adu/u.second
        print(f'Change units: ADUs to ADUs/second')
        
        frame.data = (frame.data / hd['EXPTIME'])
        frame.write(PATH / f"{hd['GTCPRGID']}_{hd['GTCOBID']}_{filt}_{(hd['DATE'].split('T')[0]).replace('-','')}_{sky}_BBI.fits",
                    overwrite=True)
        print(f"Frame generated: {hd['GTCPRGID']}_{hd['GTCOBID']}_{filt}_{(hd['DATE'].split('T')[0]).replace('-','')}_{sky}_BBI.fits")
    

if __name__ == '__main__':
    PRG = args.program
    OB = args.block

    conf = readJSON()
    
    #Reduction Recipe. Esta receta se encarga de efectuar la limpieza de la imágenes
    #haciendo una sustracción del masterbias y dividiendo por el masterflat normalizado.
    #Posteriormente se salvan las imágenes ya limpias.
    print('Start the reduction...')
    o = Reduction(PRG, OB, main_path=conf['DIRECTORIES']['PATH_DATA'],
                path_mask=conf['DIRECTORIES']['PATH_BPM'])
    o.get_imagetypes()
    o.load_BPM()
    o.sort_down_drawer()
    o.do_masterbias()
    o.do_masterflat()
    o.get_std(no_CRs=conf['REDUCTION']['no_CRs'], contrast_arg = conf['REDUCTION']['contrast'],
            cr_threshold_arg = conf['REDUCTION']['cr_threshold'],
            neighbor_threshold_arg = conf['REDUCTION']['neighbor_threshold'])
    o.get_target(no_CRs=conf['REDUCTION']['no_CRs'], contrast_arg = conf['REDUCTION']['contrast'],
            cr_threshold_arg = conf['REDUCTION']['cr_threshold'],
            neighbor_threshold_arg = conf['REDUCTION']['neighbor_threshold'])
    o.remove_fringing()
    o.sustract_sky()
    o.save_target(std=conf['REDUCTION']['save_std'])
    o.save_target(sky=conf['REDUCTION']['save_sky'])
    o.save_target(fringing=conf['REDUCTION']['save_fringing'])
    o.save_target()
    print('Finished the reduction successfully')
    print(2*"\n")
    
    #Aligned Recipe. Las imágenes de ciencia limpias son alineadas en función del
    #filtro empleado en cada caso. Luego se salva como imagen alienada.
    print("Start the alignment...")
    al = OsirisAlign(PRG, OB, conf)
    for filt in list(set(al.ic.summary['filtro'])):
        for sky in ['SKY', 'NOSKY']:
            print(f'Aligment for {filt} & {sky}')
            align = al.aligning(filt, sky=sky)
            lst = al.load_frames(filt, sky=sky)
            fr = CCDData.read(lst[0], unit='adu')
            header = fr.header
            header['exptime'] = al.total_exptime * (al.num + 1.)
            wcs = fr.wcs
            save_fits(align, header, wcs, al.PATH_REDUCED / f'aligned_result_{filt}_{sky}.fits')
            print(f'Alineado para {filt} & {sky} realizado exitosamente')

    print('El alineado para cada filtro finalizado')
    print(2*"\n")

    #Astrometry Recipe. La imagen alineada para cada filtro se aplica un proceso de
    #astrometrización de la imagen, para determinar con precisión la posición real
    #de los cuerpos celestes presentes en la escena.
    print("Start the astrometrization...")
    del filt
    ic_ast = ccdp.ImageFileCollection(al.PATH_REDUCED, keywords='*', glob_include='ali*', glob_exclude='*NOSKY*')
    lst_filt = list(ic_ast.summary['filtro'])
    for filt in lst_filt:
        print(f'Astrometrization for {filt}')
        best_wcs, new_frame = solving_astrometry(PRG, OB, filt, conf, sky='SKY', calib_std=False)
        for sky in ['SKY', 'NOSKY']:
            if sky == 'SKY':
                new_frame.write(al.PATH_REDUCED / f'ast_result_{filt}_{sky}.fits', overwrite=True)
            else:
                nosky = CCDData.read(al.PATH_REDUCED / f'aligned_result_{filt}_{sky}.fits', unit='adu')
                nosky.wcs = WCS(best_wcs)
                nosky.write(al.PATH_REDUCED / f'ast_result_{filt}_{sky}.fits', overwrite=True)
            print(f'Astrometrization done successfully for {filt}')
            time.sleep(10)
    
    print('Start astrometrization for STD star')
    time.sleep(30)
    ic_std = ccdp.ImageFileCollection(al.PATH_REDUCED, keywords='*', glob_include='*std*', glob_exclude='*NOSKY*')
    best_wcs_std, new_frame_std = solving_astrometry(PRG, OB, filt, conf, sky='SKY', calib_std=True)
    for path_to_std in ic_std.files_filtered(include_path=True):
        print(f'Astrometrization for: {path_to_std}')
        std_img = CCDData.read(path_to_std, unit='adu')
        std_img.wcs = WCS(best_wcs_std)
        std_img.write(path_to_std, overwrite=True)


    print('Astrometrization finished')
    print(2*"\n")


    #Photometry Recipe. A partir de la estrella de calibración, se estima la magnitud
    #instrumental. De este modo, podemos estimar posteriormente la magnitud aparente
    #de los cuerpos celestes presentes en la imagen de ciencia. Este proceso se hace
    # por cada filtro empleado. Además, la función Results() permite incluir esta
    #información en el header de la imagen de ciencia limpia, alineada y astrometrizada.
    print("Start the estimation of ZeroPoint...")
    del filt
    ic_pho = ccdp.ImageFileCollection(al.PATH_REDUCED, keywords='*', glob_include='red*')
    lst_filt = list(set(ic_pho.summary['filtro']))
    for filt in lst_filt:
        print(f'Filter selected is {filt}')
        ZP, eZP = photometry(PRG, OB, ic_pho.files_filtered(imgtype='STD',filter2=filt)[0], conf)
        print(f'Estimated ZP: {ZP} +- {eZP} for {filt}')
        Results(al.PATH_REDUCED, ZP, eZP, o.MASK, filt)
    
    print('End of the reduction. The results are available in reduced directory.')
