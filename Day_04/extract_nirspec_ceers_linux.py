

import matplotlib as mpl
mpl.use('TkAgg')
#mpl.use('MacOSX')
#mpl.use('Qt4Agg')
import numpy as np
from matplotlib import pyplot as pl
from astropy.io.ascii import read
import astropy.io.fits as pf
from scipy.optimize import curve_fit
from matplotlib.backend_bases import MouseButton
import os
from glob import glob
from time import gmtime

import warnings
warnings.filterwarnings("ignore")

import sys

pl.rcParams['keymap.pan'].remove('p')
pl.rcParams['keymap.save'].remove('s')
pl.rcParams['keymap.fullscreen'].remove('f')
pl.rcParams['keymap.home'].remove('r')
pl.rcParams['keymap.zoom'].remove('o')

cont = False

pointing = 'P8'
#pointing = 'P7'
#pointing = 'P9'
#pointing = 'P10'
#pointing = 'P4'
#pointing = 'P5'

pointing=sys.argv[1]

# path to 2D spectrum files
# need to change path2d to "<path-to-directory-containing-pointing-directory>/<pointing>/2D/" where pointing is, e.g., "P8"
path2d = '/media/home/team_workspaces/JWST-Heidelberg-Summer-School/Lecturers_Area/08_Thursday_Session_3/Data/%s/2D/'%pointing

# path to output directory for extracted 1D spectra
#path1d = path2d[:-3] + '1D/'
path1d = '1D/'
if not os.path.exists(path1d): os.system('mkdir %s'%path1d)

def inspect2d(id):
    global ax, origlims
    global lastpressed
    global keypresses2d
    global vmarr, im2darr, obswlarr
    global ax2drat, fig
    global z, restwlax
    global zax, zxlow, zxhigh, zylow, zyhigh, zoomcomplete
    global extracted, saved
    global plotyexp, yexparr

    print('Working on ID = %s'%id)

    czcat = read('input_redshift.dat')
    czcat = czcat[(czcat['ID_CEERS'] == id) & (czcat['pointing'] == pointing)]
    if len(czcat) == 1: print('Redshift estimates: z_CEERS = %.3f\tz_3DHST = %.3f'%(czcat['z_CEERS'], czcat['z_3DHST']))
    elif len(czcat) == 0: print('No match in the photo-z catalog!')
    else:
        print('Multiple matches in the photo-z catalog.  Something    iis weird here...')
        [print('Redshift estimates: z_CEERS = %.3f\tz_3DHST = %.3f'%(czcat['z_CEERS'].data[zz], czcat['z_3DHST'].data[zz])) for zz in range(len(czcat))]

    # read in all 2D files for this ID and parse for the program name, target ID, mask number, grating config, and detector
    list2d = glob(path2d + '%s/*%s*_2d.fits'%(id,id))
    progarr = []
    idarr = []
    masknamearr = []
    detarr = []
    gconfigarr = []
    for ll in list2d:
        splitname = ll.strip().split('/')[-1].split('_')
        splitname2 = splitname[2].split('-')
        progarr += [splitname[0]]
        idarr += [int(splitname2[0])]
        masknamearr += [splitname[1]]
        gconfigarr += [splitname2[1]+'-'+splitname2[2]]
        detarr += [splitname2[3]]
    progarr = np.array(progarr)
    idarr = np.array(idarr)
    masknamearr = np.array(masknamearr)
    gconfigarr = np.array(gconfigarr)
    detarr = np.array(detarr)

    # number of grating configurations
    gconfigs = np.unique(gconfigarr)
    ngrat = len(gconfigs)

    """
    fig = pl.figure(1, figsize=(16, 9))
    fig.suptitle('ID = %s'%id)

    top = 0.9
    bottom = 0.1
    left = 0.05
    right = 0.95

    hspace = 0.1
    height = (top-bottom-float(ngrat-1)*hspace)/float(ngrat)
    ax2drat = height/(right - left)

    #ax = [fig.add_axes([left, bottom+gg*(height+hspace), right-left, height]) for gg in range(ngrat)]
    #ax = ax[::-1]
    #ax = np.array(ax)

    ax = []
    for gg in range(ngrat):
     aa = fig.add_axes([left, bottom+gg*(height+hspace), right-left, height])
     ax += [aa]
    ax = ax[::-1]
    ax = np.array(ax)
    """

    gratarr = []
    wlarr = []
    specarr = []
    errarr = []
    hdrarr = []
    vmarr = []
    origvmarr = []
    #im2darr = []
    obswlarr = []
    crval1arr = []
    cdelt1arr = []
    shdrarr = []
    # loop through each grating config
    #for i in [2]:
    print('Number of grating configurations = %s'%ngrat)
    print(gconfigs)
    for i in range(ngrat):
        progs = progarr[gconfigarr == gconfigs[i]]
        ids = idarr[gconfigarr == gconfigs[i]]
        masknames = masknamearr[gconfigarr == gconfigs[i]]
        grats = gconfigarr[gconfigarr == gconfigs[i]]
        dets = detarr[gconfigarr == gconfigs[i]]
        # pull data from both detectors in this grating, if available
        wlsub = []
        specsub = []
        errsub = []
        hdrsub = []
        for dd in range(len(dets)):
            inname = path2d + '%s/%s_%s_%s-%s-%s_2d.fits'%(ids[dd], progs[dd], masknames[dd], ids[dd], grats[dd], dets[dd])
            hdu = pf.open(inname)
            specsub += [np.array(hdu['SCI'].data)]
            shdr = hdu['SCI'].header
            hdrsub += [shdr]
            errsub += [np.array(hdu['ERR'].data)]
            wlsub += [shdr['crval1'] + np.arange(shdr['naxis1'])*shdr['cdelt1']]
        #if len(specarr) > 1:
        # if np.shape(specarr[0]) != np.shape(specarr[1]): print('ERROR: The shape of the 2D spectrum is not the same for data in each detector!')
        # if wlarr[0] != wlarr[1]: print('ERROR: The wavelength solutions do not match for data in each detector!')

        wl = wlsub[0]
        hdr = hdrsub[0]
        minny = np.amin([len(ss) for ss in specsub])
        spec = np.sum([ss[:minny,:] for ss in specsub], axis=0)
        for ee in errsub:
            ee[np.logical_not(np.isfinite(ee))] = 0.
        err = np.sum([ee[:minny,:] for ee in errsub], axis=0)
        #spec = np.sum(specsub, axis=0)
        #err = np.sum(errsub, axis=0)
        err[err == 0.] = np.inf
        crval1arr += [hdr['crval1']]
        cdelt1arr += [hdr['cdelt1']]

        gratarr += [grats[0]]
        wlarr += [wl]
        specarr += [spec]
        errarr += [err]
        hdrarr += [hdr]
        program = progarr[0]
        maskname = masknamearr[0]


        """
        # 2D spectrum
        vm = np.sort(spec.flatten()[spec.flatten() != 0.])[int(0.97*len(spec.flatten()[spec.flatten() != 0.]))]
        vmarr += [vm]
        origvmarr += [vm]
        im2darr += [ax[i].imshow(spec, vmin=-vm, vmax=vm, cmap='viridis', origin='lower')]
        ax[i].set_title('%s'%grats[0])

        def pix2wl(i, z, norm):
         def func(pix): return (crval1arr[i]+pix*cdelt1arr[i])/(1.+z)*norm
         return func
        def wl2pix(i, z, norm):
         def func(wl): return (wl/norm*(1.+z)-crval1arr[i])/cdelt1arr[i]
         return func

        #sax2d = ax[i].secondary_xaxis('top', functions=(lambda pix: crval1arr[i]+pix*cdelt1arr[i], lambda wl: (wl-crval1arr[i])/cdelt1arr[i]))
        sax2d = ax[i].secondary_xaxis('top', functions=(pix2wl(i, 0., 1.), wl2pix(i, 0., 1.)))
        sax2d.set_xlabel(r'$\lambda_{\mathrm{obs}}$ ($\mu$m)', fontsize=12, labelpad=3)
        #ax2d.tick_params(bottom=False, labelbottom=False)
        obswlarr += [sax2d]
        """

    #vmarr = np.array(vmarr)
    #im2darr = np.array(im2darr)
    #origvmarr = np.array(origvmarr)

    # list 1D files that have been extracted
    list1d = glob(path1d + '*-%s-*-*_1d.fits'%(id))
    exgrats = []
    for ll in list1d:
        for gg in gratarr:
            if gg in ll: exgrats += [gg]
    extracted = [True if gg in exgrats else False for gg in gratarr] # if gg in exgrats else False]

    #extracted = [False for gg in range(ngrat)]
    saved = False
    boxcar = False

    #origlims = np.array([aa.axis() for aa in ax])

    # interactivity functionality
    lastpressed = None

    zax = None
    zxlow = None
    zxhigh = None
    zylow = None
    zyhigh = None
    zoomcomplete = True

    excomplete = True
    exfitlow = None
    exfithigh = None

    fitcomplete = True
    lfit = None
    lerr = None
    fxlow = None
    fxhigh = None

    #restwlax = [False for aa in range(len(ax))]
    wlax2dflag = 0  # 0 if pixels, 1 if observed wl, 2 if rest wl

    plotyexp = True

    spec1d = []
    var1d = []
    spec1dbox = []
    var1dbox = []

    z = 0.

    def onkeypress2d(event):
        global ax, origlims
        global lastpressed
        global keypresses2d
        global vmarr, im2darr, obswlarr
        global ax2drat
        global z, restwlax
        global zax, zxlow, zxhigh, zylow, zyhigh, zoomcomplete
        global saved
        global yexparr, plotyexp

        # reset the axes
        if event.key == 'r':
            pl.close(fig)
            reset2d()
            #event.inaxes.axis(origlims[ax == event.inaxes][0])
            #im2darr[ax == event.inaxes][0].set_clim(vmin=-origvmarr[ax == event.inaxes], vmax=origvmarr[ax == event.inaxes])
            #vmarr = np.copy(origvmarr)
            #pl.draw()

        # print the extraction windows, redshift catalog
        #if event.key == 'p':

        # toggle expected y position of main target on and off
        if event.key == 'y':
            if plotyexp:
                [yy.remove() for yy in yexparr]
                pl.draw(); plotyexp = False
            else:
                [ax[yy].add_line(yexparr[yy]) for yy in range(len(yexparr))]
                pl.draw(); plotyexp = True


        # extract 1D spectrum for the chosen grating
        if event.key == 'e':
            eind = np.where( ax == event.inaxes )[0][0]
            outroot = '%s-%s-%s-%s'%(program, maskname, id, gratarr[eind])

            # disconnect interactive functions
            fig.canvas.mpl_disconnect(keypresses2d)

            pl.close(fig)

            saved = extract1d(id, wlarr[eind], specarr[eind], errarr[eind], hdrarr[eind], path2d, path1d, outroot, gratarr[eind])
            extracted[eind] = saved
            print('1d extraction finished.')
            reset2d()

        # use a pixel location of a line to estimate redshift
        if event.key == 'f':
            xpix = event.xdata
            find = np.where( ax == event.inaxes )[0][0]
            lambdaobs = pix2wl(find, 0., 1.e4)(xpix)
            print('Cursor is at x=%.1f pixels, or lambda_obs=%.1f Angstroms'%(xpix,lambdaobs))

            # oii, neiii, hd, hg, hb, oiii4959, oiii5007, ha, nii6585, sii6718, sii6731, siii9069, siii9531
            linenamearr = ['Lyalpha', 'CIII]1908', '[OII]doublet', '[NeIII]3869', 'Hdelta', 'Hgamma', 'Hbeta', '[OIII]4959', '[OIII]5007', '[NII]6548', 'Halpha', '[NII]6584', '[SII]6716', '[SII]6731', '[SIII]9069', '[SIII]9531', 'HeI10830', 'Pabeta', 'Paalpha']
            linecent0arr = [1215.67, round((1906.68+1908.73)/2.,3), round((3727.092+3729.875)/2.,3), 3869.856, 4102.89, 4341.68, 4862.68, 4960.295, 5008.240, 6549.86, 6564.61, 6585.27, 6718.29, 6732.67, 9071.09, 9533.21, 10833.305, 12821.578, 18756.420]
            linenumarr = ['%s'%ll for ll in range(len(linenamearr))]
            linestrarr = ['%.0f:\t%-13s\tcent0 = %s AA\tz = %.5f\n'%(ll, linenamearr[ll], linecent0arr[ll], lambdaobs/linecent0arr[ll]-1.) for ll in range(len(linenamearr))]
            linestr = ''.join(linestrarr)

            # UV lines
            #linecents = [1548.19+1., 1550.77+1., 1640.42, 1660.81, 1666.15, 1749.67, 1752.16, 1882.47, 1892.03, 1906.68, 1908.73, 2321.66, 2325.40, 2326.93, 2328.12, 2471.03]
            #linenames = ['civ1548', 'civ1550', 'heii1640', 'oiii1661', 'oiii1666', 'niii1750', 'niii1752', 'siiii1882', 'siiii1892', 'ciii1907', 'ciii1909', 'oiii2322', 'cii2325', 'cii2327', 'cii2328', 'oii2471']
            # optical lines
            #linecents += [(3727.092+3729.875)/2., 3869.856, (3890.166+3889.750)/2., (3971.198+3968.593)/2., 4077.500, 4102.89, 4341.68, 4364.436, 4687.021, 4862.68, 4960.295, 5008.240, 5756.24, 5877.25, 6302.046, 6313.806, 6365.535, 6549.86, 6564.61, 6585.27, 6679.995, 6718.29, 6732.67, 7067.138, 7137.767, 7321.937, 7332.210, 9071.09, 9533.21]
            #linenames += ['oii3727', 'neiii3869', 'hihei3888', 'hineiii3970', 'sii4076', 'hd', 'hg', 'oiii4363', 'heii4686', 'hb', 'oiii4959', 'oiii5007', 'nii5755', 'hei5876', 'oi6300', 'siii6312', 'oi6363', 'nii6548', 'ha', 'nii6584', 'hei6680', 'sii6716', 'sii6731', 'hei7067', 'ariii7135', 'oii7320', 'oii7330', 'siii9069', 'siii9531']

            # disconnect interactive functions
            fig.canvas.mpl_disconnect(keypresses2d)

            num = input('\nIdentify the line:\n\n' + linestr + 'c: Custom line rest wavelength\nEnter selection: ')

            if num == 'c':  # enter custom line name/rest wavelength
                linecent0 = float(input('Enter rest-frame wavelength in Angstroms: '))
                z = lambdaobs/linecent0 - 1.
                print('Rest wavelength = %.3f'%linecent0)
                print('Updating redshift guess to z=%.4f'%z)
            elif num not in linenumarr:
                print('Invalid selection.  Exiting without updating redshift guess.')
            else:
                num = int(num)
                linecent0 = linecent0arr[num]; linename = linenamearr[num]
                z = lambdaobs/linecent0 - 1.
                print('Rest wavelength = %.3f'%linecent0)
                print('Updating redshift guess to z=%.4f'%z)

            # reconnect interactive functions
            keypresses2d = fig.canvas.mpl_connect('key_press_event', onkeypress2d)

        # modify redshift catalog
        #if event.key == 'm':

        if event.key == 'h' or event.key == '?':
            helpstr = '\n In the 2D science spectrum panel, click and drag to add pixels in that x-range to the spatial profile.\n\
       Key bindings:\n\
       r -- "reset": Reset the 2D spectrum images.\n\
       y -- "yposition": Toggle marker for the expected y-position of the primary target (orange dotted line).\n\
       e -- "extract": Place the cursor over a panel and press "e" to begin 1D extraction for that grating configuration.\n\
       f -- "fit": Identify a line or rest-frame wavelength at the pixel position of the cursor to estimate the redshift.  Once a redshift has been estimated, the rest-frame wavelength can be displayed for all panels with "w".\n\
       w -- "wavelength": Press "w" to toggle between observed wavelength and rest-frame wavelength for the selected panel.  This feature requires a redshift estimate using "f".\n\
       z -- "zoom": Press "z" twice in the same panel to zoom in to a range with corners defined by the cursor positions when "z" was pressed.\n\
       x -- "axis reset": Resets the plotted range to the original range for the axis in which the cursor is.\n\
       c -- "center": Move the point under the cursor to the center of the panel.\n\
       i -- "in": Zoom the selected axis in, keeping the same center coordinate.\n\
       o -- "out": Zoom the selected axis out, keeping the same center coordinate.\n\
       q -- "quit": Exit and return to the object selection menu.\n\
       arrow keys: Use the arrow keys to shift the plot.  The 2D spectrum panel can only be shifted left and right, while up and down changes the stretch.\n'
            print(helpstr)

        # toggle rest wavelength in ax1d based on z_mean from current redshift catalog
        if event.key == 'w':
            wind = np.where( ax == event.inaxes )[0][0]
            if not restwlax[wind]:  # convert to rest wavelength
                zm = z
                obswlarr[wind].remove()
                sax2d = ax[wind].secondary_xaxis('top', functions=(pix2wl(wind, z, 1.e4), wl2pix(wind, z, 1.e4)))
                sax2d.set_xlabel(r'$\lambda_{\mathrm{rest}}$ ($\AA$) $z=%.4f$'%z, fontsize=12, labelpad=3)
                obswlarr[wind] = sax2d
                pl.draw()
                restwlax[wind] = True
            elif restwlax[wind]:  # convert to observed wavelength
                zm = z
                obswlarr[wind].remove()
                sax2d = ax[wind].secondary_xaxis('top', functions=(pix2wl(wind, 0., 1.), wl2pix(wind, 0., 1.)))
                sax2d.set_xlabel(r'$\lambda_{\mathrm{obs}}$ ($\mu$m)', fontsize=12, labelpad=3)
                obswlarr[wind] = sax2d
                pl.draw()
                restwlax[wind] = False

        # zoom in an axis
        if event.key == 'z':
            if zax != event.inaxes: zoomcomplete = True
            if lastpressed != 'z' or (lastpressed == 'z' and zoomcomplete == True):  # set lower bounds
                zax = event.inaxes
                zxlow, zylow = event.xdata, event.ydata
                zoomcomplete = False
            elif lastpressed == 'z' and zoomcomplete == False:  # set higher bounds
                if event.inaxes == zax:
                    zxhigh, zyhigh = event.xdata, event.ydata
                    zxcent = (zxhigh+zxlow)/2.
                    figsize = fig.get_size_inches()
                    zxrange = float(shdr['naxis2'])/ax2drat*(figsize[0]/figsize[1])
                    event.inaxes.set_xlim([zxcent-zxrange/2., zxcent+zxrange/2.])
                    pl.draw()
                zoomcomplete = True

        # reset axis zoom
        if event.key == 'x':
            event.inaxes.axis(origlims[ax == event.inaxes][0]); pl.draw()

        # pan center of axis
        if event.key == 'c':
            newcx, newcy = event.xdata, event.ydata
            oldxlow, oldxhigh = event.inaxes.get_xlim()
            oldylow, oldyhigh = event.inaxes.get_ylim()
            oldxrange = (oldxhigh - oldxlow)
            oldyrange = (oldyhigh - oldylow)
            event.inaxes.set_xlim([newcx - oldxrange/2., newcx + oldxrange/2.])
            pl.draw()

        # pan axis with arrow keys
        if event.key == 'left':
            oldxlow, oldxhigh = event.inaxes.get_xlim()
            oldylow, oldyhigh = event.inaxes.get_ylim()
            oldxrange = (oldxhigh - oldxlow)
            oldyrange = (oldyhigh - oldylow)
            event.inaxes.axis([oldxlow-0.1*oldxrange, oldxhigh-0.1*oldxrange, oldylow, oldyhigh])
            pl.draw()
        if event.key == 'right':
            oldxlow, oldxhigh = event.inaxes.get_xlim()
            oldylow, oldyhigh = event.inaxes.get_ylim()
            oldxrange = (oldxhigh - oldxlow)
            oldyrange = (oldyhigh - oldylow)
            event.inaxes.axis([oldxlow+0.1*oldxrange, oldxhigh+0.1*oldxrange, oldylow, oldyhigh])
            pl.draw()
        # increase or decrease stretch scale of 2D spectrum
        if event.key == 'up':
            vmarr[ax == event.inaxes] *= 2.
            im2darr[ax == event.inaxes][0].set_clim(vmin=-vmarr[ax == event.inaxes], vmax=vmarr[ax == event.inaxes])
            pl.draw()
        if event.key == 'down':
            vmarr[ax == event.inaxes] /= 2.
            im2darr[ax == event.inaxes][0].set_clim(vmin=-vmarr[ax == event.inaxes], vmax=vmarr[ax == event.inaxes])
            pl.draw()

        # zoom in and out with fixed center
        if event.key == 'i':
            oldxlow, oldxhigh = event.inaxes.get_xlim()
            oldylow, oldyhigh = event.inaxes.get_ylim()
            oldxrange = (oldxhigh - oldxlow)
            oldyrange = (oldyhigh - oldylow)
            oldxcent = (oldxlow + oldxhigh)/2.
            oldycent = (oldylow + oldyhigh)/2.
            figsize = fig.get_size_inches()
            minxrange = float(shdr['naxis2'])/ax2drat*(figsize[0]/figsize[1])
            newxrange = max([0.5*oldxrange,minxrange])
            event.inaxes.axis([oldxcent-newxrange/2.,oldxcent+newxrange/2.,oldylow,oldyhigh])
            pl.draw()

        # zoom in and out with fixed center
        if event.key == 'o':
            oldxlow, oldxhigh = event.inaxes.get_xlim()
            oldylow, oldyhigh = event.inaxes.get_ylim()
            oldxrange = (oldxhigh - oldxlow)
            oldyrange = (oldyhigh - oldylow)
            oldxcent = (oldxlow + oldxhigh)/2.
            oldycent = (oldylow + oldyhigh)/2.
            figsize = fig.get_size_inches()
            minxrange = float(shdr['naxis2'])/ax2drat*(figsize[0]/figsize[1])
            newxrange = max([2.0*oldxrange,minxrange])
            event.inaxes.axis([oldxcent-newxrange/2.,oldxcent+newxrange/2.,oldylow,oldyhigh])
            pl.draw()

        # save extracted spectrum and redshift catalog
        #if event.key == 's':
        # print ('Saving output')
        # writezcat()
        # writespec1d(wl, spec1d, err1d)

        # auto extraction for remaining grating configurations
        if event.key == 'b':
            # list 1D files that have been extracted
            list1d = glob(path1d + '*-%s-*-*_1d.fits'%(id))
            manexgrats = []
            for ll in list1d:
                for gg in gratarr:
                    if gg in ll: manexgrats += [gg]
            autoexgrats = [gg for gg in gratarr if gg not in manexgrats]
            print('\nBeginning automatic extraction of remaining grating configurations.')
            print('Manually extracted 1D files: %s'%list1d)
            print('Available grating configurations from 2D files: %s'%gratarr)
            print('Automatically extracting these grating configurations: %s'%autoexgrats)
            excat = read(path1d + '00_extract_info.txt')
            excat = excat[excat['id'].data == id]
            print(id)
            print(excat)
            manypos = np.mean(excat['ypos'].data)
            manyfwhm = np.mean(excat['yfwhm'].data)
            manymin = np.amin(excat['ymin'].data)
            manymax = np.amax(excat['ymax'].data)

            # disconnect interactive functions
            fig.canvas.mpl_disconnect(keypresses2d)
            pl.close(fig)

            for ag in autoexgrats:
                print(ag)
                print(gratarr)
                aeind = np.where( np.array(gratarr) == ag )[0][0]
                print(aeind)
                outroot = '%s-%s-%s-%s'%(program, maskname, id, gratarr[aeind])
                saved = autoextract1d(id, wlarr[aeind], specarr[aeind], errarr[aeind], hdrarr[aeind], path2d, path1d, outroot, gratarr[aeind], manypos, manyfwhm, manymin, manymax)

                extracted[aeind] = saved
            print('Automatic 1d extraction finished.')
            reset2d()

        lastpressed = event.key

    #keypresses2d = fig.canvas.mpl_connect('key_press_event', onkeypress2d)

    #pl.show()

    def reset2d():
        global vmarr, origvmarr, im2darr, obswlarr, origlims, fig, ax, restwlax, ax2drat, pix2wl, wl2pix, keypresses2d, yexparr, plotyexp
        fig = pl.figure(1, figsize=(16, 9))
        fig.suptitle('ID = %s'%id)

        top = 0.9
        bottom = 0.1
        left = 0.05
        right = 0.95

        hspace = 0.1
        height = (top-bottom-float(ngrat-1)*hspace)/float(ngrat)
        ax2drat = height/(right - left)

        #ax = [fig.add_axes([left, bottom+gg*(height+hspace), right-left, height]) for gg in range(ngrat)]
        #ax = ax[::-1]
        #ax = np.array(ax)

        ax = []
        for gg in range(ngrat):
            aa = fig.add_axes([left, bottom+gg*(height+hspace), right-left, height])
            ax += [aa]
        ax = ax[::-1]
        ax = np.array(ax)

        vmarr = []
        origvmarr = []
        im2darr = []
        obswlarr = []

        def pix2wl(i, z, norm):
            def func(pix): return (crval1arr[i]+pix*cdelt1arr[i])/(1.+z)*norm
            return func
        def wl2pix(i, z, norm):
            def func(wl): return (wl/norm*(1.+z)-crval1arr[i])/cdelt1arr[i]
            return func

        yexparr = []
        for ii in range(len(ax)):
            wl = wlarr[ii]
            spec = specarr[ii]
            err = errarr[ii]
            hdr = hdrarr[ii]
            ytrace = hdr['YTRACE']

            # 2D spectrum
            vm = np.sort(spec.flatten()[spec.flatten() != 0.])[int(0.97*len(spec.flatten()[spec.flatten() != 0.]))]
            vmarr += [vm]
            origvmarr += [vm]
            im2darr += [ax[ii].imshow(spec, vmin=-vm, vmax=vm, cmap='viridis', origin='lower')]
            #im2darr += [ax[ii].imshow(spec, vmin=-vm, vmax=vm, cmap='gist_gray', origin='lower')]
            #im2d = ax[ii].imshow(spec, vmin=-vm, vmax=vm, cmap='gist_gray', origin='lower')
            #print(im2d)
            #im2darr = np.append(im2darr, im2d)
            yexparr += [ax[ii].axhline(ytrace, color='orange', linestyle=':', linewidth=1.0)]
            plotyexp = True

            if extracted[ii]: estr = 'EXTRACTION SAVED'
            else: estr = ''
            ax[ii].set_title('%s     %s'%(gratarr[ii], estr))

            #sax2d = ax[ii].secondary_xaxis('top', functions=(lambda pix: crval1arr[ii]+pix*cdelt1arr[ii], lambda wl: (wl-crval1arr[ii])/cdelt1arr[ii]))
            sax2d = ax[ii].secondary_xaxis('top', functions=(pix2wl(ii, 0., 1.), wl2pix(ii, 0., 1.)))
            sax2d.set_xlabel(r'$\lambda_{\mathrm{obs}}$ ($\mu$m)', fontsize=12, labelpad=3)
            #ax2d.tick_params(bottom=False, labelbottom=False)
            obswlarr += [sax2d]

        vmarr = np.array(vmarr)
        im2darr = np.array(im2darr)
        origvmarr = np.array(origvmarr)

        origlims = np.array([aa.axis() for aa in ax])
        restwlax = [False for aa in range(len(ax))]

        keypresses2d = fig.canvas.mpl_connect('key_press_event', onkeypress2d)

        pl.show()

    reset2d()

def extract1d(id, wl, spec, err, shdr, path2d, path1d, outroot, grating):

    print('\nBeginning 1D extraction for target %s in the %s configuration'%(id,grating))
    print('Press "h" or "?" for help.')

    #suffix = '_1d'
    #outname = outroot + suffix + '.txt'

    global ax2d, ax1d, axex
    global ax2d, ax1d, axex, origexlims, extracted
    global ix, iy
    global spat, wspat, spaterr
    global coords, xex1, xex2, xex1temp, vspans
    global ax2d, ax1d, axex, orig2dlims, orig1dlims, origexlims
    global spat, wspat, spaterr
    global lastpressed
    global spec1d, var1d, err1d, spec1dbox, var1dbox, err1dbox
    global fitcomplete, fxlow, fxhigh, lfit, lerr
    global mouseclicks, mousereleases, keypresses
    global zcatarr, lineplotarr
    global vm, im2d
    global restwlax1d, rax1d, wlax2dflag, rax2d
    global yexp, yexpex, plotyexp
    global zax, zxlow, zxhigh, zylow, zyhigh, zoomcomplete
    global saved
    global ypos, yposerr, yfwhm, yfwhmerr, exsnr, exymin, exymax
    global mex1, mex2, mspans, mex1temp, lastclicked
    global excomplete, exfitlow, exfithigh

    spec *= shdr['TOMUJY']
    err *= shdr['TOMUJY']

    coords = []
    xex1 = []
    xex2 = []
    vspans = []
    xex1temp = None

    mex1 = []
    mex2 = []
    mspans = []
    mex1temp = None
    lastclicked = None

    def gauss(x,u, sig, flux): return flux*np.exp(-((u-x)**2./(2.*sig**2.)))/np.sqrt(2*np.pi*sig**2.)
    # line with a constant continuum
    def gaussline(x, u, sig, flux, a): return gauss(x, u, sig, flux) + a

    def uJy_to_cgs(wl): return (1.e-23)*(3.e14)/(wl**2.)*1.e-6   # wl in um

    # check if redshift catalog exists, if not then create it, if so then open it
    zcatname = '00_redshift_catalog.txt'
    if os.path.exists(path1d+zcatname):
        zcat = open(path1d+zcatname, 'a')
    else:
        zcat = open(path1d+zcatname, 'w')
        zcat.write('# id\tline\tlambda0\tlambdaobs\tlambdaobserr\tz\tzerr\tflux\tfluxerr\tgrating\n')

    # check if extraction info catalog exists, if not then create it, if so then open it
    excatname = '00_extract_info.txt'
    if os.path.exists(path1d+excatname):
        excat = open(path1d+excatname, 'a')
    else:
        excat = open(path1d+excatname, 'w')
        excat.write('# id\tgrating\typos\typoserr\tyfwhm\tyfwhmerr\tysnr\tymin\tymax\tautoflag\n')

    extracted1d = False
    saved = False

    #### interactive functions ####
    def onclick(event):
        global ax2d, ax1d, axex
        global spec, err, shdr
        global mex1, mex2, mspans, mex1temp, lastclicked, wl
        global coords, xex1, xex2, xex1temp, extracted1d
        if event.inaxes == ax2d and event.button == MouseButton.LEFT:
            global ix, iy
            ix, iy = event.xdata, event.ydata
            #print('x = %.3f, y = %.3f'%(ix, iy))

            coords.append((ix, iy))
            #xex1.append(int(ix))
            xex1temp = int(ix)

        if event.inaxes == ax1d and event.button == MouseButton.LEFT and extracted1d:
            ix, iy = event.xdata, event.ydata
            #mex1.append(ix)
            mex1temp = ix

        lastclicked = event.inaxes

    def onrelease(event):
        global ax2d, ax1d, axex, origexlims, extracted1d
        global spat, wspat, spaterr
        global coords, xex1, xex2, xex1temp, vspans
        global mex1, mex2, mspans, lastclicked, wl

        if event.button is MouseButton.LEFT and event.inaxes == ax2d and lastclicked == ax2d:
            global ix, iy
            ix, iy = event.xdata, event.ydata
            #print('x = %.3f, y = %.3f'%(ix, iy))
            coords.append((ix, iy))
            xex1.append(xex1temp)
            xex2.append(int(ix))

            if xex2[-1] == xex1[-1]: xex2[-1] += 1  # can't have any 0 pixel slices!
            if xex2[-1] < xex1[-1]: tempx = xex1[-1]; xex1[-1] = xex2[-1]; xex2[-1] = tempx  # need xex1 to be the lower bound

            # plot selected ranges in 2D spectrum
            vspans += [ax2d.axvspan(xex1[-1], xex2[-1], color='r', edgecolor='none', alpha=0.3, zorder=1000)]
            pl.draw()

            ### plot collapsed spatial profile from selected regions
            #axex.clear(); axex.set_title('extraction profile')
            #ax1d.clear(); ax1d.set_title('extracted 1D spectrum')
            resetex(); reset1d()
            extracted1d = False
            # collapse the selected x ranges
            spat = np.zeros(len(spec))
            wspat = np.zeros(len(spec))
            spaterr = np.zeros(len(spec))
            for i in range(len(xex1)):
                spat += np.mean(spec[:,xex1[i]:xex2[i]], axis=1)
                wspat += np.ma.average(spec[:,xex1[i]:xex2[i]], weights=err[:,xex1[i]:xex2[i]]**-2., axis=1) #*float(xex2[i]-xex1[i])
                spaterr += np.sqrt( np.sum(err[:,xex1[i]:xex2[i]]**2., axis=1) )

            axex.plot(np.arange(len(spat)), spat, 'k-', drawstyle='steps-mid')
            axex.plot(np.arange(len(spat)), wspat, 'm-', drawstyle='steps-mid')
            #axex.axvline(yguess, color='orange', linestyle='--')

            axex.text(0.05, 0.95, 'extract x low = %s'%xex1, transform=axex.transAxes)
            axex.text(0.05, 0.9, 'extract x high = %s'%xex2, transform=axex.transAxes)
            ylim = axex.get_ylim()
            axex.set_ylim([ylim[0], 1.2*ylim[1]])
            origexlims = axex.axis()
            pl.draw()

        if event.inaxes == ax1d and event.button == MouseButton.LEFT and lastclicked == ax1d and extracted1d:
            ix, iy = event.xdata, event.ydata
            mex1.append(mex1temp)
            mex2.append(ix)
            mspans += [ax1d.axvspan(mex1[-1], mex2[-1], color='b', edgecolor='none', alpha=0.3, zorder=1000)]
            pl.draw()

        lastclicked = None

    lastpressed = None

    zax = None
    zxlow = None
    zxhigh = None
    zylow = None
    zyhigh = None
    zoomcomplete = True

    excomplete = True
    exfitlow = None
    exfithigh = None

    fitcomplete = True
    lfit = None
    lerr = None
    fxlow = None
    fxhigh = None

    restwlax1d = False
    wlax2dflag = 0  # 0 if pixels, 1 if observed wl, 2 if rest wl

    plotyexp = True

    bc = None
    boxcar = False

    spec1d = []
    var1d = []
    spec1dbox = []
    var1dbox = []

    zcatarr = []
    def printzcat():
        if len(zcatarr) > 0:
            print('\nCurrent redshift catalog for this object:')
            print(('line\t\tlambda0\t\tlambdaobs\terr\tzline\tzlineerr\tgrating'))
            [print('%-13s\t%.3f\t%.3f\t%.3f\t%.5f\t%.5f\t%s'%(zz[0],zz[1],zz[2],zz[3],zz[4],zz[5],zz[8])) for zz in zcatarr]
            print('Mean redshift = %.5f'%(np.average([zz[4] for zz in zcatarr], weights=[zz[5]**-2. for zz in zcatarr])))
        else:
            print('\nRedshift catalog is empty!')

    # write redshift catalog entry
    def writezcat():
        if len(zcatarr) > 0:
            for zz in zcatarr:
                zcat.write('%s\t%s\t%.3f\t%.3f\t%.3f\t%.5f\t%.5f\t%.3e\t%.3e\t%s\n'%(id,*zz))

    def zmean():
        if len(zcatarr) > 0: return np.average([zz[4] for zz in zcatarr], weights=[zz[5]**-2. for zz in zcatarr])
        else: return 0.
    lineplotarr = []

    # write extraction catalog entry
    def writeexcat():
        if len(excatarr) > 0:
            for ee in excatarr:
                excat.write('%s\t%s\t%.3f\t%.3f\t%.3f\t%.5f\t%.5f\t%.3e\t%.3e\n'%(id,*ee))

    def makemask1d(wl, mex1, mex2):
        mask1d = []
        for ww in wl:
            maskval = 0
            for mm in range(len(mex1)):
                if mex1[mm] < ww < mex2[mm]: maskval = 1
            mask1d += [maskval]
        mask1d = np.array(mask1d)
        return mask1d

    def writespec1d(wl, spec1d, err1d, spec1dbox, err1dbox, mex1, mex2):
        wlcgs = wl*1.e4
        cgscoeff = (1.e-23)*(3.e18)/(wlcgs**2.)*1.e-6
        spec1dcgs = spec1d*cgscoeff
        err1dcgs = err1d*cgscoeff
        spec1dboxcgs = spec1dbox*cgscoeff
        err1dboxcgs = err1dbox*cgscoeff

        mask1d = makemask1d(wl, mex1, mex2)
        mask1d[np.logical_not(np.isfinite(err1d))] = 1

        aoutname = path1d + outroot + '_1d.txt'
        foutname = path1d + outroot + '_1d.fits'

        # ascii output
        with open(aoutname, 'w') as outfile:
            outfile.write('lambda\tlambda_cgs\tflux\terr\tflux_cgs\terr_cgs\tflux_box\terr_box\tflux_box_cgs\terr_box_cgs\tmask\n')
            [outfile.write('%.7f\t%.2f\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.0f\n'%(wl[i], wlcgs[i], spec1d[i], err1d[i], spec1dcgs[i], err1dcgs[i], spec1dbox[i], err1dbox[i], spec1dboxcgs[i], err1dboxcgs[i], mask1d[i])) for i in range(len(wl))]

        # fits output file
        cwl = pf.Column(name='lambda', format='D', array=wl, unit='um')
        cwlcgs = pf.Column(name='lambda_cgs', format='D', array=wlcgs, unit='Angstrom')
        cflux = pf.Column(name='flux', format='D', array=spec1d, unit='uJy')
        cerr = pf.Column(name='err', format='D', array=err1d, unit='uJy')
        cfluxbox = pf.Column(name='flux_box', format='D', array=spec1dbox, unit='uJy')
        cerrbox = pf.Column(name='err_box', format='D', array=err1dbox, unit='uJy')
        cfluxcgs = pf.Column(name='flux_cgs', format='D', array=spec1dcgs, unit='erg/s/cm2/AA')
        cerrcgs = pf.Column(name='err_cgs', format='D', array=err1dcgs, unit='erg/s/cm2/AA')
        cfluxboxcgs = pf.Column(name='flux_box_cgs', format='D', array=spec1dboxcgs, unit='erg/s/cm2/AA')
        cerrboxcgs = pf.Column(name='err_box_cgs', format='D', array=err1dboxcgs, unit='erg/s/cm2/AA')
        cmask = pf.Column(name='mask', format='K', array=mask1d, unit='none')
        colarr = [cwl, cwlcgs, cflux, cerr, cfluxbox, cerrbox, cfluxcgs, cerrcgs, cfluxboxcgs, cerrboxcgs, cmask]
        coldefs = pf.ColDefs(colarr)
        outhdu = pf.BinTableHDU.from_columns(coldefs)
        outhdr = outhdu.header
        outhdr.set('ID', id, 'target ID', before='TTYPE1')
        program = outroot.split('-')[0]
        maskname = outroot.split('-')[1]
        outhdr.set('PROGRAM', program, 'program name', before='TTYPE1')
        outhdr.set('MASKNAME', maskname, 'maskname', before='TTYPE1')
        outhdr.set('YSCALE', shdr['CDELT2'], 'cross-dispersion scale in fraction of microshutter (0.53") per pixel')
        outhdr.set('YSCALEAS', shdr['CDELT2']*0.53, 'cross-dispersion in arcseconds per pixel')
        #outhdr.set('DETECTORS', detectors, 'detectors', before='TTYPE1')
        #[outhdr.set('FILE2D'+ff, , 'target ID', before='TTYPE1') for ff in ]

# key words to copy from 2D shdr
# 'NCOMBINE', 'EXPTIME', 'OTHRESH', 'TOMUJY', 'PROFCEN', 'PROFSIG', 'PROFSTRT', 'PROFSTOP', 'YTRACE', 'FILTER', 'GRATING'
        skeys = ['NCOMBINE', 'EXPTIME', 'OTHRESH', 'TOMUJY', 'PROFCEN', 'PROFSIG', 'PROFSTRT', 'PROFSTOP', 'YTRACE', 'FILTER', 'GRATING']
        for ss in skeys:
            outhdr.set(ss, shdr[ss], shdr.comments[ss], before='TTYPE1')

        outtime = gmtime()
        outhdr.set('EXDATE', '%02.f-%02.f-%04.f'%(outtime[2], outtime[1], outtime[0]), 'extraction date UTC', before='TTYPE1')
        outhdr.set('EXTIME', '%02.f:%02.f:%02.f'%(outtime[3], outtime[4], outtime[5]), 'extraction time UTC', before='TTYPE1')
        outhdr.set('EXYPOS', ypos, 'extracted y position', before='TTYPE1')
        outhdr.set('EXYPOSERR', yposerr, 'extracted y position error', before='TTYPE1')
        outhdr.set('EXYFWHM', yfwhm, 'extracted y fwhm', before='TTYPE1')
        outhdr.set('EXYFWHMERR', yfwhmerr, 'extracted y fwhm error', before='TTYPE1')
        outhdr.set('EXSNR', exsnr, 'extraction gaussian SNR', before='TTYPE1')
        outhdr.set('EXYMIN', exymin, 'minimum y of extraction window', before='TTYPE1')
        outhdr.set('EXYMAX', exymax, 'maximum y of extraction window', before='TTYPE1')
        outhdr.set('EXTYPE', 'manual', 'manual or auto extraction', before='TTYPE1')
        for xx in range(len(xex1)):
            outhdr.set('XLOW%s'%xx, xex1[xx], 'low x-pixel for extraction profile construction', before='TTYPE1')
            outhdr.set('XHIGH%s'%xx, xex2[xx], 'high x-pixel for extraction profile construction', before='TTYPE1')
        outhdr.set('EXWIDTH', float(exymax - exymin + 1)*shdr['CDELT2']*0.53, 'width of extraction window in arcseconds', before='TTYPE1')
        outhdr.set('EXYUP', float(exymax + 0.5) - ypos, 'pixels above trace for extraction window', before='TTYPE1')
        outhdr.set('EXYDOWN', ypos - float(exymin - 0.5), 'pixels below trace for extraction window', before='TTYPE1')
        outhdr.set('EXYUPAS', (float(exymax + 0.5) - ypos)*shdr['CDELT2']*0.53, 'arcseconds above trace for extraction window', before='TTYPE1')
        outhdr.set('EXYDOWNAS', (ypos - float(exymin - 0.5))*shdr['CDELT2']*0.53, 'arcseconds below trace for extraction window', before='TTYPE1')
        if len(zcatarr) > 0: zm = np.average([zz[4] for zz in zcatarr], weights=[zz[5]**-2. for zz in zcatarr])
        else: zm = -1.
        outhdr.set('ZINIT', zm, 'mean redshift from lines', before='TTYPE1')

        #print(outhdu.header)
        outhdu.writeto(foutname, overwrite=True)
        print('1D spectrum saved to:')
        print('%s'%aoutname)
        print('%s'%foutname)
        print('\n')


    def onkeypress(event):
        global coords, xex1, xex2, vspans
        global ax2d, ax1d, axex, orig2dlims, orig1dlims, origexlims
        global spat, wspat, spaterr
        global lastpressed
        global spec1d, var1d, err1d, spec1dbox, var1dbox, err1dbox, extracted1d, saved
        global fitcomplete, fxlow, fxhigh, lfit, lerr
        global excomplete, exfitlow, exfithigh
        global mouseclicks, mousereleases, keypresses
        global zcatarr, lineplotarr
        global vm, im2d
        global restwlax1d, rax1d, wlax2dflag, rax2d
        global yexp, yexpex, plotyexp
        global ypos, yposerr, yfwhm, yfwhmerr, exsnr, exymin, exymax
        global mex1, mex2, mspans
        global bc, boxcar

        #print('press', event.key)

        # reset the extraction profile and redshift catalog and extraction
        if event.key == 'r':
            extracted1d = False
            coords = []
            xex1 = []
            xex2 = []
            if len(vspans) > 0:
                [vv.remove() for vv in vspans]
                #axex.clear(); axex.set_title('extraction profile')
                #ax1d.clear(); ax1d.set_title('extracted 1D spectrum')
                resetex(); reset1d()
                pl.draw()
                vspans = []
            zcatarr = []

        # remove last entry from 1D mask catalog
        if event.key == 'u':
            if len(mspans) > 0:
                mspans[-1].remove()
                mspans = mspans[:-1]
                mex1 = mex1[:-1]
                mex2 = mex2[:-1]
                pl.draw()

        # toggle boxcar extraction in 1D plot
        if event.key == 'b' and extracted1d:
            if boxcar: bc.remove(); pl.draw(); boxcar = False
            else: ax1d.add_artist(bc); pl.draw(); boxcar = True

        # print the extraction windows, redshift catalog
        if event.key == 'p':
            print('\n')
            print('extraction windows x low = %s'%xex1)
            print('extraction windows x high = %s'%xex2)
            print('1D mask windows lambda_obs low (um) = %s'%[round(mm, 5) for mm in mex1])
            print('1D mask windows lambda_obs high (um) = %s'%[round(mm, 5) for mm in mex2])
            printzcat()

        # toggle expected y position of main target on and off
        if event.key == 'y':
            if plotyexp: yexp.remove(); yexpex.remove(); pl.draw(); plotyexp = False
            else: ax2d.add_line(yexp); axex.add_line(yexpex); pl.draw(); plotyexp = True


        # fit a Gaussian profile to the collapsed spatial profile
        if event.key == 'e':
            if event.inaxes != axex: excomplete = True
            else:
                if lastpressed != 'e' or (lastpressed == 'e' and excomplete == True):  # set lower bounds
                    exfitlow = event.xdata
                    excomplete = False
                elif lastpressed == 'e' and excomplete == False:  # set higher bounds and complete fit

                    if event.inaxes == axex and extracted1d == False:
                        extracted1d = True
                        #ax1d.clear(); ax1d.set_title('extracted 1D spectrum')
                        reset1d()
                        exfithigh = event.xdata
                        #yguess = (exfithigh + exfitlow)/2.
                        yguess = ytrace
                        ampguess = np.sum(spat[int(yguess-2):int(yguess+2)])
                        print('Fitting extraction profile with y position guess = %.2f'%yguess)
                        p0 = [yguess, 2., ampguess]
                        #lowbounds = [yguess - 5.5, 1., -np.inf]
                        #highbounds = [yguess + 5.5, np.inf, np.inf]
                        lowbounds = [yguess - 9, 1., -np.inf]
                        highbounds = [yguess + 9, np.inf, np.inf]
                        exfitarr = np.sort([exfitlow, exfithigh])
                        exfitind = np.where( (np.arange(len(spat)) > exfitarr[0]) & (np.arange(len(spat)) < exfitarr[1]) )[0]
                        exfityrange = np.arange(len(spat))[exfitind][spat[exfitind] >= 0.]
                        exfitspat = spat[exfitind][spat[exfitind] >= 0.]
                        pfit, pcov = curve_fit(gauss, exfityrange, exfitspat, p0=p0, bounds=(lowbounds, highbounds))
                        #pfit, pcov = curve_fit(gauss, np.arange(len(spat))[spat>=0.], spat[spat>=0.], p0=p0)
                        #pfit, pcov = curve_fit(gauss, np.arange(len(spat))[spat>=0.], spat[spat>=0.], sigma=spaterr[spat>=0.], p0=p0)
                        pfit = np.abs(pfit)
                        perr = np.sqrt(np.diag(pcov))
                        perr[np.logical_not(np.isfinite(perr))] = -1.
                        ypos = pfit[0]
                        yposerr = perr[0]
                        yfwhm = pfit[1]*2.355
                        yfwhmerr = perr[1]*2.355
                        exsnr = pfit[2]/perr[2]
                        if np.logical_not(np.isfinite(exsnr)): exsnr = -1.
                        print('Fit parameters: y center=%.2f\ty sigma=%.2f\tSNR=%.2f'%(pfit[0],pfit[1],pfit[2]/perr[2]))
                        plrange = np.linspace(np.amin(np.arange(len(spat))), np.amax(np.arange(len(spat))), 200)
                        axex.plot(plrange, gauss(plrange, *pfit), 'r-')
                        #axex.axvline(yguess, color='orange', linestyle='--')
                        axex.axvline(pfit[0], color='r', linestyle='--')

                        # select only the y-pixels that are part of the positive trace
                        exycent = pfit[0]
                        exycenti = int(exycent)
                        exyinds = [exycenti]
                        i = exycenti+1
                        while spat[i] > 0.:
                            exyinds += [i]
                            i += 1
                        i = exycenti-1
                        while spat[i] > 0.:
                            exyinds += [i]
                            i -= 1
                        exyinds = np.sort(exyinds)
                        #print(exyinds)

                        exymin = np.amin(exyinds)
                        exymax = np.amax(exyinds)

                        axex.plot(np.arange(len(spat))[exyinds], spat[exyinds], 'b-')
                        axex.plot(plrange[(plrange >= min(exyinds)) & (plrange <= max(exyinds))], gauss(plrange[(plrange >= min(exyinds)) & (plrange <= max(exyinds))], pfit[0], pfit[1], pfit[2]), 'c-')

                        axex.text(0.05, 0.85, 'fit y center = %.2f'%exycent, transform=axex.transAxes)
                        axex.text(0.05, 0.8, 'fit y fwhm = %.2f'%(pfit[1]*2.355), transform=axex.transAxes)
                        axex.text(0.05, 0.75, 'fit SNR = %.2f'%(pfit[2]/perr[2]), transform=axex.transAxes)

                        spec1d = []
                        var1d = []
                        spec1dbox = []
                        var1dbox = []
                        for xx in range(nx):
                            spec1d += [np.sum( gauss(exyinds, pfit[0], pfit[1], 1.)*spec[exyinds,xx]/(err[exyinds,xx]**2.) )/np.sum( (gauss(exyinds, pfit[0], pfit[1], 1.)**2.)/(err[exyinds,xx]**2.) )]
                            var1d += [np.sum( gauss(exyinds, pfit[0], pfit[1], 1.) )/np.sum( (gauss(exyinds, pfit[0], pfit[1], 1.)**2.)/(err[exyinds,xx]**2.) )]

                            spec1dbox += [np.sum(spec[exyinds,xx])]
                            var1dbox += [np.sum(err[exyinds,xx]**2.)]

                        spec1d = np.array(spec1d)
                        var1d = np.array(var1d)
                        err1d = np.sqrt(var1d)
                        spec1dbox = np.array(spec1dbox)
                        var1dbox = np.array(var1dbox)
                        err1dbox = np.sqrt(var1dbox)

                        spec1dbox[np.isnan(spec1d)] = np.nan

                        [ax1d.axvspan(wl[xex1[ii]], wl[xex2[ii]], color='r', alpha=0.15, zorder=-5) for ii in range(len(xex1))]
                        ax1d.plot(wl, spec1d, 'k-', drawstyle='steps-mid', label='optimal')
                        bc, = ax1d.plot(wl, spec1dbox, 'g-', drawstyle='steps-mid', label='boxcar', zorder=-1)
                        bc.remove(); boxcar = False
                        ax1d.fill_between(wl, 0., err1d, color='0.7')
                        #ax1d.fill_between(wl, 0., err1dbox, color='r', alpha=0.3)
                        #ax1d.legend(loc='upper right')
                        orig1dlims = ax1d.axis()

                        pl.draw()

        # fit an emission line with a Gaussian profile
        if event.key == 'f':
            if event.inaxes != ax1d: fitcomplete = True
            else:
                if lastpressed != 'f' or (lastpressed == 'f' and fitcomplete == True):  # set lower bounds
                    fxlow = event.xdata
                    fitcomplete = False
                elif lastpressed == 'f' and fitcomplete == False:  # set higher bounds and complete fit
                    fxhigh = event.xdata
                    fxlow, fxhigh = np.sort([fxlow, fxhigh])
                    ax1d.set_xlim(np.sort([fxlow, fxhigh]))
                    print('\nFitting a line between %.0f AA and %.0f AA'%(fxlow*1.e4,fxhigh*1.e4))

                    fitinds = np.where( (wl >= fxlow) & (wl <= fxhigh) & np.isfinite(spec1d) & (np.isfinite(err1d)) )[0]

                    contguess = np.median(spec1d[fitinds])
                    fluxguess = np.sum(spec1d[fitinds] - contguess)
                    p0 = [wl[fitinds][np.argmax(spec1d[fitinds])], 50./3.e5*np.mean(wl[fitinds]), fluxguess, contguess]
                    lfit, lcov = curve_fit(gaussline, wl[fitinds], spec1d[fitinds], sigma=err1d[fitinds], p0=p0)
                    lerr = np.sqrt(np.diag(lcov))
                    #print(lfit)
                    print('Fit parameters:\tcent=%.1f AA\tsigma=%.2f AA\tSNR=%.2f'%(lfit[0]*1.e4,lfit[1]*1.e4,lfit[2]/lerr[2]))
                    frange = np.linspace(wl[fitinds[0]], wl[fitinds[-1]], 500)
                    fc, = ax1d.plot(frange, [lfit[-1] for ww in frange], 'b--')
                    fg, = ax1d.plot(frange, gaussline(frange, *lfit), 'c-')
                    fu = ax1d.axvline(lfit[0], color='c', linestyle=':')
                    fparr = [fc, fg, fu]

                    fymin = np.amin(gaussline(frange, *lfit))
                    fymax = np.amax(gaussline(frange, *lfit))
                    ax1d.set_ylim(np.sort([fymin - 0.40*(fymax - fymin), fymax + 0.20*(fymax - fymin)]))
                    pl.draw()
                    fitcomplete = True

                    # oii, neiii, hd, hg, hb, oiii4959, oiii5007, ha, nii6585, sii6718, sii6731, siii9069, siii9531
                    linenamearr = ['Lyalpha', 'CIII]1908', '[OII]doublet', '[NeIII]3869', 'Hdelta', 'Hgamma', 'Hbeta', '[OIII]4959', '[OIII]5007', '[NII]6548', 'Halpha', '[NII]6584', '[SII]6716', '[SII]6731', '[SIII]9069', '[SIII]9531', 'HeI10830', 'Pabeta', 'Paalpha']
                    linecent0arr = [1215.67, round((1906.68+1908.73)/2.,3), round((3727.092+3729.875)/2.,3), 3869.856, 4102.89, 4341.68, 4862.68, 4960.295, 5008.240, 6549.86, 6564.61, 6585.27, 6718.29, 6732.67, 9071.09, 9533.21, 10833.305, 12821.578, 18756.420]
                    linenumarr = ['%s'%ll for ll in range(len(linenamearr))]
                    linestrarr = ['%.0f:\t%-13s\tcent0 = %s AA\tz = %.5f\n'%(ll, linenamearr[ll], linecent0arr[ll], lfit[0]*1.e4/linecent0arr[ll]-1.) for ll in range(len(linenamearr))]
                    linestr = ''.join(linestrarr)

            # UV lines
            #linecents = [1548.19+1., 1550.77+1., 1640.42, 1660.81, 1666.15, 1749.67, 1752.16, 1882.47, 1892.03, 1906.68, 1908.73, 2321.66, 2325.40, 2326.93, 2328.12, 2471.03]
            #linenames = ['civ1548', 'civ1550', 'heii1640', 'oiii1661', 'oiii1666', 'niii1750', 'niii1752', 'siiii1882', 'siiii1892', 'ciii1907', 'ciii1909', 'oiii2322', 'cii2325', 'cii2327', 'cii2328', 'oii2471']
            # optical lines
            #linecents += [(3727.092+3729.875)/2., 3869.856, (3890.166+3889.750)/2., (3971.198+3968.593)/2., 4077.500, 4102.89, 4341.68, 4364.436, 4687.021, 4862.68, 4960.295, 5008.240, 5756.24, 5877.25, 6302.046, 6313.806, 6365.535, 6549.86, 6564.61, 6585.27, 6679.995, 6718.29, 6732.67, 7067.138, 7137.767, 7321.937, 7332.210, 9071.09, 9533.21]
            #linenames += ['oii3727', 'neiii3869', 'hihei3888', 'hineiii3970', 'sii4076', 'hd', 'hg', 'oiii4363', 'heii4686', 'hb', 'oiii4959', 'oiii5007', 'nii5755', 'hei5876', 'oi6300', 'siii6312', 'oi6363', 'nii6548', 'ha', 'nii6584', 'hei6680', 'sii6716', 'sii6731', 'hei7067', 'ariii7135', 'oii7320', 'oii7330', 'siii9069', 'siii9531']

                    # disconnect interactive functions
                    efig.canvas.mpl_disconnect(mouseclicks)
                    efig.canvas.mpl_disconnect(mousereleases)
                    efig.canvas.mpl_disconnect(keypresses)

                    num = input('\nIdentify the line:\n\n' + linestr + 'c: Custom line name/rest wavelength\nEnter selection: ')

                    if num == 'c':  # enter custom line name/rest wavelength
                        linename = input('Enter line name (e.g., CIII]1908): ')
                        linecent0 = float(input('Enter rest-frame wavelength in Angstroms: '))
                        linecent = lfit[0]*1.e4
                        linez = linecent/linecent0 - 1.
                        linecenterr = lerr[0]*1.e4; linezerr = linecenterr/linecent0
                        lineflux = lfit[2]*uJy_to_cgs(lfit[0]); linefluxerr = lerr[2]*uJy_to_cgs(lfit[0])
                        print('Line identified as %s: z = %.5f'%(linename, linez))
                        print('Saving line %s to redshift catalog.'%linename)
                        zcatarr += [[linename, linecent0, linecent, linecenterr, linez, linezerr, lineflux, linefluxerr, grating]]
                        printzcat()
                        lineplotarr += [fparr]
                    elif num not in linenumarr:
                        print('Invalid selection.  Exiting without saving line fit.')
                        [ff.remove() for ff in fparr]
                        pl.draw()
                    else:
                        lineplotarr += [fparr]
                        num = int(num)
                        linecent = lfit[0]*1.e4; linecent0 = linecent0arr[num]; linename = linenamearr[num]
                        linez = linecent/linecent0 - 1.
                        linecenterr = lerr[0]*1.e4; linezerr = linecenterr/linecent0
                        lineflux = lfit[2]*uJy_to_cgs(lfit[0]); linefluxerr = lerr[2]*uJy_to_cgs(lfit[0])
                        print('Line identified as %s: z = %.5f'%(linename, linez))
                        print('Saving line %s to redshift catalog.'%linename)
                        zcatarr += [[linename, linecent0, linecent, linecenterr, linez, linezerr, lineflux, linefluxerr, grating]]
                        printzcat()

                    # reconnect interactive functions
                    mouseclicks = efig.canvas.mpl_connect('button_press_event', onclick)
                    mousereleases = efig.canvas.mpl_connect('button_release_event', onrelease)
                    keypresses = efig.canvas.mpl_connect('key_press_event', onkeypress)

        # modify redshift catalog
        if event.key == 'm':
            # disconnect interactive functions
            efig.canvas.mpl_disconnect(mouseclicks)
            efig.canvas.mpl_disconnect(mousereleases)
            efig.canvas.mpl_disconnect(keypresses)

            print('\nCurrent redshift catalog for this object:')
            print(('\tline\t\tlambda0\t\tlambdaobs\terr\tzline\tzlineerr'))
            [print('%s\t%-13s\t%.3f\t%.3f\t%.3f\t%.5f\t%.5f'%(zz,zcatarr[zz][0],zcatarr[zz][1],zcatarr[zz][2],zcatarr[zz][3],zcatarr[zz][4],zcatarr[zz][5])) for zz in range(len(zcatarr))]
            num = input('Select a number to remove from the redshift catalog: ')
            if num not in ['%s'%rr for rr in range(len(zcatarr))]:
                print('Invalid selection.  Exiting without modifying the redshift catalog.')
            else:
                num = int(num)
                zcatarr = zcatarr[:num] + zcatarr[num+1:]
                [ff.remove() for ff in lineplotarr[num]]
                lineplotarr = lineplotarr[:num] + lineplotarr[num+1:]
                pl.draw()
                printzcat()

            # reconnect interactive functions
            mouseclicks = efig.canvas.mpl_connect('button_press_event', onclick)
            mousereleases = efig.canvas.mpl_connect('button_release_event', onrelease)
            keypresses = efig.canvas.mpl_connect('key_press_event', onkeypress)

        # toggle rest wavelength in ax1d based on z_mean from current redshift catalog
        if event.key == 'w' and event.inaxes == ax1d:
            if len(zcatarr) > 0 and not restwlax1d:
                zm = zmean()
                rax1d = ax1d.secondary_xaxis('bottom', functions=(lambda wl: wl*1.e4/(1.+zm), lambda wl: wl/1.e4*(1.+zm)))
                rax1d.set_xlabel(r'$\lambda_{\mathrm{rest}}$ ($\AA$)', fontsize=12, labelpad=-32)
                rax1d.tick_params(axis='x', direction='in', pad=-15)
                restwlax1d = True
                pl.draw()
            elif restwlax1d:
                rax1d.remove()
                restwlax1d = False
                pl.draw()

        # toggle from pixels to observed wl to rest wl in ax2d
        if event.key == 'w' and event.inaxes == ax2d:
            if wlax2dflag == 0:  # convert to obs
                rax2d = ax2d.secondary_xaxis('bottom', functions=(lambda pix: shdr['crval1']+pix*shdr['cdelt1'], lambda wl: (wl-shdr['crval1'])/shdr['cdelt1']))
                rax2d.set_xlabel(r'$\lambda_{\mathrm{obs}}$ ($\mu$m)', fontsize=12, labelpad=0)
                ax2d.tick_params(bottom=False, labelbottom=False)
                wlax2dflag = 1
                pl.draw()
            elif len(zcatarr) == 0 and wlax2dflag == 1:  # convert from obs to pix because zmean isn't defined
                rax2d.remove()
                ax2d.tick_params(bottom=True, labelbottom=True)
                wlax2dflag = 0
                pl.draw()
            elif len(zcatarr) > 0 and wlax2dflag == 1:  # convert from obs to rest based on zmean
                rax2d.remove()
                zm = zmean()
                rax2d = ax2d.secondary_xaxis('bottom', functions=(lambda pix: (shdr['crval1']+pix*shdr['cdelt1'])*1.e4/(1.+zm), lambda wl: (wl*(1.+zm)/1.e4-shdr['crval1'])/shdr['cdelt1']))
                rax2d.set_xlabel(r'$\lambda_{\mathrm{rest}}$ ($\AA$)', fontsize=12, labelpad=0)
                ax2d.tick_params(bottom=False, labelbottom=False)
                wlax2dflag = 2
                pl.draw()
            elif wlax2dflag == 2:  # convert from rest to pix
                rax2d.remove()
                wlax2dflag = 0
                ax2d.tick_params(bottom=True, labelbottom=True)
                pl.draw()

        # zoom in an axis
        if event.key == 'z':
            global zax, zxlow, zxhigh, zylow, zyhigh, zoomcomplete
            if zax != event.inaxes: zoomcomplete = True
            if lastpressed != 'z' or (lastpressed == 'z' and zoomcomplete == True):  # set lower bounds
                zax = event.inaxes
                zxlow, zylow = event.xdata, event.ydata
                zoomcomplete = False
            elif lastpressed == 'z' and zoomcomplete == False:  # set higher bounds
                if event.inaxes == zax:
                    zxhigh, zyhigh = event.xdata, event.ydata
                    if zax != ax2d:
                        zax.set_xlim(np.sort([zxlow, zxhigh]))
                        zax.set_ylim(np.sort([zylow, zyhigh]))  # never change the 2D y range
                    else:
                        zxcent = (zxhigh+zxlow)/2.
                        figsize = efig.get_size_inches()
                        zxrange = float(shdr['naxis2'])/ax2drat*(figsize[0]/figsize[1])
                        event.inaxes.set_xlim([zxcent-zxrange/2., zxcent+zxrange/2.])
                    pl.draw()
                zoomcomplete = True

        # reset axis zoom
        if event.key == 'x':
            if event.inaxes == ax2d: ax2d.axis(orig2dlims); pl.draw()
            if event.inaxes == ax1d: ax1d.axis(orig1dlims); pl.draw()
            if event.inaxes == axex: axex.axis(origexlims); pl.draw()

        if event.key == 'h' or event.key == '?':
            helpstr = '\n In the 2D science spectrum panel, click and drag to add pixels in that x-range to the spatial profile.\n\
       In the extracted 1D spectrum panel, click and drag to add wavelength ranges to the 1D mask.\n\
       Key bindings:\n\
       r -- "reset": Reset the extraction profile, redshift catalog, and 1D extraction.\n\
       p -- "print": Print useful information including the extraction windows and current redshift catalog.\n\
       y -- "yposition": Toggle marker for the expected y-position of the primary target (orange dotted line).\n\
       e -- "extract": With the cursor in the extraction profile panel, press "e" twice to fit a Gaussian profile to the spatial profile.  Hit "e" once on the at the cursor y-position and extract the 1D spectrum.\n\
       b -- "boxcar": Toggle a plot of the boxcar extraction on and off.\n\
       f -- "fit": Fit an emission line with a Gaussian profile.  Press "f" twice, once on each side of the line, to fit the line.  The centroid guess is the highest data point in the selected wavelength range.  Follow the prompt to estimate the redshift and add the line to the redshift catalog.  Line fits currently stored in the redshift catalog will remain plotted in the 1D spectrum panel.\n\
       m -- "modify": Modify the redshift catalog.  This function allows you to remove an entry from the redshift catalog.\n\
       u -- "undo mask": Remove the last masked region from the 1D mask catalog.\n\
       w -- "wavelength": Press "w" in the 2D spectrum panel to switch between pixel and observed wavelength on the x-axis.  If a redshift has been estimated from emission lines, "w" will switch between pixel, observed, and rest-frame wavelength in the 2D spectrum panel.  In the 1D spectrum panel, "w" will toggle rest-frame wavelength tick marks on the x-axis.\n\
       z -- "zoom": Press "z" twice in the same panel to zoom in to a range with corners defined by the cursor positions when "z" was pressed.\n\
       x -- "axis reset": Resets the plotted range to the original range for the axis in which the cursor is.\n\
       c -- "center": Move the point under the cursor to the center of the panel.\n\
       i -- "in": Zoom the selected axis in, keeping the same center coordinate.\n\
       o -- "out": Zoom the selected axis out, keeping the same center coordinate.\n\
       s -- "save": Save the extracted 1D spectrum, extraction information, and redshift catalog.\n\
       q -- "quit": Exit the 1D extraction.  Make sure to save first!\n\
       arrow keys: Use the arrow keys to shift the plot.  The extraction profile and 1D spectrum panels can be shifted up/down/left/right.  The 2D spectrum panel can only be shifted left and right, while up and down changes the stretch.\n'

            print(helpstr)


        # pan center of axis
        if event.key == 'c':
            newcx, newcy = event.xdata, event.ydata
            oldxlow, oldxhigh = event.inaxes.get_xlim()
            oldylow, oldyhigh = event.inaxes.get_ylim()
            oldxrange = (oldxhigh - oldxlow)
            oldyrange = (oldyhigh - oldylow)
            event.inaxes.set_xlim([newcx - oldxrange/2., newcx + oldxrange/2.])
            if event.inaxes != ax2d: event.inaxes.set_ylim([newcy - oldyrange/2., newcy + oldyrange/2.])
            pl.draw()

        # pan axis with arrow keys
        if event.key == 'left':
            oldxlow, oldxhigh = event.inaxes.get_xlim()
            oldylow, oldyhigh = event.inaxes.get_ylim()
            oldxrange = (oldxhigh - oldxlow)
            oldyrange = (oldyhigh - oldylow)
            event.inaxes.axis([oldxlow-0.1*oldxrange, oldxhigh-0.1*oldxrange, oldylow, oldyhigh])
            pl.draw()
        if event.key == 'right':
            oldxlow, oldxhigh = event.inaxes.get_xlim()
            oldylow, oldyhigh = event.inaxes.get_ylim()
            oldxrange = (oldxhigh - oldxlow)
            oldyrange = (oldyhigh - oldylow)
            event.inaxes.axis([oldxlow+0.1*oldxrange, oldxhigh+0.1*oldxrange, oldylow, oldyhigh])
            pl.draw()
        if event.key == 'down' and event.inaxes != ax2d:
            oldxlow, oldxhigh = event.inaxes.get_xlim()
            oldylow, oldyhigh = event.inaxes.get_ylim()
            oldxrange = (oldxhigh - oldxlow)
            oldyrange = (oldyhigh - oldylow)
            event.inaxes.axis([oldxlow, oldxhigh, oldylow-0.1*oldyrange, oldyhigh-0.1*oldyrange])
            pl.draw()
        if event.key == 'up' and event.inaxes != ax2d:
            oldxlow, oldxhigh = event.inaxes.get_xlim()
            oldylow, oldyhigh = event.inaxes.get_ylim()
            oldxrange = (oldxhigh - oldxlow)
            oldyrange = (oldyhigh - oldylow)
            event.inaxes.axis([oldxlow, oldxhigh, oldylow+0.1*oldyrange, oldyhigh+0.1*oldyrange])
            pl.draw()

        # increase or decrease stretch scale of 2D spectrum
        if event.key == 'up' and event.inaxes == ax2d:
            vm *= 2.
            im2d.set_clim(vmin=-vm, vmax=vm)
            pl.draw()
        if event.key == 'down' and event.inaxes == ax2d:
            vm /= 2.
            im2d.set_clim(vmin=-vm, vmax=vm)
            pl.draw()

        # zoom in and out with fixed center
        if event.key == 'i':
            oldxlow, oldxhigh = event.inaxes.get_xlim()
            oldylow, oldyhigh = event.inaxes.get_ylim()
            oldxrange = (oldxhigh - oldxlow)
            oldyrange = (oldyhigh - oldylow)
            oldxcent = (oldxlow + oldxhigh)/2.
            oldycent = (oldylow + oldyhigh)/2.
            if event.inaxes != ax2d: event.inaxes.axis([oldxcent-0.5*oldxrange/2.,oldxcent+0.5*oldxrange/2.,oldycent-0.5*oldyrange/2.,oldycent+0.5*oldyrange/2.])
            else:
                figsize = efig.get_size_inches()
                minxrange = float(shdr['naxis2'])/ax2drat*(figsize[0]/figsize[1])
                newxrange = max([0.5*oldxrange,minxrange])
                event.inaxes.axis([oldxcent-newxrange/2.,oldxcent+newxrange/2.,oldylow,oldyhigh])
            pl.draw()

        # zoom in and out with fixed center
        if event.key == 'o':
            oldxlow, oldxhigh = event.inaxes.get_xlim()
            oldylow, oldyhigh = event.inaxes.get_ylim()
            oldxrange = (oldxhigh - oldxlow)
            oldyrange = (oldyhigh - oldylow)
            oldxcent = (oldxlow + oldxhigh)/2.
            oldycent = (oldylow + oldyhigh)/2.
            if event.inaxes != ax2d: event.inaxes.axis([oldxcent-2.0*oldxrange/2.,oldxcent+2.0*oldxrange/2.,oldycent-2.0*oldyrange/2.,oldycent+2.0*oldyrange/2.])
            else:
                figsize = efig.get_size_inches()
                minxrange = float(shdr['naxis2'])/ax2drat*(figsize[0]/figsize[1])
                newxrange = max([2.0*oldxrange,minxrange])
                event.inaxes.axis([oldxcent-newxrange/2.,oldxcent+newxrange/2.,oldylow,oldyhigh])
            pl.draw()

        # save extracted spectrum, and redshift and extract catalogs
        if event.key == 's' and extracted1d:
            print ('Saving output')
            writezcat()
            writespec1d(wl, spec1d, err1d, spec1dbox, err1dbox, mex1, mex2)
            excat.write('%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.1f\t%.1f\t%.0f\n'%(id, grating, ypos, yposerr, yfwhm, yfwhmerr, exsnr, exymin, exymax, 0))  # write extraction cat
            saved = True

        lastpressed = event.key

    def onkeyrelease(event):
        global coords, xex1, xex2, vspans
        global ax2d, ax1d, axex
        print('release', event.key)

        # zoom in an axis
        if event.key == 'z':
            global zax, zxlow, zxhigh, zylow, zyhigh
            if event.inaxes == zax:
                zxhigh, zyhigh = event.xdata, event.ydata
                zax.set_xlim(np.sort([zxlow, zxhigh]))
                zax.set_ylim(np.sort([zylow, zyhigh]))
                pl.draw()

    #### end interactive functions ####

    nx = shdr['NAXIS1']
    ny = shdr['NAXIS2']
    #ypos = shdr['SRCYPOS']
    #yguess = 20.
    #if ypos < 0.: ypos = 1.+ypos
    #yguess = float(ny) - ypos*float(ny)
    #yguess = float(ny)/2.*(1 + ypos)
    ytrace = shdr['YTRACE']

    wl = shdr['crval1'] + np.arange(shdr['naxis1'])*shdr['cdelt1']

    efig = pl.figure(2, figsize=(16, 9))
    efig.suptitle('id=%s     grating=%s     manual extraction'%(id, grating))


    # initialize axes, size and spacing figure setup
    lbuffer = 0.05
    rbuffer = 0.02
    bbuffer = 0.07
    tbuffer = 0.04
    wspace = 0.05
    hspace = 0.06

    h2d = 0.15

    w1d = 0.6

    ax2d = efig.add_axes([lbuffer, 1.-tbuffer-h2d, 1.-rbuffer-lbuffer, h2d])
    yexp = ax2d.axhline(ytrace, color='orange', linestyle=':', linewidth=1.0)
    ax2drat = h2d/(1.-rbuffer-lbuffer)
    def reset2d():
        ax2d.clear()
        ax2d.set_title('2D science spectrum')
        if plotyexp: ax2d.add_line(yexp)
    reset2d()

    ax1d = efig.add_axes([lbuffer+(1.-rbuffer-w1d-wspace-lbuffer)+wspace, bbuffer, w1d, 1.-tbuffer-h2d-hspace-bbuffer])
    #ax1d.yaxis.get_ticklocs(minor=True)
    #ax1d.minorticks_on()
    #ax1d.tick_params(which='minor', length=3, direction='out')
    #ax1d.tick_params(which='major', length=6)
    def reset1d():
        global mex1, mex2, mspans
        ax1d.clear()
        ax1d.set_title('extracted 1D spectrum')
        ax1d.set_xlabel(r'$\lambda_{\mathrm{obs}}$ ($\mu$m)', fontsize=15)
        ax1d.set_ylabel(r'F$_\nu$ ($\mu$Jy)', fontsize=15)
        mex1 = []
        mex2 = []
        if len(mspans) > 0:
            [mm.remove() for mm in mspans]
            pl.draw()
            mspans = []
    reset1d()

    axex = efig.add_axes([lbuffer, bbuffer, 1.-rbuffer-w1d-wspace-lbuffer, 1.-tbuffer-h2d-hspace-bbuffer])
    yexpex = axex.axvline(ytrace, color='orange', linestyle=':', linewidth=1.0)
    def resetex():
        axex.clear()
        axex.set_title('extraction profile')
        axex.set_xlabel('y pixels', fontsize=15)
        axex.set_ylabel('signal collapsed over x ranges', fontsize=15)
        if plotyexp: axex.add_line(yexpex)
    resetex()

    # 2D spectrum
    vm = np.sort(spec.flatten()[spec.flatten() != 0.])[int(0.97*len(spec.flatten()[spec.flatten() != 0.]))]
    origvm = vm
    im2d = ax2d.imshow(spec, vmin=-vm, vmax=vm, cmap='viridis', origin='lower')

    orig2dlims = ax2d.axis()
    orig1dlims = ax1d.axis()
    origexlims = axex.axis()

    mouseclicks = efig.canvas.mpl_connect('button_press_event', onclick)
    mousereleases = efig.canvas.mpl_connect('button_release_event', onrelease)
    keypresses = efig.canvas.mpl_connect('key_press_event', onkeypress)
    #efig.canvas.mpl_connect('key_release_event', onkeyrelease)

# pl.rcParams['keymap.pan'].remove('p')
# pl.rcParams['keymap.save'].remove('s')
# pl.rcParams['keymap.fullscreen'].remove('f')
# pl.rcParams['keymap.home'].remove('r')
# pl.rcParams['keymap.zoom'].remove('o')

    pl.show()
    pl.pause(10000)

    # disconnect interactive functions
    efig.canvas.mpl_disconnect(mouseclicks)
    efig.canvas.mpl_disconnect(mousereleases)
    efig.canvas.mpl_disconnect(keypresses)

    pl.close(efig)

    zcat.close()
    excat.close()

    #pl.figure(1)
    #print('finished 1d extraction')

    return saved

### automatic (blind) extraction

def autoextract1d(id, wl, spec, err, shdr, path2d, path1d, outroot, grating, ypos, yfwhm, ymin, ymax):

    print('\nBeginning automatic 1D extraction for target %s in the %s configuration'%(id,grating))
    print('Press "h" or "?" for help.')

    global ax2d, ax1d, axex
    global ax2d, ax1d, axex, origexlims, extracted, extracted1d
    global ix, iy
    global spat, wspat, spaterr
    global coords, xex1, xex2, xex1temp, vspans
    global ax2d, ax1d, axex, orig2dlims, orig1dlims, origexlims
    global spat, wspat, spaterr
    global lastpressed
    global spec1d, var1d, err1d, spec1dbox, var1dbox, err1dbox
    global fitcomplete, fxlow, fxhigh, lfit, lerr
    global mouseclicks, mousereleases, keypresses
    global zcatarr, lineplotarr
    global vm, im2d
    global restwlax1d, rax1d, wlax2dflag, rax2d
    global yexp, yexpex, plotyexp
    global zax, zxlow, zxhigh, zylow, zyhigh, zoomcomplete
    global saved
    global mex1, mex2, mspans, mex1temp, lastclicked
    global zmean, boxcar, bc

    zaexcat = read(path1d + '00_redshift_catalog.txt')
    zaexcat = zaexcat[zaexcat['id'].data == id]
    zmean = np.average(zaexcat['z'].data, weights=zaexcat['zerr'].data**-2.)
    print('zmean = %.5f'%zmean)

    spec *= shdr['TOMUJY']
    err *= shdr['TOMUJY']

    yposerr = -1.
    yfwhmerr = -1.
    exsnr = -1.

    nx = shdr['NAXIS1']
    ny = shdr['NAXIS2']
    #ypos = shdr['SRCYPOS']
    #yguess = 20.
    #if ypos < 0.: ypos = 1.+ypos
    #yguess = float(ny) - ypos*float(ny)
    #yguess = float(ny)/2.*(1 + ypos)
    ytrace = shdr['YTRACE']

    wl = shdr['crval1'] + np.arange(shdr['naxis1'])*shdr['cdelt1']

    mex1 = []
    mex2 = []
    mspans = []
    mex1temp = None
    lastclicked = None

    efig = pl.figure(2, figsize=(16, 9))
    efig.suptitle('id=%s     grating=%s     automatic extraction'%(id, grating))

    # initialize axes, size and spacing figure setup
    lbuffer = 0.05
    rbuffer = 0.02
    bbuffer = 0.07
    tbuffer = 0.04
    wspace = 0.05
    hspace = 0.06

    h2d = 0.15

    w1d = 0.6

    plotyexp = True

    ax2d = efig.add_axes([lbuffer, 1.-tbuffer-h2d, 1.-rbuffer-lbuffer, h2d])
    yexp = ax2d.axhline(ytrace, color='orange', linestyle=':', linewidth=1.0)
    ax2drat = h2d/(1.-rbuffer-lbuffer)
    def reset2d():
        ax2d.clear()
        ax2d.set_title('2D science spectrum')
        if plotyexp: ax2d.add_line(yexp)
    reset2d()

    ax1d = efig.add_axes([lbuffer+(1.-rbuffer-w1d-wspace-lbuffer)+wspace, bbuffer, w1d, 1.-tbuffer-h2d-hspace-bbuffer])
    #ax1d.yaxis.get_ticklocs(minor=True)
    #ax1d.minorticks_on()
    #ax1d.tick_params(which='minor', length=3, direction='out')
    #ax1d.tick_params(which='major', length=6)
    def reset1d():
        global mex1, mex2, mspans
        ax1d.clear()
        ax1d.set_title('extracted 1D spectrum')
        ax1d.set_xlabel(r'$\lambda_{\mathrm{obs}}$ ($\mu$m)', fontsize=15)
        ax1d.set_ylabel(r'F$_\nu$ ($\mu$Jy)', fontsize=15)
        mex1 = []
        mex2 = []
        if len(mspans) > 0:
            [mm.remove() for mm in mspans]
            pl.draw()
            mspans = []
    reset1d()

    axex = efig.add_axes([lbuffer, bbuffer, 1.-rbuffer-w1d-wspace-lbuffer, 1.-tbuffer-h2d-hspace-bbuffer])
    yexpex = axex.axvline(ytrace, color='orange', linestyle=':', linewidth=1.0)
    def resetex():
        axex.clear()
        axex.set_title('extraction profile')
        axex.set_xlabel('y pixels', fontsize=15)
        axex.set_ylabel('signal collapsed over x ranges', fontsize=15)
        if plotyexp: axex.add_line(yexpex)
    resetex()

    # 2D spectrum
    vm = np.sort(spec.flatten()[spec.flatten() != 0.])[int(0.97*len(spec.flatten()[spec.flatten() != 0.]))]
    origvm = vm
    im2d = ax2d.imshow(spec, vmin=-vm, vmax=vm, cmap='viridis', origin='lower')

    orig2dlims = ax2d.axis()
    orig1dlims = ax1d.axis()
    origexlims = axex.axis()

# pl.rcParams['keymap.pan'].remove('p')
# pl.rcParams['keymap.save'].remove('s')
# pl.rcParams['keymap.fullscreen'].remove('f')
# pl.rcParams['keymap.home'].remove('r')
# pl.rcParams['keymap.zoom'].remove('o')

    def gauss(x,u, sig, flux): return flux*np.exp(-((u-x)**2./(2.*sig**2.)))/np.sqrt(2*np.pi*sig**2.)
    # line with a constant continuum
    def gaussline(x, u, sig, flux, a): return gauss(x, u, sig, flux) + a

    def uJy_to_cgs(wl): return (1.e-23)*(3.e14)/(wl**2.)*1.e-6   # wl in um

    pfit = [ypos, yfwhm/2.355, 1.]
    extracted1d = True
    print('Automatic profile parameters: y center=%.2f\ty sigma=%.2f'%(ypos,yfwhm/2.355))
    plrange = np.linspace(np.amin(np.arange(ny)), np.amax(np.arange(ny)), 200)
    axex.plot(plrange, gauss(plrange, *pfit), 'r-')
    #axex.axvline(yguess, color='orange', linestyle='--')
    axex.axvline(pfit[0], color='r', linestyle='--')

    exyinds = np.array([rr for rr in range(int(ymin), int(ymax)+1)])
    exymin = np.amin(exyinds)
    exymax = np.amax(exyinds)

    axex.plot(plrange[(plrange >= min(exyinds)) & (plrange <= max(exyinds))], gauss(plrange[(plrange >= min(exyinds)) & (plrange <= max(exyinds))], pfit[0], pfit[1], pfit[2]), 'c-')

    axex.text(0.05, 0.85, 'fit y center = %.2f'%pfit[0], transform=axex.transAxes)
    axex.text(0.05, 0.8, 'fit y fwhm = %.2f'%(pfit[1]*2.355), transform=axex.transAxes)

    spec1d = []
    var1d = []
    spec1dbox = []
    var1dbox = []
    for xx in range(nx):
        spec1d += [np.sum( gauss(exyinds, pfit[0], pfit[1], 1.)*spec[exyinds,xx]/(err[exyinds,xx]**2.) )/np.sum( (gauss(exyinds, pfit[0], pfit[1], 1.)**2.)/(err[exyinds,xx]**2.) )]
        var1d += [np.sum( gauss(exyinds, pfit[0], pfit[1], 1.) )/np.sum( (gauss(exyinds, pfit[0], pfit[1], 1.)**2.)/(err[exyinds,xx]**2.) )]

        spec1dbox += [np.sum(spec[exyinds,xx])]
        var1dbox += [np.sum(err[exyinds,xx]**2.)]

    spec1d = np.array(spec1d)
    var1d = np.array(var1d)
    err1d = np.sqrt(var1d)
    spec1dbox = np.array(spec1dbox)
    var1dbox = np.array(var1dbox)
    err1dbox = np.sqrt(var1dbox)

    spec1dbox[np.isnan(spec1d)] = np.nan

    ax1d.plot(wl, spec1d, 'k-', drawstyle='steps-mid', label='optimal')

    bc = None
    boxcar = False
    bc, = ax1d.plot(wl, spec1dbox, 'g-', drawstyle='steps-mid', label='boxcar', zorder=-1)
    bc.remove(); boxcar = False
    #ax1d.fill_between(wl, 0., err1dbox, color='g', alpha=0.3)

    ax1d.fill_between(wl, 0., err1d, color='0.7')
    #ax1d.legend(loc='upper right')
    orig1dlims = ax1d.axis()

    origexlims = axex.axis()
    #pl.draw()


# # check if redshift catalog exists, if not then create it, if so then open it
# zcatname = '00_redshift_catalog.txt'
# if os.path.exists(path1d+zcatname):
#  zcat = open(path1d+zcatname, 'a')
# else:
#  zcat = open(path1d+zcatname, 'w')
#  zcat.write('# id\tline\tlambda0\tlambdaobs\tlambdaobserr\tz\tzerr\tflux\tfluxerr\tgrating\n')

    # check if extraction info catalog exists, if not then create it, if so then open it
    excatname = '00_extract_info.txt'
    if os.path.exists(path1d+excatname):
        excat = open(path1d+excatname, 'a')
    else:
        excat = open(path1d+excatname, 'w')
        excat.write('# id\tgrating\typos\typoserr\tyfwhm\tyfwhmerr\tysnr\tymin\tymax\tautoflag\n')

    saved = False

    #### interactive functions ####
    def onclick(event):
        global ax2d, ax1d, axex
        global spec, err, shdr
        global mex1, mex2, mspans, mex1temp, lastclicked, wl
        global coords, xex1, xex2, xex1temp, extracted1d

        if event.inaxes == ax1d and event.button == MouseButton.LEFT and extracted1d:
            ix, iy = event.xdata, event.ydata
            #mex1.append(ix)
            mex1temp = ix

        lastclicked = event.inaxes

    def onrelease(event):
        global ax2d, ax1d, axex, origexlims, extracted1d
        global spat, wspat, spaterr
        global coords, xex1, xex2, xex1temp, vspans
        global mex1, mex2, mspans, lastclicked, wl

        if event.inaxes == ax1d and event.button == MouseButton.LEFT and lastclicked == ax1d and extracted1d:
            ix, iy = event.xdata, event.ydata
            mex1.append(mex1temp)
            mex2.append(ix)
            mspans += [ax1d.axvspan(mex1[-1], mex2[-1], color='b', edgecolor='none', alpha=0.3, zorder=1000)]
            pl.draw()

        lastclicked = None

    lastpressed = None

    zax = None
    zxlow = None
    zxhigh = None
    zylow = None
    zyhigh = None
    zoomcomplete = True

    fitcomplete = True
    lfit = None
    lerr = None
    fxlow = None
    fxhigh = None

    restwlax1d = False
    wlax2dflag = 0  # 0 if pixels, 1 if observed wl, 2 if rest wl

    plotyexp = True

# zcatarr = []
# def printzcat():
#  if len(zcatarr) > 0:
#   print('\nCurrent redshift catalog for this object:')
#   print(('line\t\tlambda0\t\tlambdaobs\terr\tzline\tzlineerr\tgrating'))
#   [print('%-13s\t%.3f\t%.3f\t%.3f\t%.5f\t%.5f\t%s'%(zz[0],zz[1],zz[2],zz[3],zz[4],zz[5],zz[8])) for zz in zcatarr]
#   print('Mean redshift = %.5f'%(np.average([zz[4] for zz in zcatarr], weights=[zz[5]**-2. for zz in zcatarr])))
#  else:
#   print('\nRedshift catalog is empty!')

# def zmean():
#  if len(zcatarr) > 0: return np.average([zz[4] for zz in zcatarr], weights=[zz[5]**-2. for zz in zcatarr])
#  else: return 0.
# lineplotarr = []

    # write extraction catalog entry
    def writeexcat():
        if len(excatarr) > 0:
            for ee in excatarr:
                excat.write('%s\t%s\t%.3f\t%.3f\t%.3f\t%.5f\t%.5f\t%.3e\t%.3e\n'%(id,*ee))

    def makemask1d(wl, mex1, mex2):
        mask1d = []
        for ww in wl:
            maskval = 0
            for mm in range(len(mex1)):
                if mex1[mm] < ww < mex2[mm]: maskval = 1
            mask1d += [maskval]
        mask1d = np.array(mask1d)
        return mask1d

    def writespec1d(wl, spec1d, err1d, spec1dbox, err1dbox, mex1, mex2):
        global zmean
        wlcgs = wl*1.e4
        cgscoeff = (1.e-23)*(3.e18)/(wlcgs**2.)*1.e-6
        spec1dcgs = spec1d*cgscoeff
        err1dcgs = err1d*cgscoeff
        spec1dboxcgs = spec1dbox*cgscoeff
        err1dboxcgs = err1dbox*cgscoeff

        mask1d = makemask1d(wl, mex1, mex2)
        mask1d[np.logical_not(np.isfinite(err1d))] = 1

        aoutname = path1d + outroot + '_1d.txt'
        foutname = path1d + outroot + '_1d.fits'

        # ascii output
        with open(aoutname, 'w') as outfile:
            outfile.write('lambda\tlambda_cgs\tflux\terr\tflux_cgs\terr_cgs\tflux_box\terr_box\tflux_box_cgs\terr_box_cgs\tmask\n')
            [outfile.write('%.7f\t%.2f\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.0f\n'%(wl[i], wlcgs[i], spec1d[i], err1d[i], spec1dcgs[i], err1dcgs[i], spec1dbox[i], err1dbox[i], spec1dboxcgs[i], err1dboxcgs[i], mask1d[i])) for i in range(len(wl))]

        # fits output file
        cwl = pf.Column(name='lambda', format='D', array=wl, unit='um')
        cwlcgs = pf.Column(name='lambda_cgs', format='D', array=wlcgs, unit='Angstrom')
        cflux = pf.Column(name='flux', format='D', array=spec1d, unit='uJy')
        cerr = pf.Column(name='err', format='D', array=err1d, unit='uJy')
        cfluxbox = pf.Column(name='flux_box', format='D', array=spec1dbox, unit='uJy')
        cerrbox = pf.Column(name='err_box', format='D', array=err1dbox, unit='uJy')
        cfluxcgs = pf.Column(name='flux_cgs', format='D', array=spec1dcgs, unit='erg/s/cm2/AA')
        cerrcgs = pf.Column(name='err_cgs', format='D', array=err1dcgs, unit='erg/s/cm2/AA')
        cfluxboxcgs = pf.Column(name='flux_box_cgs', format='D', array=spec1dboxcgs, unit='erg/s/cm2/AA')
        cerrboxcgs = pf.Column(name='err_box_cgs', format='D', array=err1dboxcgs, unit='erg/s/cm2/AA')
        cmask = pf.Column(name='mask', format='K', array=mask1d, unit='none')
        colarr = [cwl, cwlcgs, cflux, cerr, cfluxbox, cerrbox, cfluxcgs, cerrcgs, cfluxboxcgs, cerrboxcgs, cmask]
        coldefs = pf.ColDefs(colarr)
        outhdu = pf.BinTableHDU.from_columns(coldefs)
        outhdr = outhdu.header
        outhdr.set('ID', id, 'target ID', before='TTYPE1')
        program = outroot.split('-')[0]
        maskname = outroot.split('-')[1]
        outhdr.set('PROGRAM', program, 'program name', before='TTYPE1')
        outhdr.set('MASKNAME', maskname, 'maskname', before='TTYPE1')
        outhdr.set('YSCALE', shdr['CDELT2'], 'cross-dispersion scale in fraction of microshutter (0.53") per pixel')
        outhdr.set('YSCALEAS', shdr['CDELT2']*0.53, 'cross-dispersion in arcseconds per pixel')
        #outhdr.set('DETECTORS', detectors, 'detectors', before='TTYPE1')
        #[outhdr.set('FILE2D'+ff, , 'target ID', before='TTYPE1') for ff in ]

# key words to copy from 2D shdr
# 'NCOMBINE', 'EXPTIME', 'OTHRESH', 'TOMUJY', 'PROFCEN', 'PROFSIG', 'PROFSTRT', 'PROFSTOP', 'YTRACE', 'FILTER', 'GRATING'
        skeys = ['NCOMBINE', 'EXPTIME', 'OTHRESH', 'TOMUJY', 'PROFCEN', 'PROFSIG', 'PROFSTRT', 'PROFSTOP', 'YTRACE', 'FILTER', 'GRATING']
        for ss in skeys:
            outhdr.set(ss, shdr[ss], shdr.comments[ss], before='TTYPE1')

        outtime = gmtime()
        outhdr.set('EXDATE', '%02.f-%02.f-%04.f'%(outtime[2], outtime[1], outtime[0]), 'extraction date UTC', before='TTYPE1')
        outhdr.set('EXTIME', '%02.f:%02.f:%02.f'%(outtime[3], outtime[4], outtime[5]), 'extraction time UTC', before='TTYPE1')
        outhdr.set('EXYPOS', ypos, 'extracted y position', before='TTYPE1')
        outhdr.set('EXYPOSERR', yposerr, 'extracted y position error', before='TTYPE1')
        outhdr.set('EXYFWHM', yfwhm, 'extracted y fwhm', before='TTYPE1')
        outhdr.set('EXYFWHMERR', yfwhmerr, 'extracted y fwhm error', before='TTYPE1')
        outhdr.set('EXSNR', exsnr, 'extraction gaussian SNR', before='TTYPE1')
        outhdr.set('EXYMIN', exymin, 'minimum y of extraction window', before='TTYPE1')
        outhdr.set('EXYMAX', exymax, 'maximum y of extraction window', before='TTYPE1')
        outhdr.set('EXTYPE', 'auto', 'manual or auto extraction', before='TTYPE1')
        outhdr.set('EXWIDTH', float(exymax - exymin + 1)*shdr['CDELT2']*0.53, 'width of extraction window in arcseconds', before='TTYPE1')
        outhdr.set('EXYUP', float(exymax + 0.5) - ypos, 'pixels above trace for extraction window', before='TTYPE1')
        outhdr.set('EXYDOWN', ypos - float(exymin - 0.5), 'pixels below trace for extraction window', before='TTYPE1')
        outhdr.set('EXYUPAS', (float(exymax + 0.5) - ypos)*shdr['CDELT2']*0.53, 'arcseconds above trace for extraction window', before='TTYPE1')
        outhdr.set('EXYDOWNAS', (ypos - float(exymin - 0.5))*shdr['CDELT2']*0.53, 'arcseconds below trace for extraction window', before='TTYPE1')
        #if len(zcatarr) > 0: zm = np.average([zz[4] for zz in zcatarr], weights=[zz[5]**-2. for zz in zcatarr])
        #else: zm = -1.
        outhdr.set('ZINIT', zmean, 'mean redshift from lines', before='TTYPE1')

        #print(outhdu.header)
        outhdu.writeto(foutname, overwrite=True)
        print('Automatically extracted 1D spectrum saved to:')
        print('%s'%aoutname)
        print('%s'%foutname)
        print('\n')


    def onkeypress(event):
        global coords, xex1, xex2, vspans
        global ax2d, ax1d, axex, orig2dlims, orig1dlims, origexlims
        global spat, wspat, spaterr
        global lastpressed
        global spec1d, var1d, err1d, spec1dbox, var1dbox, err1dbox, extracted1d, saved
        global fitcomplete, fxlow, fxhigh, lfit, lerr
        global mouseclicks, mousereleases, keypresses
        global zcatarr, lineplotarr
        global vm, im2d
        global restwlax1d, rax1d, wlax2dflag, rax2d
        global yexp, yexpex, plotyexp
        #global yposerr, yfwhmerr, exsnr, exymin, exymax
        #global ypos, yposerr, yfwhm, yfwhmerr, exsnr, exymin, exymax
        global mex1, mex2, mspans
        global bc, boxcar
        global zmean

        #print('press', event.key)

        # remove last entry from 1D mask catalog
        if event.key == 'u':
            if len(mspans) > 0:
                mspans[-1].remove()
                mspans = mspans[:-1]
                mex1 = mex1[:-1]
                mex2 = mex2[:-1]
                pl.draw()

        # toggle boxcar extraction in 1D plot
        if event.key == 'b' and extracted1d:
            if boxcar: bc.remove(); pl.draw(); boxcar = False
            else: ax1d.add_artist(bc); pl.draw(); boxcar = True

#     # print the extraction windows, redshift catalog
        if event.key == 'p':
            print('\n')
#      print('extraction windows x low = %s'%xex1)
#      print('extraction windows x high = %s'%xex2)
            print('1D mask windows lambda_obs low (um) = %s'%[round(mm, 5) for mm in mex1])
            print('1D mask windows lambda_obs high (um) = %s'%[round(mm, 5) for mm in mex2])
#      printzcat()

        # toggle expected y position of main target on and off
        if event.key == 'y':
            if plotyexp: yexp.remove(); yexpex.remove(); pl.draw(); plotyexp = False
            else: ax2d.add_line(yexp); axex.add_line(yexpex); pl.draw(); plotyexp = True

#     # fit an emission line with a Gaussian profile
#     if event.key == 'f':
#      if event.inaxes != ax1d: fitcomplete = True
#      else:
#       if lastpressed != 'f' or (lastpressed == 'f' and fitcomplete == True):  # set lower bounds
#        fxlow = event.xdata
#        fitcomplete = False
#       elif lastpressed == 'f' and fitcomplete == False:  # set higher bounds and complete fit
#        fxhigh = event.xdata
#        fxlow, fxhigh = np.sort([fxlow, fxhigh])
#        ax1d.set_xlim(np.sort([fxlow, fxhigh]))
#        print('\nFitting a line between %.0f AA and %.0f AA'%(fxlow*1.e4,fxhigh*1.e4))
#
#        fitinds = np.where( (wl >= fxlow) & (wl <= fxhigh) & np.isfinite(spec1d) & (np.isfinite(err1d)) )[0]
#
#        contguess = np.median(spec1d[fitinds])
#        fluxguess = np.sum(spec1d[fitinds] - contguess)
#        p0 = [wl[fitinds][np.argmax(spec1d[fitinds])], 50./3.e5*np.mean(wl[fitinds]), fluxguess, contguess]
#        lfit, lcov = curve_fit(gaussline, wl[fitinds], spec1d[fitinds], sigma=err1d[fitinds], p0=p0)
#        lerr = np.sqrt(np.diag(lcov))
#        #print(lfit)
#        print('Fit parameters:\tcent=%.1f AA\tsigma=%.2f AA\tSNR=%.2f'%(lfit[0]*1.e4,lfit[1]*1.e4,lfit[2]/lerr[2]))
#        frange = np.linspace(wl[fitinds[0]], wl[fitinds[-1]], 500)
#        fc, = ax1d.plot(frange, [lfit[-1] for ww in frange], 'b--')
#        fg, = ax1d.plot(frange, gaussline(frange, *lfit), 'c-')
#        fu = ax1d.axvline(lfit[0], color='c', linestyle=':')
#        fparr = [fc, fg, fu]
#
#        fymin = np.amin(gaussline(frange, *lfit))
#        fymax = np.amax(gaussline(frange, *lfit))
#        ax1d.set_ylim(np.sort([fymin - 0.40*(fymax - fymin), fymax + 0.20*(fymax - fymin)]))
#        pl.draw()
#        fitcomplete = True
#
#        # oii, neiii, hd, hg, hb, oiii4959, oiii5007, ha, nii6585, sii6718, sii6731, siii9069, siii9531
#        linenamearr = ['[OII]doublet', '[NeIII]3869', 'Hdelta', 'Hgamma', 'Hbeta', '[OIII]4959', '[OIII]5007', '[NII]6548', 'Halpha', '[NII]6584', '[SII]6716', '[SII]6731', '[SIII]9069', '[SIII]9531']
#        linecent0arr = [round((3727.092+3729.875)/2.,3), 3869.856, 4102.89, 4341.68, 4862.68, 4960.295, 5008.240, 6549.86, 6564.61, 6585.27, 6718.29, 6732.67, 9071.09, 9533.21]
#        linenumarr = ['%s'%ll for ll in range(len(linenamearr))]
#        linestrarr = ['%.0f:\t%-13s\tcent0 = %s AA\tz = %.5f\n'%(ll, linenamearr[ll], linecent0arr[ll], lfit[0]*1.e4/linecent0arr[ll]-1.) for ll in range(len(linenamearr))]
#        linestr = ''.join(linestrarr)
#
#      # UV lines
#      #linecents = [1548.19+1., 1550.77+1., 1640.42, 1660.81, 1666.15, 1749.67, 1752.16, 1882.47, 1892.03, 1906.68, 1908.73, 2321.66, 2325.40, 2326.93, 2328.12, 2471.03]
#      #linenames = ['civ1548', 'civ1550', 'heii1640', 'oiii1661', 'oiii1666', 'niii1750', 'niii1752', 'siiii1882', 'siiii1892', 'ciii1907', 'ciii1909', 'oiii2322', 'cii2325', 'cii2327', 'cii2328', 'oii2471']
#      # optical lines
#      #linecents += [(3727.092+3729.875)/2., 3869.856, (3890.166+3889.750)/2., (3971.198+3968.593)/2., 4077.500, 4102.89, 4341.68, 4364.436, 4687.021, 4862.68, 4960.295, 5008.240, 5756.24, 5877.25, 6302.046, 6313.806, 6365.535, 6549.86, 6564.61, 6585.27, 6679.995, 6718.29, 6732.67, 7067.138, 7137.767, 7321.937, 7332.210, 9071.09, 9533.21]
#      #linenames += ['oii3727', 'neiii3869', 'hihei3888', 'hineiii3970', 'sii4076', 'hd', 'hg', 'oiii4363', 'heii4686', 'hb', 'oiii4959', 'oiii5007', 'nii5755', 'hei5876', 'oi6300', 'siii6312', 'oi6363', 'nii6548', 'ha', 'nii6584', 'hei6680', 'sii6716', 'sii6731', 'hei7067', 'ariii7135', 'oii7320', 'oii7330', 'siii9069', 'siii9531']
#
#        # disconnect interactive functions
#        efig.canvas.mpl_disconnect(mouseclicks)
#        efig.canvas.mpl_disconnect(mousereleases)
#        efig.canvas.mpl_disconnect(keypresses)
#
#        num = input('\nIdentify the line:\n\n' + linestr + 'c: Custom line name/rest wavelength\nEnter selection: ')
#
#        if num == 'c':  # enter custom line name/rest wavelength
#         linename = input('Enter line name (e.g., CIII]1908): ')
#         linecent0 = float(input('Enter rest-frame wavelength in Angstroms: '))
#         linecent = lfit[0]*1.e4
#         linez = linecent/linecent0 - 1.
#         linecenterr = lerr[0]*1.e4; linezerr = linecenterr/linecent0
#         lineflux = lfit[2]*uJy_to_cgs(lfit[0]); linefluxerr = lerr[2]*uJy_to_cgs(lfit[0])
#         print('Line identified as %s: z = %.5f'%(linename, linez))
#         print('Saving line %s to redshift catalog.'%linename)
#         zcatarr += [[linename, linecent0, linecent, linecenterr, linez, linezerr, lineflux, linefluxerr, grating]]
#         printzcat()
#         lineplotarr += [fparr]
#        elif num not in linenumarr:
#         print('Invalid selection.  Exiting without saving line fit.')
#         [ff.remove() for ff in fparr]
#         pl.draw()
#        else:
#         lineplotarr += [fparr]
#         num = int(num)
#         linecent = lfit[0]*1.e4; linecent0 = linecent0arr[num]; linename = linenamearr[num]
#         linez = linecent/linecent0 - 1.
#         linecenterr = lerr[0]*1.e4; linezerr = linecenterr/linecent0
#         lineflux = lfit[2]*uJy_to_cgs(lfit[0]); linefluxerr = lerr[2]*uJy_to_cgs(lfit[0])
#         print('Line identified as %s: z = %.5f'%(linename, linez))
#         print('Saving line %s to redshift catalog.'%linename)
#         zcatarr += [[linename, linecent0, linecent, linecenterr, linez, linezerr, lineflux, linefluxerr, grating]]
#         printzcat()
#
#        # reconnect interactive functions
#        mouseclicks = efig.canvas.mpl_connect('button_press_event', onclick)
#        mousereleases = efig.canvas.mpl_connect('button_release_event', onrelease)
#        keypresses = efig.canvas.mpl_connect('key_press_event', onkeypress)

#     # modify redshift catalog
#     if event.key == 'm':
#      # disconnect interactive functions
#      efig.canvas.mpl_disconnect(mouseclicks)
#      efig.canvas.mpl_disconnect(mousereleases)
#      efig.canvas.mpl_disconnect(keypresses)
#
#      print('\nCurrent redshift catalog for this object:')
#      print(('\tline\t\tlambda0\t\tlambdaobs\terr\tzline\tzlineerr'))
#      [print('%s\t%-13s\t%.3f\t%.3f\t%.3f\t%.5f\t%.5f'%(zz,zcatarr[zz][0],zcatarr[zz][1],zcatarr[zz][2],zcatarr[zz][3],zcatarr[zz][4],zcatarr[zz][5])) for zz in range(len(zcatarr))]
#      num = input('Select a number to remove from the redshift catalog: ')
#      if num not in ['%s'%rr for rr in range(len(zcatarr))]:
#       print('Invalid selection.  Exiting without modifying the redshift catalog.')
#      else:
#       num = int(num)
#       zcatarr = zcatarr[:num] + zcatarr[num+1:]
#       [ff.remove() for ff in lineplotarr[num]]
#       lineplotarr = lineplotarr[:num] + lineplotarr[num+1:]
#       pl.draw()
#       printzcat()
#
#      # reconnect interactive functions
#      mouseclicks = efig.canvas.mpl_connect('button_press_event', onclick)
#      mousereleases = efig.canvas.mpl_connect('button_release_event', onrelease)
#      keypresses = efig.canvas.mpl_connect('key_press_event', onkeypress)

        # toggle rest wavelength in ax1d based on z_mean from current redshift catalog
        if event.key == 'w' and event.inaxes == ax1d:
            if not restwlax1d:
                zm = zmean
                rax1d = ax1d.secondary_xaxis('bottom', functions=(lambda wl: wl*1.e4/(1.+zm), lambda wl: wl/1.e4*(1.+zm)))
                rax1d.set_xlabel(r'$\lambda_{\mathrm{rest}}$ ($\AA$)', fontsize=12, labelpad=-32)
                rax1d.tick_params(axis='x', direction='in', pad=-15)
                restwlax1d = True
                pl.draw()
            elif restwlax1d:
                rax1d.remove()
                restwlax1d = False
                pl.draw()

        # toggle from pixels to observed wl to rest wl in ax2d
        if event.key == 'w' and event.inaxes == ax2d:
            if wlax2dflag == 0:  # convert to obs
                rax2d = ax2d.secondary_xaxis('bottom', functions=(lambda pix: shdr['crval1']+pix*shdr['cdelt1'], lambda wl: (wl-shdr['crval1'])/shdr['cdelt1']))
                rax2d.set_xlabel(r'$\lambda_{\mathrm{obs}}$ ($\mu$m)', fontsize=12, labelpad=0)
                ax2d.tick_params(bottom=False, labelbottom=False)
                wlax2dflag = 1
                pl.draw()
            elif wlax2dflag == 1:  # convert from obs to rest based on zmean
                rax2d.remove()
                zm = zmean
                rax2d = ax2d.secondary_xaxis('bottom', functions=(lambda pix: (shdr['crval1']+pix*shdr['cdelt1'])*1.e4/(1.+zm), lambda wl: (wl*(1.+zm)/1.e4-shdr['crval1'])/shdr['cdelt1']))
                rax2d.set_xlabel(r'$\lambda_{\mathrm{rest}}$ ($\AA$)', fontsize=12, labelpad=0)
                ax2d.tick_params(bottom=False, labelbottom=False)
                wlax2dflag = 2
                pl.draw()
            elif wlax2dflag == 2:  # convert from rest to pix
                rax2d.remove()
                wlax2dflag = 0
                ax2d.tick_params(bottom=True, labelbottom=True)
                pl.draw()

        # zoom in an axis
        if event.key == 'z':
            global zax, zxlow, zxhigh, zylow, zyhigh, zoomcomplete
            if zax != event.inaxes: zoomcomplete = True
            if lastpressed != 'z' or (lastpressed == 'z' and zoomcomplete == True):  # set lower bounds
                zax = event.inaxes
                zxlow, zylow = event.xdata, event.ydata
                zoomcomplete = False
            elif lastpressed == 'z' and zoomcomplete == False:  # set higher bounds
                if event.inaxes == zax:
                    zxhigh, zyhigh = event.xdata, event.ydata
                    if zax != ax2d:
                        zax.set_xlim(np.sort([zxlow, zxhigh]))
                        zax.set_ylim(np.sort([zylow, zyhigh]))  # never change the 2D y range
                    else:
                        zxcent = (zxhigh+zxlow)/2.
                        figsize = efig.get_size_inches()
                        zxrange = float(shdr['naxis2'])/ax2drat*(figsize[0]/figsize[1])
                        event.inaxes.set_xlim([zxcent-zxrange/2., zxcent+zxrange/2.])
                    pl.draw()
                zoomcomplete = True

        # reset axis zoom
        if event.key == 'x':
            if event.inaxes == ax2d: ax2d.axis(orig2dlims); pl.draw()
            if event.inaxes == ax1d: ax1d.axis(orig1dlims); pl.draw()
            if event.inaxes == axex: axex.axis(origexlims); pl.draw()

        if event.key == 'h' or event.key == '?':
            helpstr = '\n In the 2D science spectrum panel, click and drag to add pixels in that x-range to the spatial profile.\n\
       In the extracted 1D spectrum panel, click and drag to add wavelength ranges to the 1D mask.\n\
       Key bindings:\n\
      # p -- "print": Print useful information including the extraction windows and current redshift catalog.\n\
       y -- "yposition": Toggle marker for the expected y-position of the primary target (orange dotted line).\n\
       b -- "boxcar": Toggle a plot of the boxcar extraction on and off.\n\
      # f -- "fit": Fit an emission line with a Gaussian profile.  Press "f" twice, once on each side of the line, to fit the line.  The centroid guess is the highest data point in the selected wavelength range.  Follow the prompt to estimate the redshift and add the line to the redshift catalog.  Line fits currently stored in the redshift catalog will remain plotted in the 1D spectrum panel.\n\
      # m -- "modify": Modify the redshift catalog.  This function allows you to remove an entry from the redshift catalog.\n\
       u -- "undo mask": Remove the last masked region from the 1D mask catalog.\n\
       w -- "wavelength": Press "w" in the 2D spectrum panel to switch between pixel and observed wavelength on the x-axis.  If a redshift has been estimated from emission lines, "w" will switch between pixel, observed, and rest-frame wavelength in the 2D spectrum panel.  In the 1D spectrum panel, "w" will toggle rest-frame wavelength tick marks on the x-axis.\n\
       z -- "zoom": Press "z" twice in the same panel to zoom in to a range with corners defined by the cursor positions when "z" was pressed.\n\
       x -- "axis reset": Resets the plotted range to the original range for the axis in which the cursor is.\n\
       c -- "center": Move the point under the cursor to the center of the panel.\n\
       i -- "in": Zoom the selected axis in, keeping the same center coordinate.\n\
       o -- "out": Zoom the selected axis out, keeping the same center coordinate.\n\
       s -- "save": Save the extracted 1D spectrum, extraction information, and redshift catalog.\n\
       q -- "quit": Exit the 1D extraction.  Make sure to save first!\n\
       arrow keys: Use the arrow keys to shift the plot.  The extraction profile and 1D spectrum panels can be shifted up/down/left/right.  The 2D spectrum panel can only be shifted left and right, while up and down changes the stretch.\n'

            print(helpstr)


        # pan center of axis
        if event.key == 'c':
            newcx, newcy = event.xdata, event.ydata
            oldxlow, oldxhigh = event.inaxes.get_xlim()
            oldylow, oldyhigh = event.inaxes.get_ylim()
            oldxrange = (oldxhigh - oldxlow)
            oldyrange = (oldyhigh - oldylow)
            event.inaxes.set_xlim([newcx - oldxrange/2., newcx + oldxrange/2.])
            if event.inaxes != ax2d: event.inaxes.set_ylim([newcy - oldyrange/2., newcy + oldyrange/2.])
            pl.draw()

        # pan axis with arrow keys
        if event.key == 'left':
            oldxlow, oldxhigh = event.inaxes.get_xlim()
            oldylow, oldyhigh = event.inaxes.get_ylim()
            oldxrange = (oldxhigh - oldxlow)
            oldyrange = (oldyhigh - oldylow)
            event.inaxes.axis([oldxlow-0.1*oldxrange, oldxhigh-0.1*oldxrange, oldylow, oldyhigh])
            pl.draw()
        if event.key == 'right':
            oldxlow, oldxhigh = event.inaxes.get_xlim()
            oldylow, oldyhigh = event.inaxes.get_ylim()
            oldxrange = (oldxhigh - oldxlow)
            oldyrange = (oldyhigh - oldylow)
            event.inaxes.axis([oldxlow+0.1*oldxrange, oldxhigh+0.1*oldxrange, oldylow, oldyhigh])
            pl.draw()
        if event.key == 'down' and event.inaxes != ax2d:
            oldxlow, oldxhigh = event.inaxes.get_xlim()
            oldylow, oldyhigh = event.inaxes.get_ylim()
            oldxrange = (oldxhigh - oldxlow)
            oldyrange = (oldyhigh - oldylow)
            event.inaxes.axis([oldxlow, oldxhigh, oldylow-0.1*oldyrange, oldyhigh-0.1*oldyrange])
            pl.draw()
        if event.key == 'up' and event.inaxes != ax2d:
            oldxlow, oldxhigh = event.inaxes.get_xlim()
            oldylow, oldyhigh = event.inaxes.get_ylim()
            oldxrange = (oldxhigh - oldxlow)
            oldyrange = (oldyhigh - oldylow)
            event.inaxes.axis([oldxlow, oldxhigh, oldylow+0.1*oldyrange, oldyhigh+0.1*oldyrange])
            pl.draw()

        # increase or decrease stretch scale of 2D spectrum
        if event.key == 'up' and event.inaxes == ax2d:
            vm *= 2.
            im2d.set_clim(vmin=-vm, vmax=vm)
            pl.draw()
        if event.key == 'down' and event.inaxes == ax2d:
            vm /= 2.
            im2d.set_clim(vmin=-vm, vmax=vm)
            pl.draw()

        # zoom in and out with fixed center
        if event.key == 'i':
            oldxlow, oldxhigh = event.inaxes.get_xlim()
            oldylow, oldyhigh = event.inaxes.get_ylim()
            oldxrange = (oldxhigh - oldxlow)
            oldyrange = (oldyhigh - oldylow)
            oldxcent = (oldxlow + oldxhigh)/2.
            oldycent = (oldylow + oldyhigh)/2.
            if event.inaxes != ax2d: event.inaxes.axis([oldxcent-0.5*oldxrange/2.,oldxcent+0.5*oldxrange/2.,oldycent-0.5*oldyrange/2.,oldycent+0.5*oldyrange/2.])
            else:
                figsize = efig.get_size_inches()
                minxrange = float(shdr['naxis2'])/ax2drat*(figsize[0]/figsize[1])
                newxrange = max([0.5*oldxrange,minxrange])
                event.inaxes.axis([oldxcent-newxrange/2.,oldxcent+newxrange/2.,oldylow,oldyhigh])
            pl.draw()

        # zoom in and out with fixed center
        if event.key == 'o':
            oldxlow, oldxhigh = event.inaxes.get_xlim()
            oldylow, oldyhigh = event.inaxes.get_ylim()
            oldxrange = (oldxhigh - oldxlow)
            oldyrange = (oldyhigh - oldylow)
            oldxcent = (oldxlow + oldxhigh)/2.
            oldycent = (oldylow + oldyhigh)/2.
            if event.inaxes != ax2d: event.inaxes.axis([oldxcent-2.0*oldxrange/2.,oldxcent+2.0*oldxrange/2.,oldycent-2.0*oldyrange/2.,oldycent+2.0*oldyrange/2.])
            else:
                figsize = efig.get_size_inches()
                minxrange = float(shdr['naxis2'])/ax2drat*(figsize[0]/figsize[1])
                newxrange = max([2.0*oldxrange,minxrange])
                event.inaxes.axis([oldxcent-newxrange/2.,oldxcent+newxrange/2.,oldylow,oldyhigh])
            pl.draw()

        # save extracted spectrum, and redshift and extract catalogs
        if event.key == 's' and extracted1d:
            print ('Saving output')
            #writezcat()
            writespec1d(wl, spec1d, err1d, spec1dbox, err1dbox, mex1, mex2)
            excat.write('%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.1f\t%.1f\t%.0f\n'%(id, grating, ypos, yposerr, yfwhm, yfwhmerr, exsnr, exymin, exymax, 1))  # write extraction cat
            saved = True

        lastpressed = event.key

    #### end interactive functions ####

    mouseclicks = efig.canvas.mpl_connect('button_press_event', onclick)
    mousereleases = efig.canvas.mpl_connect('button_release_event', onrelease)
    keypresses = efig.canvas.mpl_connect('key_press_event', onkeypress)

    pl.show()
    pl.pause(10000)

    # disconnect interactive functions
    efig.canvas.mpl_disconnect(mouseclicks)
    efig.canvas.mpl_disconnect(mousereleases)
    efig.canvas.mpl_disconnect(keypresses)

    pl.close(efig)

# zcat.close()
    excat.close()

    #pl.figure(1)
    #print('finished 1d extraction')

    return saved

### end auto extraction

def picktarget():
    filelist = glob(path2d + '*/*_2d.fits')
    progarr = []
    idarr = []
    masknamearr = []
    detarr = []
    gratarr = []
    for ll in filelist:
        splitname = ll.strip().split('/')[-1].split('_')
        splitname2 = splitname[2].split('-')
        progarr += [splitname[0]]
        idarr += [int(splitname2[0])]
        masknamearr += [splitname[1]]
        gratarr += [splitname2[1]+'-'+splitname2[2]]
        detarr += [splitname2[3]]
    progarr = np.array(progarr)
    idarr = np.array(idarr)
    masknamearr = np.array(masknamearr)
    gratarr = np.array(gratarr)
    detarr = np.array(detarr)

    idlist = np.unique(idarr)
    proglist = np.unique(progarr)
    gratlist = np.unique(gratarr)

    if len(proglist) > 1: print('WARNING: more than one program present!')

    ntarg = len(idlist)

    print('Program = %s'%proglist[0])
    print('Grating configurations = %s'%gratlist)
    print('Number of targets = %s'%ntarg)
    print('\nSelect a target to extract:')
    optstr = 'option\tid number\t\t'
    for gg in gratlist: optstr += '%s\t'%gg
    print(optstr)
    #idstrarr = ['%.0f:\t%.0f\n'%(i, idlist[i]) for i in range(len(idlist))]
    idstrarr = []
    for i in range(len(idlist)):
        istr = '%.0f:\t%-15s'%(i, idlist[i])

        # list 1D files that have been extracted
        list1d = glob(path1d + '*-%s-*-*_1d.fits'%(idlist[i]))
        exgrats = []
        for ll in list1d:
            for gg in gratlist:
                if gg in ll: exgrats += [gg]
        extractstatus = [True if gg in exgrats else False for gg in gratlist] # if gg in exgrats else False]
        for ee in extractstatus:
            if ee: istr += '\t\tyes'
            else: istr += '\t\tno'
        istr += '\n'
        idstrarr += [istr]

    validoptions = ['%s'%i for i in range(len(idlist))]
    idstr = ''.join(idstrarr) + 'q: Quit'
    print(idstr)
    num = input('Enter selection: ')
    if num not in validoptions:
        print('Exiting.')
    else:
        num = int(num)
        id = idlist[num]
        inspect2d(id)
    return num

selection = 0
while selection != 'q':
    selection = picktarget()
