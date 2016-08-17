import numpy as np
import pylab as pl
import subprocess as sp
import uuid
import pdb
import os
import string as ss
# get rid of annoying fits table warnings
import warnings as warn
with warn.catch_warnings():
    warn.simplefilter("ignore")
    import atpy as ap
# import some kind of pyfits
try:
    import astropy.io.fits as pf
except:
    try:
        import pyfits as pf
    except:
        raise ImportError("astropy.io.fits nor pyfits found.")


# setup slow (hard disk) and fast (RAM disk) methods, then choose one
molecfitinstall_slow = '/Volumes/Data3/rcwh/data2/software/mf/'
molecfitinstall_fast = '/Volumes/RAMDisk/molecfit1.1.1/'
molecfitinstall_fasthsim = '/mnt/ramdisk/molecfit1.1.1/'
molecfitinstall_laptop = '/Users/rcwh/software/molecfit1.1.1/'
molecfittmp_slow = '/tmp/'
molecfittmp_fast = '/Volumes/RAMDisk/tmp/'
molecfittmp_fasthsim = '/mnt/ramdisk/tmp/'
molecfittmp_laptop = '/tmp/'

molecfitinstall = molecfitinstall_slow # molecfitinstall_laptop # _fasthsim #
molecfittmp     =     molecfittmp_slow # molecfittmp_laptop     # _fasthsim #


def createParFile(basedir, objfile, listname, outfile, outdir='./', trans=True, parfile=None, \
                  wrange_include=None, wrange_exclude=None, prange_exclude=None, \
                  columns=["LAMBDA", "FLUX", None, None], \
                  default_error=0.01, wlgtomicron=1.0, vac_air='vac', \
                  plot_type="N", plot_range=False, ftol=1e-2, xtol=1e-2, \
                  list_molec=['H2O', 'CH4', 'O3'], fit_molec=[True, True, True], relcol=[1.0,1.0,1.0], \
                  flux_unit=0, fit_back=True, telback=0.1, fit_cont=True, cont_n=3, cont_const=1.0, \
                  fit_wlc=True, wlc_n=3, wlc_const=0.0, fitresbox=False, kernmode=False, \
                  fit_res_box=False, relres_box=1.0, fit_res_gauss=True, res_gauss=1.0, fit_res_lorentz=False, res_lorentz=0.5, \
                  kernfac=10.0, varkern=False, kernel_file=None, obsdate=None, obsdate_key="MJD-OBS", \
                  utc=None, utc_key="UTC", telalt=None, telalt_key="ESO TEL ALT", rhum=None, \
                  rhum_key="ESO TEL AMBI RHUM", pres=None, pres_key="ESO TEL AMBI PRES START", \
                  temp=None, temp_key="ESO TEL AMBI TEMP", m1temp=None, m1temp_key="ESO TEL TH M1 TEMP", \
                  geoelev=None, geoelev_key="ESO TEL GEOELEV", longitude=None, \
                  longitude_key="ESO TEL GEOLON", latitude=None, latitude_key="ESO TEL GEOLAT", \
                  slitw=0.4, slitw_key="ESO INS SLIT1 WID", pixsc=0.086, pixsc_key=None, \
                  ref_atm="equ.atm", \
                  gdas_dir=molecfitinstall+"data/profiles/grib", gdas_prof="auto", layers=1, emix=5.0, pwv=-1.0):
    """
    RH 27/12/15

    Script to create a Molecfit1.1.1 parameter file, given various input files

    columns: wave, flux, fluxerr, mask

    #vac_air:   vac/air
    #dateKey:   MJD-OBS / DATE_VAL
    #timeKey:   TM-START / TIME_VAL
    #telAltKey: ESO TEL ALT / TELALT_VAL
    #VARFWHM:    0 / 1 
    #FLUXLIM:    -1 or cut value
    #PLOT_TYPE:  W / X / N
    
    """

    # create output parfilename
    if parfile==None:
        parfile = outdir+objfile[objfile.rfind("/"):objfile.rfind(".fits")]+".molecfit.par"
    # create outdir if not present
    if (not os.path.isdir(outdir)):
        os.mkdir(outdir)
        
    # open file handler
    fp = open(parfile, "w")
    # start writing info
    fp.write("\n")
    # dir and files
    fp.write("\nbasedir: "+basedir) # base directory
    fp.write("\nfilename: "+objfile) # object spectrum
    if listname==None:
        fp.write("\nlistname: "+str(listname).lower())# subsequent files to be corrected in list mode
    else:
        fp.write("\nlistname: "+listname)# subsequent files to be corrected in list mode
    fp.write("\ntrans: "+str(int(trans)))# type of spectrum: transmission = True/1, emission=False/0

    if wrange_include==None:
        fp.write("\nwrange_include: "+str(wrange_include).lower()) # none
    else:
        fp.write("\nwrange_include: "+wrange_include) # range to INclude in fit (ASCII or FITS table)
    if wrange_exclude==None:
        fp.write("\nwrange_exclude: "+str(wrange_exclude).lower()) # none
    else:
        fp.write("\nwrange_exclude: "+wrange_exclude) # range to EXclude in fit (")
    if prange_exclude==None:
        fp.write("\nprange_exclude: "+str(prange_exclude).lower()) # none
    else:
        fp.write("\nprange_exclude: "+prange_exclude) # range to EXclude in fit (")

    # generate colnames
    ncolumns = len(columns)
    if ncolumns>4: raise ValueError('Columns must have no more than 4 elements')
    if ncolumns<4:
        nnone=4-ncolumns
        extra = ["NULL" for x in xrange(nnone)]
        columns.extend(extra)
    if len(columns)!=4: raise ValueError("This should never happen")
    cns=''
    for cn in columns:
        cns=cns+str(cn)+" "
    fp.write("\ncolumns: "+cns)
    # input structure
    fp.write("\ndefault_error: "+str(default_error))
    fp.write("\nwlgtomicron: "+str(wlgtomicron))
    fp.write("\nvac_air: "+vac_air)
    # results
    fp.write("\noutput_dir: "+outdir)             # output dir
    fp.write("\noutput_name: "+outfile)           # corrected object spectrum
    # plotting
    fp.write("\nplot_creation: "+plot_type)
    fp.write("\nplot_range: "+str(int(plot_range)))
    # fit precision
    # fitting
    fp.write("\nftol: "+str(ftol))
    fp.write("\nxtol: "+str(xtol))
    # molec columns: define species present, and which to fit
    fp.write("\nlist_molec: "+ss.join(list_molec))
    fp.write("\nfit_molec: "+ss.join(np.array(np.array(fit_molec,dtype=np.int),dtype=str).tolist()))
    fp.write("\nrelcol: "+ss.join(np.array(np.array(relcol,dtype=np.float),dtype=str).tolist()))
    # background and continuum
    fp.write("\nflux_unit: "+str(int(flux_unit))) # what flux units to use. see manual, 0 most likely
    fp.write("\nfit_back: "+str(int(fit_back))) # fit telescope background?
    fp.write("\ntelback: "+str(telback)) # initial guess for telescope background fit
    fp.write("\nfit_cont: "+str(int(fit_cont))) # fit continuum?
    fp.write("\ncont_n: "+str(cont_n)) # polynomial degree for fit
    fp.write("\ncont_const: "+str(cont_const)) # initial term for continuum fit
    # wavelength solution
    fp.write("\nfit_wlc: "+str(int(fit_wlc))) # fit polynomial for new wavelength solution?
    fp.write("\nwlc_n: "+str(wlc_n)) # order of wavelength fit
    fp.write("\nwlc_const: "+str(wlc_const)) # initial guess for const term
    # resolution fitting
    fp.write("\nfit_res_box: "+str(int(fit_res_box))) # fit resolution by boxcar?
    fp.write("\nrelres_box: "+str(relres_box)) # initial guess of boxcar relative to slitwidth (0-2)
    fp.write("\nkernmode: "+str(int(kernmode))) # Voigt profile approx instead of Gauss and Lorent
    fp.write("\nfit_res_gauss: "+str(int(fit_res_gauss))) # fit resolution with Gaussian?
    fp.write("\nres_gauss: "+str(res_gauss)) # initial FWHM of Gauss in pixels
    fp.write("\nfit_res_lorentz: "+str(int(fit_res_lorentz))) # fit resolution with Lorentian?
    fp.write("\nres_lorentz: "+str(res_lorentz))
    fp.write("\nkernfac: "+str(kernfac)) # size of kernal in FWHM (?!)
    fp.write("\nvarkern: "+str(int(varkern))) # vary kernel size - linear with wavelength
    fp.write("\nkernel_file: "+str(kernel_file).lower()) # instead of parametric kernel, give ascii profile here
    
    # Ambient params

    # make helper function - either use Keyword val or passed value.
    def oneOrOther(thing,key,tag, upper=False):
        if thing==None:
            # use keyword to get value
            fp.write("\n"+tag)
            fp.write("\n"+tag+"_key: "+key)
        elif((key==None) & (thing!=None)):
            # force use of value, not keyword
            fp.write("\n"+tag+": "+str(thing))
            if (not upper):
                fp.write("\n"+tag+"_key: "+str(key).lower())
            else:
                fp.write("\n"+tag+"_key: "+str(key).upper())
        else:
            raise IOError("Inputs not understood")

    oneOrOther(obsdate,obsdate_key,"obsdate")
    oneOrOther(utc, utc_key, "utc")
    oneOrOther(telalt, telalt_key, "telalt")
    oneOrOther(rhum, rhum_key, "rhum")
    oneOrOther(pres, pres_key, "pres")
    oneOrOther(temp, temp_key, "temp")
    oneOrOther(m1temp, m1temp_key, "m1temp")
    oneOrOther(geoelev, geoelev_key, "geoelev")
    oneOrOther(longitude, longitude_key, "longitude")
    oneOrOther(latitude, latitude_key, "latitude")

    # Instrument params
    oneOrOther(slitw, slitw_key, "slitw", upper=True)
    oneOrOther(pixsc, pixsc_key, "pixsc", upper=True)

    # Atmospheric profiles
    fp.write("\nref_atm: "+ref_atm)
    fp.write("\ngdas_dir: "+gdas_dir)
    fp.write("\ngdas_prof: "+gdas_prof)
    fp.write("\nlayers: "+str(int(layers)))
    fp.write("\nemix: "+str(emix))
    fp.write("\npwv: "+str(pwv))
    
    fp.write("\nend\n")
    fp.close()

    return parfile

def runMolecfit(parfile, installdir=molecfitinstall, bin="bin/molecfit", options="m", quiet=True):
    """
    RH 4/1/16

    Wrapper for running the molecfit exec from the command line

    """

    if quiet:
        FNULL = open(os.devnull, 'w')
        sp.call([installdir+bin, parfile, options], stdout=FNULL, stderr=FNULL)
        FNULL.close()
    else:
        sp.call([installdir+bin, parfile, options])

def runCalctrans(parfile, installdir=molecfitinstall, bin="bin/calctrans", quiet=True):
    """
    RH 4/1/16

    Wrapper for running the molecfit exec from the command line

    """

    if quiet:
        FNULL = open(os.devnull, 'w')
        sp.call([installdir+bin, parfile], stdout=FNULL, stderr=FNULL)
        FNULL.close()
    else:
        sp.call([installdir+bin, parfile])


def fitstabMolecfitWrapper(wave, flux, \
                           obsdate, utc, telalt, rhum, pres, temp, m1temp, geoelev, longitude, latitude, slitw, pixsc, \
                           trans=True, default_error=1e-2, wlgtomicron=1.0, vac_air="vac", wrange_include=None, plot_range=True, \
                           ftol=1e-2, xtol=1e-2, list_molec=['H2O', 'O2', 'CH4', 'CO2'], fit_molec=[True, True, True, True], \
                           flux_unit=0, fit_back=False, fit_cont=1, cont_n=6, fit_wlc=True, wlc_n=6, wlc_const=0.0, \
                           fit_res_gauss=True, res_gauss=2.5, kernfac=10.0, varkern=False, \
                           outname="molecfit_fitsTab", tmp=molecfittmp, plot_type="N"):
    """
    RH 4/1/16

    Given wave and flux numpy arrays, create FITS files and run molecfit on them

    """

    dirname = tmp+"molecfit"+str(uuid.uuid4())
    objfile = dirname+"/OBJ"+str(uuid.uuid4())+".txt"

    sp.call(["mkdir", dirname])
    objtab = ap.Table()
    objtab.add_column('LAMBDA', wave)
    objtab.add_column('FLUX', flux)
    objtab.table_name="Obj"
    objtab.write(objfile, type="fits", verbose=False)

    # write wrange include file
    if wrange_include!=None:
        wrfile = dirname+"/wranges.txt"
        fp=open(wrfile, "w")
        for wr in wrange_include:
            fp.write(str(wr[0])+" "+str(wr[1])+"\n")
        fp.write("\n")
        fp.close()
        wrange_include = wrfile
            
    columns = ["LAMBDA", "FLUX"]
    parfile = createParFile(dirname, objfile, None, outname, outdir=dirname, \
                            columns=columns, \
                            obsdate=obsdate, obsdate_key=None, \
                            utc=utc, utc_key=None, \
                            telalt=telalt, telalt_key=None, \
                            rhum=rhum, rhum_key=None, \
                            pres=pres, pres_key=None, \
                            temp=temp, temp_key=None, \
                            m1temp=m1temp, m1temp_key=None, \
                            geoelev=geoelev, geoelev_key=None, \
                            longitude=longitude, longitude_key=None, \
                            latitude=latitude, latitude_key=None,\
                            slitw=slitw, slitw_key=None, \
                            pixsc=pixsc, pixsc_key=None, \
                            trans=trans, default_error=default_error, wlgtomicron=wlgtomicron, \
                            vac_air=vac_air, wrange_include=wrange_include, plot_range=plot_range, \
                            ftol=ftol, xtol=xtol, list_molec=list_molec, fit_molec=fit_molec, flux_unit=flux_unit, \
                            fit_back=fit_back, fit_cont=fit_cont, cont_n=cont_n, fit_wlc=fit_wlc, wlc_n=wlc_n, wlc_const=wlc_const, \
                            fit_res_gauss=fit_res_gauss, res_gauss=res_gauss, kernfac=10.0, varkern=False, plot_type=plot_type)
    
    runMolecfit(parfile)
    runCalctrans(parfile)
    
    
def asciiMolecfitWrapper(wave, flux, \
                         obsdate, utc, telalt, rhum, pres, temp, m1temp, geoelev, longitude, latitude, slitw, pixsc, \
                         err=None, mask=None, trans=True, default_error=1e-2, wlgtomicron=1.0, vac_air="vac", \
                         wrange_include=None, wrange_exclude=None, prange_exclude=None, plot_range=True, \
                         ftol=1e-2, xtol=1e-2, list_molec=['H2O', 'O2', 'CH4', 'CO2'], fit_molec=[True, True, True, True], \
                         relcol=[1.0,1.0,1.0,1.0],flux_unit=0, fit_back=False, fit_cont=1, cont_n=6, \
                         cont_const=0.0, fit_wlc=True, wlc_n=6, wlc_const=0.0, \
                         fit_res_gauss=True, res_gauss=2.5, kernfac=10.0, varkern=False, relres_box=0.0, \
                         outname="molecfit_results", tmp=molecfittmp,\
                         plot_type="N", calcNewWave=True, cleanTmp=True, basedir=molecfitinstall, doPlot=True, quiet=False):
    """
    RH 19/11/15

    Given wave and flux numpy arrays, create an ascii files and run molecfit on them
    
    """

    dirname = tmp+"molecfit"+str(uuid.uuid4())
    objfile = dirname+"/OBJ"+str(uuid.uuid4())+".txt"
    molecfitfile = dirname+"/"+outname+"_fit.fits"
    calctransfile = dirname+"/"+outname+"_tac.fits"###objfile[:objfile.rfind(".txt")]+"_TAC.fits"
    if calcNewWave: resFile = outfile+"_fit.res"

    # do some simple checks...
    # no nans
    good = np.where(np.isfinite(wave) & np.isfinite(flux) )
    wave = wave[good]
    flux = flux[good]
    if type(err)!=type(None): err = err[good]
    if type(mask)!=type(None): mask = mask[good]
        
    # wave is sorted (monotonically increasing?)
    sw = wave.argsort()
    wave = wave[sw]
    flux = flux[sw]
    if type(err)!=type(None): err = err[sw]
    if type(mask)!=type(None): mask = mask[sw]

    # now make the ascii files
    sp.call(["mkdir", dirname])
    if (type(mask)!=type(None)) & (type(err)==type(None)):
        np.savetxt(objfile, np.array([wave,flux,mask]).T)
        columns = ["LAMBDA", "FLUX", "NULL", "QUAL"]
    elif (type(mask)==type(None)) & (type(err)!=type(None)):
        np.savetxt(objfile, np.array([wave,flux,err]).T)
        columns = ["LAMBDA", "FLUX", "ERR", "NULL"]
    elif (type(mask)!=type(None)) & (type(err)!=type(None)):
        np.savetxt(objfile, np.array([wave,flux,err,mask]).T)
        columns = ["LAMBDA", "FLUX", "ERR", "QUAL"]
    else:
        np.savetxt(objfile, np.array([wave,flux]).T)
        columns = ["LAMBDA", "FLUX"]

    # write wrange include file if given list
    if type(wrange_include)==type(list()):
        wrifile = dirname+"/wrange_incl.txt"
        fp=open(wrifile, "w")
        for wr in wrange_include:
            fp.write(str(wr[0])+" "+str(wr[1])+"\n")
        fp.write("\n")
        fp.close()
        wrange_include = wrifile
    # write wrange exclude file if given list
    if type(wrange_exclude)==type(list()):
        wrefile = dirname+"/wrange_excl.txt"
        fp=open(wrefile, "w")
        for wr in wrange_exclude:
            fp.write(str(wr[0])+" "+str(wr[1])+"\n")
        fp.write("\n")
        fp.close()
        wrange_exclude = wrefile
    # write prange exclude file if given list
    if type(prange_exclude)==type(list()):
        prefile = dirname+"/prange_excl.txt"
        fp=open(prefile, "w")
        for pr in prange_exclude:
            fp.write(str(pr[0])+" "+str(pr[1])+"\n")
        fp.write("\n")
        fp.close()
        prange_exclude = prefile

    # write the par file 
    parfile = createParFile(basedir, objfile, None, outname, outdir=dirname, \
                            columns=columns, \
                            obsdate=obsdate, obsdate_key=None, \
                            utc=utc, utc_key=None, \
                            telalt=telalt, telalt_key=None, \
                            rhum=rhum, rhum_key=None, \
                            pres=pres, pres_key=None, \
                            temp=temp, temp_key=None, \
                            m1temp=m1temp, m1temp_key=None, \
                            geoelev=geoelev, geoelev_key=None, \
                            longitude=longitude, longitude_key=None, \
                            latitude=latitude, latitude_key=None,\
                            slitw=slitw, slitw_key=None, \
                            pixsc=pixsc, pixsc_key=None, \
                            trans=trans, default_error=default_error, wlgtomicron=wlgtomicron, \
                            vac_air=vac_air, wrange_include=wrange_include, \
                            wrange_exclude=wrange_exclude, prange_exclude=prange_exclude, plot_range=plot_range, \
                            ftol=ftol, xtol=xtol, list_molec=list_molec, \
                            fit_molec=fit_molec, flux_unit=flux_unit, relcol=relcol, relres_box=relres_box, \
                            fit_back=fit_back, fit_cont=fit_cont, cont_n=cont_n, cont_const=cont_const, \
                            fit_wlc=fit_wlc, wlc_n=wlc_n, wlc_const=wlc_const, \
                            fit_res_gauss=fit_res_gauss, res_gauss=res_gauss, \
                            kernfac=kernfac, varkern=varkern, plot_type=plot_type)
    # run skycor
    runMolecfit(parfile, quiet=quiet)
    runCalctrans(parfile, quiet=quiet)

    # collect molecfit results
    try:
        MFresults = ap.Table()
        MFresults.read(molecfitfile, type='fits', verbose=False)
    except:
        warn.warn("Failed to open molecfit results file.")
        MFresults=None
    # collect calctans results
    try:
        CTresults = ap.Table()
        CTresults.read(calctransfile, type='fits', verbose=False)
    except:
        warn.warn("Failed to open calctrans results file.")

    if doPlot: plotMolecfitResults(MFresults, CTresults, displPC=99.0)
    pdb.set_trace()

    # get new wave axis if asked
    if calcNewWave & (type(MFresults)!=type(None)):
        coefs = readWaveInfo(dirname, resFile)
        if coefs!=None:
            x = 2.0*(wave-wave.min())/(wave.max()-wave.min()) - 1.0
            newx = cheby(x,coefs)
            newwave = 0.5*(newx + 1.0)*(wave.max()-wave.min())+wave.min()
        else:
            newwave = wave
        MFresults.add_column('lambda_new', newwave)

    if cleanTmp:
        sp.call(["rm", "-rf", dirname])

    return MFresults,CTresults

def readWaveInfo(dir, resFile, verbose=False):
    """
    RH 20/11/15

    Helper function to read in the (recalibrated) wavelength info in the molecfit RES file.

    """

    waveFile = "WAVEINFO_"+resFile # init
    # extract wave info from res file
    fp = open(dir+waveFile, "w")
    cat = sp.Popen(('cat', dir+resFile), stdout=sp.PIPE)
    sp.call(('grep', '^w'), stdin=cat.stdout, stdout=fp)
    cat.wait()
    fp.close()
    # read in wave info
    fp2=open(dir+waveFile, "r")
    info= []
    for line in fp2:
        info.append(line.split())
    if len(info)==0:
        if verbose: print "WARNING: no wave info found!"
        coefs = None
    else:
        # extract cheby coefs, sorted
        info = np.array(info)
        sa = info[:,1].argsort()
        coefs = np.array(info[sa,3],dtype=np.float)

    fp2.close()
    
    return coefs

def plotMolecfitResults(MFresults, CTresults=None, figsize=(20,12), displPC=99.0):
    """
    RH 5/1/16

    Plot the molecfit results: original and fit for all fit regions and a zoom on just the fitted regions underneath
    
    """

    pl.figure(figsize=figsize)
    pl.subplot(211)
    pl.plot(MFresults['lambda'], MFresults['flux'], "k-", alpha=0.7)
    pl.plot(MFresults['lambda'], MFresults['mtrans'], "b-")
    if type(CTresults)!=type(None):
        pl.plot(CTresults['lambda'], CTresults['flux'], "k:", alpha=0.7)
        pl.plot(CTresults['lambda'], CTresults['mtrans'], "c-")
        pl.plot(CTresults['lambda'], CTresults['flux']/CTresults['mtrans'], "g-")
    # find model regions
    ureg = list(set(list(MFresults['mrange'])))
    ureg.remove([0])
    nreg=len(ureg)
    for r in ureg:
        loc = np.where(MFresults['mrange']==r)
        pl.plot(MFresults['lambda'][loc], MFresults['mflux'][loc], "r-")
    pl.axis([MFresults['lambda'].min(), MFresults['lambda'].max(), 0.0, np.percentile(MFresults['flux'], displPC)/displPC*100.0])
    for nri,r in enumerate(ureg):
        pl.subplot(2,nreg,nri+nreg+1)
        loc = np.where(MFresults['mrange']==r)
        pl.plot(MFresults['lambda'][loc], MFresults['flux'][loc], "k-", alpha=0.7)
        pl.plot(MFresults['lambda'][loc], MFresults['mflux'][loc], "r-")
        pl.plot(MFresults['lambda'][loc], MFresults['flux'][loc]-MFresults['mflux'][loc], "g.")
        
        
def cheby(x, coefs, autonorm=True):
    """
    RH 20/11/15

    Generate Chebyshev Polynomials of the first kind, multiplied by a coefficient array and summed.
    AUTONORM makes sure that we evaluate the Chebys with x normalised between -1 and 1

    """

    # autonorm x array
    if autonorm:
        xmin = x.min()
        xmax = x.max()
        xx   = 2.0*(x-xmin)/(xmax-xmin) - 1.0
    else:
        xx=x

    # check I'm doing what i think I'm doing
    if len(coefs.shape) > 1:
        raise IOError("Array of coeffients has more than one dimension?")

    # shortcut using numpy
    poly = np.polynomial.chebyshev.chebval(xx,coefs)
    
    return poly

    
def test_example(input="examples/input/xshoo_vis_1.fits", installdir=molecfitinstall):
    """
    RH 19/11/15

    Test this code by running the standard molecfit example

    """
    # test if outdir already exists and quit if so
    tag = input[input.rfind("/")+1:input.rfind(".fits")]
    outdir=installdir+'examples/test_'+tag
    if os.path.isdir(outdir):
        raise IOError("Output directory already exists. Please delete and run again.")
    # create par file
    output = "OUT_"+tag
    
    parfile = createParFile(installdir, input, None, output, \
                            outdir=outdir, columns=["LAMBDA", "FLUX", "ERR", "QUAL"], \
                            trans=True, wlgtomicron=1e-3, vac_air='air', \
                            wrange_include=installdir+"examples/config/include_xshoo_vis.dat", \
                            wrange_exclude=None, \
                            prange_exclude=installdir+"examples/config/exclude_xshoo_vis.dat", \
                            list_molec=["H2O", "O2"], fit_molec=[True,False], relcol=[1.0,1.0], \
                            flux_unit=False, fit_back=False, fit_cont=True, cont_n=3, cont_const=1.0, \
                            fit_wlc=True, wlc_n=0, wlc_const=0.0, fit_res_gauss=True, res_gauss=1.0, \
                            kernfac=30.0, varkern=True, slitw=0.9, slitw_key=None, pixsc=0.16, pixsc_key=None)
    runMolecfit(parfile, quiet=False)
    if os.path.isfile(outdir+'/'+output+'_fit.fits'):
        print "Run successful"
    else:
        print "Run FAILED."


def test_ascii(objfile=molecfitinstall+"examples/input/xshoo_vis_1.fits", \
               installdir=molecfitinstall, calcNewWave=False, plot_type="XP", ftol=1e-2, xtol=1e-2):
    """
    RH 19/11/15

    Test the ascii wrapper around a std example

    """

    otab = ap.Table()
    otab.read(objfile, type="fits", verbose=False)

    # get some info from the header
    ofh = pf.open(objfile)
    hdr = ofh[1].header
    obsdate = hdr['MJD-OBS']
    utc     = hdr['UTC']
    telalt  = hdr['ESO TEL ALT']
    rhum    = hdr['ESO TEL AMBI RHUM']
    pres    = hdr['ESO TEL AMBI PRES START']
    temp    = hdr['ESO TEL AMBI TEMP']
    m1temp  = hdr['ESO TEL TH M1 TEMP']
    geoelev = hdr['ESO TEL GEOELEV']
    longitude=hdr['ESO TEL GEOLON']
    latitude= hdr['ESO TEL GEOLAT']
    slitw   = 0.9
    pixsc   = 0.16
    # close the fits file
    ofh.close()

    # pass info to ascii wrapper
    results = asciiMolecfitWrapper(otab['LAMBDA'], otab['FLUX'], \
                                   obsdate, utc, telalt, rhum, pres, temp, m1temp, geoelev, longitude, latitude, \
                                   err=otab['ERR'], mask=otab['QUAL'], slitw=slitw, pixsc=pixsc, \
                                   trans=True, wlgtomicron=1e-3, vac_air='air', \
                                   wrange_include=installdir+"examples/config/include_xshoo_vis.dat", \
                                   wrange_exclude=None, xtol=xtol, ftol=ftol, \
                                   prange_exclude=installdir+"examples/config/exclude_xshoo_vis.dat", \
                                   list_molec=["H2O", "O2"], fit_molec=[True,False], relcol=[1.0,1.0], \
                                   flux_unit=False, fit_back=False, fit_cont=True, cont_n=6, cont_const=1.0, \
                                   fit_wlc=True, wlc_n=1, wlc_const=0.0, fit_res_gauss=True, res_gauss=0.9, relres_box=0.0, \
                                   kernfac=10.0, varkern=True, calcNewWave=calcNewWave, cleanTmp=False, basedir=installdir, \
                                   plot_type=plot_type)

    
    print 'Test successful!'
    pdb.set_trace()
    
    

def KMOS1DFITSspec(fitsfile, extname="IFU.22.DATA", \
                   wrange_include=[[1.11,1.16],[1.203,1.209],[1.26,1.275],[1.309,1.34]], \
                   wrange_exclude=None, prange_exclude=[[0,90],[1970,2200]], \
                   plot_type="N", ftol=1e-2, xtol=1e-2, \
                   list_molec=['H2O', 'O2', 'CH4', 'CO2'], fit_molec=[True, True, True, True], relcol=[1.0,1.0,1.0,1.0],
                   cont_n=3, cont_const=1.0, wlc_n=1, wlc_const=0.0, varkern=False, \
                   installdir=molecfitinstall, \
                   outdir="/Volumes/data4/rcwh/data7/KMOS/XMMUJ2235/sci/2013-10-30/KMOSspec", quiet=False, doPlot=True):
    """
    RH 8/1/15

    Run molecfit & calctrans on KMOS spectrum

    """

    outfile=fitsfile[fitsfile.rfind("/"):fitsfile.rfind(".fits")]+"_out"
    molecfitfile = outdir+outfile+"_fit.fits"
    calctransfile = outdir+outfile+"_tac.fits"

    # make dir if doesn't exist
    if (not os.path.isdir(outdir)):
        os.mkdir(outdir)

    # write wrange include file if given list
    if type(wrange_include)==type(list()):
        wrifile = outdir+"/wrange_incl.txt"
        fp=open(wrifile, "w")
        for wr in wrange_include:
            fp.write(str(wr[0])+" "+str(wr[1])+"\n")
        fp.write("\n")
        fp.close()
        wrange_include = wrifile
    # write wrange exclude file if given list
    if type(wrange_exclude)==type(list()):
        wrefile = outdir+"/wrange_excl.txt"
        fp=open(wrefile, "w")
        for wr in wrange_exclude:
            fp.write(str(wr[0])+" "+str(wr[1])+"\n")
        fp.write("\n")
        fp.close()
        wrange_exclude = wrefile
    # write prange exclude file if given list
    if type(prange_exclude)==type(list()):
        prefile = outdir+"/prange_excl.txt"
        fp=open(prefile, "w")
        for pr in prange_exclude:
            fp.write(str(pr[0])+" "+str(pr[1])+"\n")
        fp.write("\n")
        fp.close()
        prange_exclude = prefile

    parfile = createParFile(installdir, fitsfile, None, outfile, outdir=outdir, trans=True, parfile=None, \
                            wrange_include=wrange_include, wrange_exclude=None, prange_exclude=prange_exclude, \
                            columns=["NULL", extname, "NULL", "NULL"], \
                            default_error=0.01, wlgtomicron=1.0, vac_air='vac', \
                            plot_type=plot_type, plot_range=False, ftol=ftol, xtol=xtol, \
                            list_molec=list_molec, fit_molec=fit_molec, relcol=relcol, \
                            flux_unit=0, fit_back=True, telback=0.1, fit_cont=True, cont_n=cont_n, cont_const=cont_const, \
                            fit_wlc=True, wlc_n=wlc_n, wlc_const=wlc_const, fitresbox=False, kernmode=False, \
                            fit_res_box=False, relres_box=1.0, fit_res_gauss=True, res_gauss=1.0, fit_res_lorentz=False, res_lorentz=0.5, \
                            kernfac=10.0, varkern=varkern, kernel_file=None, obsdate=None, obsdate_key="MJD-OBS", \
                            utc=None, utc_key="UTC", telalt=None, telalt_key="ESO TEL ALT", rhum=None, \
                            rhum_key="ESO TEL AMBI RHUM", pres=None, pres_key="ESO TEL AMBI PRES START", \
                            temp=None, temp_key="ESO TEL AMBI TEMP", m1temp=None, m1temp_key="ESO TEL TH M1 TEMP", \
                            geoelev=None, geoelev_key="ESO TEL GEOELEV", longitude=None, \
                            longitude_key="ESO TEL GEOLON", latitude=None, latitude_key="ESO TEL GEOLAT", \
                            slitw=0.2, slitw_key=None, pixsc=0.2, pixsc_key=None, \
                            gdas_prof="auto", layers=1, emix=5.0, pwv=-1.0)

    
    # run skycor
    runMolecfit(parfile, quiet=quiet)
    runCalctrans(parfile, quiet=quiet)
    
    # collect molecfit results
    try:
        MFresults = ap.Table()
        MFresults.read(molecfitfile, type='fits', verbose=False)
    except:
        warn.warn("Failed to open molecfit results file.")
        MFresults=None
    # collect calctans results
    try:
        CTresults = ap.Table()
        CTresults.read(calctransfile, type='fits', verbose=False)
    except:
        warn.warn("Failed to open calctrans results file.")

    if doPlot: plotMolecfitResults(MFresults, CTresults, displPC=99.0)

    return MFresults, CTresults

    
