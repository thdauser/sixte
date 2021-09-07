import subprocess
import sys
import os
import astropy.io.fits as fits

# TODO: automate output (test and ref) file naming
# TODO: runsixt sim with off-axis psf and sources, slew, and vignetting
# TODO: makeref argument for test routines

# PYTHON Version control
SYSVERMAJOR_MIN = 3
SYSVERMINOR_MIN = 6


def check_pythonversion(vmajor=SYSVERMAJOR_MIN,
                        vminor=SYSVERMINOR_MIN):
    assert (sys.version_info.major >= vmajor and sys.version_info.minor >= vminor), \
        ' *** error : Calling with python version {}.{}, but >= {}.{} is required!' \
            .format(sys.version_info.major, sys.version_info.minor, vmajor, vminor)

def enable_sixte_slurm_use():
    os.environ['HEADASNOQUERY'] = ''
    os.environ['HEADASPROMPT'] = '/dev/null'

check_pythonversion()


class defvar:

    def __init__(self):
        self.fname_pholist = "pho.fits"
        self.fname_implist = "imp.fits"
        self.fname_rawlist = "raw.fits"
        self.fname_evtlist = "evt.fits"
        self.fname_expmap = "expmap.fits"
        self.fname_expmap_fov = "expmap_fov.fits"  # for fov_diameter testing
        self.fname_img = "img.fits"

        self.inst = "dummy_inst"
        self.mode = "dummy_mode"
        self.filt = "dummy_filt"
        self.miss = "dummy_miss"

        self.xml_path = "data/dummy/"
        self.xml_name = "default_inst.xml"
        self.xml_name_eventmode = "default_eventreadout.xml"

        self.xml = self.xml_path + self.xml_name
        self.xml_eventmode = self.xml_path + self.xml_name_eventmode

        self.RA = 0.0
        self.Dec = 0.0

        self.expos = 100
        self.dt = 1

        self.coordsys = 0
        self.projection = "TAN"
        # Standard values taken from default_inst.xml
        self.naxis1 = 9
        self.naxis2 = 9
        self.cunit1 = "deg"
        self.cunit2 = "deg"
        self.crval1 = 0
        self.crval2 = 0
        self.crpix1 = 5
        self.crpix2 = 5
        self.cdelt1 = -100e-6
        self.cdelt2 = 100e-6

        self.simput_name = "dummy.simput"
        self.simput_dir = "data/dummy/"
        self.simput = self.simput_dir + self.simput_name

        self.attitude = "data/dummy/dummy_e2e.att"
        self.vign = "data/dummy/dummy_vign.fits"

        self.seed = 0

        self.background = "no"

        self.log = "logfile.log"

        self.ref_path = "data/refdata/"

        self.prefix_refdata = "ref"
        self.prefix_testdata = "test"
        self.prefix_dummy = "dummy"

        self.mjdref = 55000
        self.tstart = 0

        def simput(self):
            self.simput = self.simput_dir + self.simput_name

        def set_simput_dir(self, simput_dir):
            self.simput_dir = simput_dir


class defpath(defvar):

    def __init__(self, subtestname='default', testname=None):

        defvar.__init__(self)

        self.log = 'logfile_' + subtestname + '.log'

        self.subtestname = subtestname
        if testname is None:
            self.testname = self.get_testname_from_path()
        else:
            self.testname = testname
        self.fullname = self.testname + "_" + self.subtestname

        self.testname_pholist = self.get_testfile_prefix(subtestname) + self.fname_pholist
        self.testname_implist = self.get_testfile_prefix(subtestname) + self.fname_implist
        self.testname_rawlist = self.get_testfile_prefix(subtestname) + self.fname_rawlist
        self.testname_evtlist = self.get_testfile_prefix(subtestname) + self.fname_evtlist
        self.testname_expmap = self.get_testfile_prefix(subtestname) + self.fname_expmap
        self.testname_img = self.get_testfile_prefix(subtestname) + self.fname_img
        
        self.refname_pholist = self.get_reffile_prefix(subtestname) + self.fname_pholist
        self.refname_implist = self.get_reffile_prefix(subtestname) + self.fname_implist
        self.refname_rawlist = self.get_reffile_prefix(subtestname) + self.fname_rawlist
        self.refname_evtlist = self.get_reffile_prefix(subtestname) + self.fname_evtlist
        # There are two exposure maps (fov_diameter=-1, 1/30)
        self.refname_expmap = self.get_reffile_prefix(subtestname) + self.fname_expmap
        self.refname_expmap_fov = self.get_reffile_prefix(subtestname) + self.fname_expmap_fov
        self.refname_img = self.get_reffile_prefix(subtestname) + self.fname_img

    def get_testname_from_path(self):
        """
        Get name of test from last/working directory. Expects
        that test is run in this directory!
        """
        return os.getcwd().split('/')[-1]

    def get_outfile_prefix(self, prefix, subtestname):
        """
        """
        return f"{prefix}_{self.testname}_{subtestname}_"

    def get_testfile_prefix(self, subtestname=None):
        """ Return automatically generated prefix for test output files"""
        if subtestname == None:
            subtestname = self.subtestname
        return self.get_outfile_prefix(self.prefix_testdata, subtestname)

    def get_reffile_prefix(self, subtestname=None):
        """ Return automatically generated prefix for reference data"""
        if subtestname == None:
            subtestname = self.subtestname
        return self.get_outfile_prefix(self.prefix_refdata, subtestname)

    def get_pointing_string(self, ra, dec, attitude):
        if (attitude == 0):
            return f"RA={ra} Dec={dec}"
        else:
            return f"Attitude={attitude}"

    def get_exposure_string(self, tstart, exposure, gtifile):
        if (gtifile == 0 or gtifile is None):
            return f"TSTART={tstart} Exposure={exposure}"
        else:
            return f"GTIFile={gtifile}"


# initialize one instance of defvar
STDTEST = defpath()


def append_to_log_file(logfile, cmd):
    fp = open(logfile, 'a+')
    fp.write(f" *** RUNNING: '{cmd}'\n" + "=" * 80 + "\n\n")
    fp.close()


def phogen(phlist, xmlfile, simput, expos=STDTEST.expos,
           ra=STDTEST.RA, dec=STDTEST.Dec,
           instrument=STDTEST.inst,
           attitude=0,
           mission=STDTEST.miss,
           mode=STDTEST.mode,
           seed=STDTEST.seed,
           clobber="yes",
           prefix="sim_",
           logfile=None,
           test=False):
    if test:
        os.environ["SIXTE_USE_PSEUDO_RNG"] = "1"

    pointing_string = STDTEST.get_pointing_string(ra, dec, attitude)

    str = f"""phogen \
        PhotonList={prefix}{phlist} \
        Instrument={instrument} \
        Mission={mission} \
        Mode={mode} \
        XMLFile={xmlfile} \
        {pointing_string} \
        Simput={simput} \
        Exposure={expos} \
        Seed={seed} \
        clobber={clobber}"""

    if logfile:
        fp = open(logfile, 'w')
        fp.write(f" *** RUNNING: '{str}'\n" + "=" * 80 + "\n\n")
        fp.close()
        str = f"{str} >> {logfile}"

    ret_val = subprocess.run(str, shell=True,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    if os.getenv('SIXTE_USE_PSEUDO_RNG'):
        del os.environ['SIXTE_USE_PSEUDO_RNG']

    return ret_val


def phoimg(phlist, implist, xmlfile, expos=STDTEST.expos,
           ra=STDTEST.RA, dec=STDTEST.Dec,
           instrument=STDTEST.inst,
           mission=STDTEST.miss,
           mode=STDTEST.mode,
           seed=STDTEST.seed,
           mjdref=STDTEST.mjdref,
           tstart=STDTEST.tstart,
           clobber="yes",
           logfile=None,
           test=False):
    if test:
        os.environ["SIXTE_USE_PSEUDO_RNG"] = "1"

    str = f"""phoimg \
        PhotonList={phlist} \
        ImpactList={implist} \
        Instrument={instrument} \
        Mission={mission} \
        Mode={mode} \
        XMLFile={xmlfile} \
        RA={ra} Dec={dec} \
        MJDREF={mjdref} \
        TSTART={tstart} \
        Exposure={expos} \
        Seed={seed} \
        clobber={clobber}"""

    if logfile:
        fp = open(logfile, 'w')
        fp.write(f" *** RUNNING: '{str}'\n" + "=" * 80 + "\n\n")
        fp.close()
        str = f"{str} >> {logfile}"

    ret_val = subprocess.run(str, shell=True,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    if os.getenv('SIXTE_USE_PSEUDO_RNG'):
        del os.environ['SIXTE_USE_PSEUDO_RNG']

    return ret_val


def gendetsim(implist, rawdata, xmlfile, expos=STDTEST.expos,
              instrument=STDTEST.inst,
              mission=STDTEST.miss,
              mode=STDTEST.mode,
              seed=STDTEST.seed,
              mjdref=STDTEST.mjdref,
              tstart=STDTEST.tstart,
              clobber="yes",
              background=STDTEST.background,
              logfile=None,
              test=False):
    if test:
        os.environ["SIXTE_USE_PSEUDO_RNG"] = "1"

    str = f"""gendetsim \
        ImpactList={implist} \
        RawData={rawdata} \
        Instrument={instrument} \
        Mission={mission} \
        Mode={mode} \
        XMLFile={xmlfile} \
        MJDREF={mjdref} \
        TSTART={tstart} \
        Exposure={expos} \
        Seed={seed} \
        Background={background} \
        clobber={clobber}"""

    if logfile:
        fp = open(logfile, 'w')
        fp.write(f" *** RUNNING: '{str}'\n" + "=" * 80 + "\n\n")
        fp.close()
        str = f"{str} >> {logfile}"

    ret_val = subprocess.run(str, shell=True,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    if os.getenv('SIXTE_USE_PSEUDO_RNG'):
        del os.environ['SIXTE_USE_PSEUDO_RNG']

    return ret_val


def runsixt(xmlfile,
            expos=STDTEST.expos,
            evtfile=STDTEST.fname_evtlist,
            prefix=STDTEST.prefix_dummy,
            rawdata=0,
            implist=0,
            ra=STDTEST.RA, dec=STDTEST.Dec,
            instrument=STDTEST.inst,
            attitude=0,
            simput=STDTEST.simput,
            mission=STDTEST.miss,
            mode=STDTEST.mode,
            seed=STDTEST.seed,
            mjdref=STDTEST.mjdref,
            tstart=STDTEST.tstart,
            clobber="yes",
            background=STDTEST.background,
            precmd="",
            logfile=None,
            test=False,
            ):
    if test:
        os.environ["SIXTE_USE_PSEUDO_RNG"] = "1"

    pointing_string = STDTEST.get_pointing_string(ra, dec, attitude)
    impfile_string = ("" if implist == 0 else f"ImpactList={implist} ")
    rawdata_string = ("" if implist == 0 else f"RawData={rawdata} ")

    str = precmd + f"""runsixt {pointing_string} \
Prefix={prefix} \
{impfile_string} {rawdata_string} \
EvtFile={evtfile}  \
Instrument={instrument} \
Mission={mission} \
Mode={mode} \
XMLFile={xmlfile} \
MJDREF={mjdref} \
Simput={simput}  \
TSTART={tstart} \
Exposure={expos} \
Seed={seed} \
Background={background} \
clobber={clobber}"""

    if logfile:
        fp = open(logfile, 'w')
        fp.write(f" *** RUNNING: '{str}'\n" + "=" * 80 + "\n\n")
        fp.close()
        str = f"{str} >> {logfile}"

    ret_val = subprocess.run(str, shell=True, executable='/bin/bash',
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    if os.getenv('SIXTE_USE_PSEUDO_RNG'):
        del os.environ['SIXTE_USE_PSEUDO_RNG']

    return ret_val


def imgev(evtfile=STDTEST.fname_evtlist, 
          image=STDTEST.fname_img,
          projection=STDTEST.projection,
          naxis1=STDTEST.naxis1, naxis2=STDTEST.naxis2,
          cunit1=STDTEST.cunit1, cunit2=STDTEST.cunit2,
          crval1=STDTEST.crval1, crval2=STDTEST.crval2,
          crpix1=STDTEST.crpix1, crpix2=STDTEST.crpix2,
          cdelt1=STDTEST.cdelt1, cdelt2=STDTEST.cdelt2,
          chatter=3, clobber="yes", history="true",
          precmd="", logfile=None, test=False):

    if test:
        os.environ["SIXTE_USE_PSEUDO_RNG"] = "1"
        
    str = precmd + f"""imgev EvtFile={evtfile} Image={image} \
Projection={projection} \
NAXIS1={naxis1} NAXIS2={naxis2} \
CUNIT1={cunit1} CUNIT2={cunit2} \
CRVAL1={crval1} CRVAL2={crval2} \
CRPIX1={crpix1} CRPIX2={crpix2} \
CDELT1={cdelt1} CDELT2={cdelt2} \
chatter={chatter} clobber={clobber} history={history}"""

    if logfile:
        fp = open(logfile, 'w')
        fp.write(f" *** RUNNING: '{str}'\n" + "=" * 80 + "\n\n")
        fp.close()
        str = f"{str} >> {logfile}"

    ret_val = subprocess.run(str, shell=True, executable='/bin/bash',
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    if os.getenv('SIXTE_USE_PSEUDO_RNG'):
        del os.environ['SIXTE_USE_PSEUDO_RNG']

    return ret_val


def exposure_map(exposuremap, xmlfile,
                 ra=STDTEST.RA, dec=STDTEST.Dec,
                 vignetting=STDTEST.vign,
                 tstart=STDTEST.tstart,
                 timespan=STDTEST.expos,
                 dt=STDTEST.dt,
                 fov_diameter=-1,
                 coordsystem=STDTEST.coordsys, projection=STDTEST.projection,
                 naxis1=STDTEST.naxis1, naxis2=STDTEST.naxis2,
                 cunit1=STDTEST.cunit1, cunit2=STDTEST.cunit2,
                 crval1=STDTEST.crval1, crval2=STDTEST.crval2,
                 crpix1=STDTEST.crpix1, crpix2=STDTEST.crpix2,
                 cdelt1=STDTEST.cdelt1, cdelt2=STDTEST.cdelt2,
                 clobber="yes",
                 logfile=None,
                 test=False):
    if test:
        os.environ["SIXTE_USE_PSEUDO_RNG"] = "1"

    str = f"""exposure_map \
        RA={ra} Dec={dec} \
        Vignetting={vignetting} \
        XMLFile={xmlfile} \
        Exposuremap={exposuremap} \
        TSTART={tstart} \
        timespan={timespan} \
        dt={dt} \
        fov_diameter={fov_diameter} \
        CoordinateSystem={coordsystem} projection_type={projection} \
        NAXIS1={naxis1} NAXIS2={naxis2} \
        CUNIT1={cunit1} CUNIT2={cunit2} \
        CRVAL1={crval1} CRVAL2={crval2} \
        CRPIX1={crpix1} CRPIX2={crpix2} \
        CDELT1={cdelt1} CDELT2={cdelt2} \
        clobber={clobber}
        """

    if logfile:
        fp = open(logfile, 'w')
        fp.write(f" *** RUNNING: '{str}'\n" + "=" * 80 + "\n\n")
        fp.close()
        str = f"{str} >> {logfile}"

    ret_val = subprocess.run(str, shell=True, executable='/bin/bash',
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    if os.getenv('SIXTE_USE_PSEUDO_RNG'):
        del os.environ['SIXTE_USE_PSEUDO_RNG']

    return ret_val


def erosim(xmlfile1, xmlfile2, xmlfile3, xmlfile4, xmlfile5, xmlfile6, xmlfile7,
           expos=STDTEST.expos,
           evtfile=STDTEST.fname_evtlist,
           prefix=STDTEST.prefix_dummy,
           rawdata=None, implist=None,
           ra=STDTEST.RA, dec=STDTEST.Dec,
           attitude=None, gtifile=None,
           simput=STDTEST.simput,
           seed=STDTEST.seed,
           background=STDTEST.background,
           mjdref=STDTEST.mjdref, tstart=STDTEST.tstart, dt=1.0,
           clobber="yes", chatter=3, precmd="",
           logfile=None, test=False
           ):
    if test:
        os.environ["SIXTE_USE_PSEUDO_RNG"] = "1"

    pointing_string = STDTEST.get_pointing_string(ra, dec, attitude)
    expos_string = STDTEST.get_exposure_string(tstart, expos, gtifile)
    impfile_string = f"ImpactList={implist}" if implist else ""
    rawdata_string = f"RawData={rawdata}" if rawdata else ""

    xmlfilelist=[xmlfile1, xmlfile2, xmlfile3, xmlfile4, xmlfile5, xmlfile6, xmlfile7]
    xmls = " ".join([f"XMLFile{tm+1}={xmlfile}" for tm, xmlfile in enumerate(xmlfilelist)])

    cmd = precmd + f"""erosim Prefix={prefix} \
EvtFile={evtfile} \
{xmls} \
{pointing_string} {expos_string} \
{impfile_string} {rawdata_string} \
MJDREF={mjdref} \
Simput={simput}  \
dt={dt} \
Seed={seed} \
Background={background} \
chatter={chatter} clobber={clobber}"""

    if logfile:
        append_to_log_file(logfile, cmd)
        cmd = f"{cmd} >> {logfile}"

    ret_val = subprocess.run(cmd, shell=True, executable='/bin/bash',
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    check_returncode(ret_val, "erosim")

    if os.getenv('SIXTE_USE_PSEUDO_RNG'):
        del os.environ['SIXTE_USE_PSEUDO_RNG']

    return ret_val


def ero_vis(attitude, gtifile, simput=False,
            ra=STDTEST.RA, dec=STDTEST.Dec, tstart=STDTEST.tstart,
            exposure=STDTEST.expos, dt=STDTEST.dt, visibility_range=1.02,
            chatter=3, clobber="yes", history="true",
            precmd="", logfile=None, test=False):
    if test:
        os.environ["SIXTE_USE_PSEUDO_RNG"] = "1"

    pointing_string = f"Simput={simput}" if simput else f"SrcRA={ra} SrcDec={dec}"

    cmd = precmd + f"""ero_vis \
Attitude={attitude} \
{pointing_string} \
GTIfile={gtifile} \
TSTART={tstart} \
Exposure={exposure} dt={dt} \
visibility_range={visibility_range} \
chatter={chatter} clobber={clobber} history={history}"""

    if logfile:
        append_to_log_file(logfile, cmd)
        cmd = f"{cmd} >> {logfile}"

    ret_val = subprocess.run(cmd, shell=True, executable='/bin/bash',
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    check_returncode(ret_val, "ero_vis")
    
    if os.getenv('SIXTE_USE_PSEUDO_RNG'):
        del os.environ['SIXTE_USE_PSEUDO_RNG']

    return ret_val


def ero_calevents(evtfile, eroevtfile, ccdnr, attitude=None,
                  ra=STDTEST.RA, dec=STDTEST.Dec, refra=0.0, refdec=0.0,
                  projection="SIN", rollangle=0.0, usepha=0, seed=STDTEST.seed,
                  chatter=3, clobber="yes", history="true",
                  precmd="", logfile=None, test=False):
    if test:
        os.environ["SIXTE_USE_PSEUDO_RNG"] = "1"

    attitude_string = f"Attitude={attitude}" if attitude else ""

    cmd = precmd + f"""ero_calevents \
EvtFile={evtfile} eroEvtFile={eroevtfile} \
CCDNr={ccdnr} Projection={projection} \
usepha={usepha} Seed={seed} \
RA={ra} Dec={dec} RefRA={refra} RefDec={refdec} rollangle={rollangle} {attitude_string} \
chatter={chatter} clobber={clobber} history={history}"""

    if logfile:
        append_to_log_file(logfile, cmd)
        cmd = f"{cmd} >> {logfile}"

    ret_val = subprocess.run(cmd, shell=True, executable='/bin/bash',
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    check_returncode(ret_val, f"ero_calevents on TM{ccdnr} ")

    if os.getenv('SIXTE_USE_PSEUDO_RNG'):
        del os.environ['SIXTE_USE_PSEUDO_RNG']

    return ret_val


def makelc(evtfile, lightcurve, tstart, length, dt=1.0,
           emin=-1.0, emax=-1.0, chanmin=-1, chanmax=-1,
           gtifile=None, chatter=3, clobber="yes", history="true",
           precmd="", logfile=None, test=False):

    if test:
        os.environ["SIXTE_USE_PSEUDO_RNG"] = "1"

    gtistring = f"GTIFile={gtifile}" if gtifile else ""

    cmd = precmd + f"""makelc \
EvtFile={evtfile} LightCurve={lightcurve} TSTART={tstart} length={length} dt={dt} \
Emin={emin} Emax={emax} Chanmin={chanmin} Chanmax={chanmax} \
{gtistring} chatter={chatter} clobber={clobber} history={history}"""

    if logfile:
        append_to_log_file(logfile, cmd)
        cmd = f"{cmd} >> {logfile}"

    ret_val = subprocess.run(cmd, shell=True, executable='/bin/bash',
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    check_returncode(ret_val, "makelc")

    if os.getenv('SIXTE_USE_PSEUDO_RNG'):
        del os.environ['SIXTE_USE_PSEUDO_RNG']
        
    return ret_val


def check_returncode(ret_val, tool):
    if (ret_val.returncode == 0):
        print(f"{tool}: run SUCCESSFUL")
    else:
        print('{}: run FAILED with ReturnCode {}'.format(tool, ret_val.returncode))
        out = ("{}".format(ret_val.stdout)).replace('\\n', '\n')
        out = out.replace('\'', '')
        out = out.replace('b', '')

        err = ("{}".format(ret_val.stderr)).replace('\\n', '\n')
        err = err.replace('\'', '')
        err = err.replace('b', '')

        print(f"{tool}: stdout: {out}")
        print(f"{tool}: stderr: {err}")
        exit(ret_val.returncode)


def check_returnvalue(ret_val, tool):
    if (ret_val == 0):
        print(f"{tool}: run SUCCESSFUL")
    else:
        print('{}: run FAILED with ReturnCode {}'.format(tool, ret_val))
        exit(ret_val)


def check_fdiff(file1, file2, tool):
    output = fits.FITSDiff(file1, file2, ignore_comments=['*'], ignore_keywords=['DATE', 'HISTORY', 'CHECKSUM'])

    if (output.identical):
        print(f'{tool}: comparison to refdata SUCCESSFUL')
    elif (output.identical == False):
        print(output.report())
        print(f'*** error *** {tool}: comparison to refdata FAILED')
        exit(1)


def clean_output():
    subprocess.call('rm -f *.fits', shell=True)
