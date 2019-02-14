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

def check_pythonversion( vmajor=SYSVERMAJOR_MIN,
                         vminor=SYSVERMINOR_MIN ):
    assert (sys.version_info.major >= vmajor and sys.version_info.minor >= vminor),\
        ' *** error : Calling with python version {}.{}, but >= {}.{} is required!'\
        .format(sys.version_info.major,sys.version_info.minor,vmajor,vminor)

check_pythonversion()


class defvar:

    def __init__(self):
        
        self.fname_pholist="pho.fits"
        self.fname_implist="imp.fits"
        self.fname_rawlist="raw.fits"
        self.fname_evtlist="evt.fits"

        self.inst="dummy_inst"
        self.mode="dummy_mode"
        self.filt="dummy_filt"
        self.miss="dummy_miss"

        self.xml_path="data/dummy/"
        self.xml_name="default_inst.xml"
        self.xml_name_eventmode="default_eventreadout.xml"

        self.xml=self.xml_path+self.xml_name
        self.xml_eventmode=self.xml_path+self.xml_name_eventmode

        self.RA=0.0
        self.Dec=0.0

        self.expos=100

        self.simput_name="dummy.simput"
        self.simput_dir="data/dummy/"
        self.simput=self.simput_dir+self.simput_name

        self.attitude = "data/dummy/dummy_e2e.att"

        self.seed=0

        self.background="no"

        self.log="logfile.log"
        
        self.ref_path="data/refdata/"

        self.prefix_refdata="ref"
        self.prefix_testdata="test"
        self.prefix_dummy="dummy"

        self.mjdref=55000
        self.tstart=0


class defpath(defvar):
    
    def __init__(self,subtestname='default',testname=None):
        
        defvar.__init__(self)
        
        self.log = 'logfile_'+subtestname+'.log'
        
        self.subtestname = subtestname
        if testname is None:
            self.testname = self.get_testname_from_path()
        else:
            self.testname = testname
        self.fullname = self.testname+"_"+self.subtestname


        self.testname_pholist = self.get_testfile_prefix(subtestname)+self.fname_pholist
        self.testname_implist = self.get_testfile_prefix(subtestname)+self.fname_implist
        self.testname_rawlist = self.get_testfile_prefix(subtestname)+self.fname_rawlist
        self.testname_evtlist = self.get_testfile_prefix(subtestname)+self.fname_evtlist

        self.refname_pholist = self.get_reffile_prefix(subtestname)+self.fname_pholist
        self.refname_implist = self.get_reffile_prefix(subtestname)+self.fname_implist
        self.refname_rawlist = self.get_reffile_prefix(subtestname)+self.fname_rawlist
        self.refname_evtlist = self.get_reffile_prefix(subtestname)+self.fname_evtlist

    def get_testname_from_path(self):
        """
        Get name of test from last/working directory. Expects
        that test is run in this directory!
        """
        return os.getcwd().split('/')[-1]

    def get_outfile_prefix(self,prefix,subtestname):
        """
        """
        return f"{prefix}_{self.testname}_{subtestname}_"
    
    def get_testfile_prefix(self,subtestname=None):
        """ Return automatically generated prefix for test output files"""
        if subtestname == None:
            subtestname = self.subtestname
        return self.get_outfile_prefix(self.prefix_testdata,subtestname)
    
    def get_reffile_prefix(self,subtestname=None):
        """ Return automatically generated prefix for reference data"""
        if subtestname == None:
            subtestname = self.subtestname
        return self.get_outfile_prefix(self.prefix_refdata,subtestname)

    def get_pointing_string(self,ra,dec,attitude):
        if (attitude == 0):
            return f"RA={ra} Dec={dec}"
        else:
            return f"Attitude={attitude}"


# initialize one instance of defvar
STDTEST = defpath()

def phogen(phlist,xmlfile,simput,expos=STDTEST.expos,
           ra=STDTEST.RA,dec=STDTEST.Dec,
           instrument=STDTEST.inst,
           attitude=0,
           mission=STDTEST.miss,
           mode=STDTEST.mode,
           seed=STDTEST.seed,
           clobber="yes",
           prefix="sim_",
           logfile=-1,
           test=-1):

    if (test):
        os.environ["SIXTE_USE_PSEUDO_RNG"] = "1"

    pointing_string = STDTEST.get_pointing_string(ra,dec,attitude)

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

    if (logfile!=-1):
        fp  = open(logfile,'w')
        fp.write(f" *** RUNNING: '{str}'\n"+"="*80+"\n\n")
        fp.close()
        str = f"{str} >> {logfile}"
        
    ret_val = subprocess.run(str,shell=True,
                             stdout=subprocess.PIPE,stderr=subprocess.PIPE)

    del os.environ['SIXTE_USE_PSEUDO_RNG']

    return ret_val


def phoimg(phlist,implist,xmlfile,expos=STDTEST.expos,
           ra=STDTEST.RA,dec=STDTEST.Dec,
           instrument=STDTEST.inst,
           mission=STDTEST.miss,
           mode=STDTEST.mode,
           seed=STDTEST.seed,
           mjdref=STDTEST.mjdref,
           tstart=STDTEST.tstart,
           clobber="yes",
           logfile=-1,
           test=-1):

    if (test):
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

    if (logfile!=-1):
        fp  = open(logfile,'w')
        fp.write(f" *** RUNNING: '{str}'\n"+"="*80+"\n\n")
        fp.close()
        str = f"{str} >> {logfile}"

    ret_val = subprocess.run(str,shell=True,
                             stdout=subprocess.PIPE,stderr=subprocess.PIPE)

    del os.environ['SIXTE_USE_PSEUDO_RNG']

    return ret_val


def gendetsim(implist,rawdata,xmlfile,expos=STDTEST.expos,
           instrument=STDTEST.inst,
           mission=STDTEST.miss,
           mode=STDTEST.mode,
           seed=STDTEST.seed,
           mjdref=STDTEST.mjdref,
           tstart=STDTEST.tstart,
           clobber="yes",
           background=STDTEST.background,
           logfile=-1,
           test=-1):

    if (test):
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

    if (logfile!=-1):
        fp  = open(logfile,'w')
        fp.write(f" *** RUNNING: '{str}'\n"+"="*80+"\n\n")
        fp.close()
        str = f"{str} >> {logfile}"

    ret_val = subprocess.run(str,shell=True,
                             stdout=subprocess.PIPE,stderr=subprocess.PIPE)

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
           logfile=-1,
           test=-1):

    if (test):
        os.environ["SIXTE_USE_PSEUDO_RNG"] = "1"


    pointing_string = STDTEST.get_pointing_string(ra,dec,attitude)
    impfile_string = ("" if implist==0 else f"ImpactList={implist} ")
    rawdata_string = ("" if implist==0 else f"RawData={rawdata} ")

    str = f"""runsixt \
    {pointing_string} \
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


    if (logfile!=-1):
        fp  = open(logfile,'w')
        fp.write(f" *** RUNNING: '{str}'\n"+"="*80+"\n\n")
        fp.close()
        str = f"{str} >> {logfile}"

    ret_val = subprocess.run(str,shell=True,
                             stdout=subprocess.PIPE,stderr=subprocess.PIPE)

    del os.environ['SIXTE_USE_PSEUDO_RNG']

    return ret_val



def check_returncode(ret_val,tool):
    if (ret_val.returncode==0):
        print(f"{tool}: run SUCCESSFUL")
    else:
        print('{}: run FAILED with ReturnCode {}'.format(tool,ret_val.returncode))
        out = ("{}".format(ret_val.stdout)).replace('\\n','\n')
        out = out.replace('\'','')
        out = out.replace('b','')


        err = ("{}".format(ret_val.stderr)).replace('\\n','\n')
        err = err.replace('\'','')
        err = err.replace('b','')

        print(f"{tool}: stdout: {out}")
        print(f"{tool}: stderr: {err}")
        exit(ret_val.returncode)


def check_returnvalue(ret_val,tool):
    if (ret_val==0):
        print(f"{tool}: run SUCCESSFUL")
    else:
        print('{}: run FAILED with ReturnCode {}'.format(tool,ret_val))
        exit(ret_val)


def check_fdiff(file1,file2,tool):
    output = fits.FITSDiff(file1,file2,ignore_comments=['*'],ignore_keywords=['DATE','HISTORY','CHECKSUM'])

    if (output.identical):
        print(f'{tool}: comparison to refdata SUCCESSFUL')
    elif (output.identical == False):
        print(output.report())
        print(f'*** error *** {tool}: comparison to refdata FAILED')
        exit(1)

def clean_output():
    subprocess.call('rm -f *.fits', shell=True)

