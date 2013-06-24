SIMPLE  = T
BITPIX  = 8
ORIGIN  = ECAP
CREATOR = SIXTE

XTENSION= BINTABLE        / Binary table extension
EXTNAME = EVENTS          / Extension name
HDUCLASS= OGIP            /
HDUCLAS1= EVENTS          /
ORIGIN  = ECAP            / Origin of FITS File
CREATOR = SIXTE           / Program that created this FITS file
MISSION = SRG             / Mission name
TELESCOP= eROSITA         / Telescope name
INSTRUME= FM4             / Instrument name
OBS_MODE= String          / Observation mode
DATAMODE= String          / Instrument data mode
FRAMETIM= 50.0            / [ms] nominal frame time
FILTER  = OPEN            / CCD filter used
OBS_ID  = String          / Observation Identifier
EXP_ID  = String          / Exposure Identifier
OBSERVER= String          / Name of PI
OBJECT  = NA              / Name of observed Object
RA_OBJ  = NA              / [deg] J2000
DEC_OBJ = NA              / [deg] J2000
DATE-OBS= UTC_format      / Date of the start of the observation
TSTART  = Real            / Start time of the observation
TSTOP   = Real            / End time of the observation
TEND    = Real            / End time of the observation
MJDREF  = 54101           / Modified Julian Date of time origin
TIMEZERO= 0.0             / Time correction
TIMEUNIT= s               / Time unit
TIMESYS = TT              / Time system (Terrestial Time)
TIMEDEL = 0.05            / [s] Length of exposure interval
RA_PNT  = NA              / Actual pointing RA
DEC_PNT = NA              / Actual pointing DEC
RADECSYS= FK5             / Stellar reference frame
EQUINOX = 2000.0          / Coordinate system equinox
LONGSTR = OGIP 1.0        / Support multi-line COMMENTs or HISTORY records
NXDIM   = 384             / CCD dimension
NYDIM   = 384             / CCD dimension
PIXLEN_X= 75.0            / [micron] Pixel size
PIXLEN_Y= 75.0            / [micron] Pixel size
#Column definitions
	TTYPE#  TIME            / Time of event detection
	TFORM#  D               / Format of column TIME
	TUNIT#  s               / Unit of column TIME
	TTYPE#  FRAME           / Frame counter
	TFORM#  J               / Format of column FRAME
	TTYPE#  PHA             / Uncorrected Event Energy
	TFORM#  I               / Format of column PHA
	TUNIT#  ADU             / Unit of column PHA
	TLMIN#  0
	TLMAX#  4095
	TTYPE#  ENERGY          / Recombined photon energy [eV]
	TFORM#  E               
	TUNIT#  eV              
	TLMIN#  0.0
	TLMAX#  20.48
	TTYPE#  RAWX            / Raw x-coordinate of pixel
	TFORM#  I               / Format of column RAWX
	TUNIT#  pixel           / Unit of column RAWX
	TTYPE#  RAWY            / Raw y-coordinate of pixel
	TFORM#  I               / Format of column RAWY
	TUNIT#  pixel           / Unit of column RAWY
	TTYPE#  RA              / Back-projected right ascension
	TFORM#  J               / Format of column RA
	TTYPE#  DEC             / Back-projected declination
	TFORM#  J               / Format of column DEC
	TTYPE#  X               / Back-projected right ascension
	TFORM#  J               / Format of column X
	TCTYP#  RA---SIN        / WCS Coord. type
	TCUNI#  deg             / WCS  physical unit of X axis
	TCRPX#  0.0             / WCS axis reference pixel
	TCRVL#  0.0             / [deg] WCS coord. at X axis ref. pixel
	TCDLT#  1.388889e-5     / [deg/pix] WCS X increment at ref. pixel
	TLMIN#  long            / Minimum for X
	TLMAX#  long            / Maximum for X
	TTYPE#  Y               / Back-projected declination
	TFORM#  J               / Format of column DEC
	TCTYP#  DEC--SIN        / WCS Coord. type
	TCUNI#  deg             / WCS  physical unit of Y axis 
	TCRPX#  0.0             / WCS axis reference pixel
	TCRVL#  0.0             / [deg] WCS coord. at Y axis ref. pixel
	TCDLT#  1.388889e-5     / [deg/pix] WCS Y increment at ref. pixel
	TLMIN#  long            / Minimum for Y
	TLMAX#  long            / Maximum for Y
	REFXCTYP RA---SIN
	REFXCUNI deg
	REFXCRPX 0.0
	REFXCRVL 0.0
	REFXCDLT 1.388889e-5
	REFXLMIN TODO
	REFXLMAX TODO
	REFXDMIN TODO
	REFXDMAX TODO
	REFYCTYP DEC--SIN
	REFYCUNI deg
	REFYCRPX 0.0
	REFYCRVL 0.0
	REFYCDLT 1.388889e-5
	REFYLMIN TODO
	REFYLMAX TODO
	REFYDMIN TODO
	REFYDMAX TODO
	TTYPE# 	CCDNR           / 1-7 for eROSITA
	TFORM# 	I   		/ Format of column CCDNR
	TTYPE#  FLAG            / flag
	TFORM#  J               / Format of column FLAG
	TTYPE#  PAT_TYP         / pattern type
	TFORM#  I               / Format of column PAT_TYP
	TTYPE#  PAT_INF         / pattern information
	TFORM#  U               / Format of column PAT_INF
	TTYPE#  EV_WEIGHT       / Inverse vignetting correction factor
	TFORM#  E               / Format of column EV_WEIGHT
