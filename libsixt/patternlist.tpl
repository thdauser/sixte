SIMPLE  = T
BITPIX  = 8
ORIGIN  = ECAP
CREATOR = Schmid
TELESCOP= AthenaPlus
INSTRUME= XIFU
FILTER  = Normal

XTENSION= BINTABLE        / Binary table extension
EXTNAME = EVENTS          / Extension name
ORIGIN  = ECAP            / Origin of FITS File
CREATOR = Schmid          / Program that created this FITS file
TELESCOP= String          / Telescope name
INSTRUME= String          / Instrument name
OBS_MODE= String          / Observation mode
DATAMODE= String          / Instrument data mode
FILTER  = String          / CCD filter used
MJDREF  = 54101           / Modified Julian Date of time origin
TIMEZERO= 0.0             / Time correction
TIMEUNIT= s               / Time unit
TIMESYS = TT              / Time system (Terrestial Time)
RADECSYS= FK5             / Stellar reference frame
EQUINOX = 2000.0          / Coordinate system equinox
LONGSTR = OGIP 1.0        / support multi-line COMMENTs or HISTORY records
#Column definitions
	TTYPE#  TIME            / Time of event detection
	TFORM#  D               / Format of column TIME
	TUNIT#  s               / Unit of column TIME
	TTYPE#  FRAME           / Frame counter
	TFORM#  J               / Format of column FRAME
	TUNIT#                  / Unit of column FRAME
	TTYPE#  PHA             / Uncorrected Event Energy
	TFORM#  J               / Format of column PHA
	TUNIT#  ADU             / Unit of column PHA
	TTYPE#  SIGNAL          / Pixel SIGNAL
	TFORM#  E               / Format of column SIGNAL
	TUNIT#  keV             / Unit of column SIGNAL
	TTYPE#  RAWX            / Raw x-coordinate of pixel
	TFORM#  I               / Format of column RAWX
	TUNIT#  pixel           / Unit of column RAWX
	TTYPE#  RAWY            / Raw y-coordinate of pixel
	TFORM#  I               / Format of column RAWY
	TUNIT#  pixel           / Unit of column RAWY
	TTYPE#  RA              / Back-projected right ascension
	TFORM#  D               / Format of column RA
	TUNIT#  deg             / Unit of column RA
	TTYPE#  DEC             / Back-projected declination
	TFORM#  D               / Format of column DEC
	TUNIT#  deg             / Unit of column DEC
	TTYPE# 	PH_ID           / Photon ID
	TFORM# 	2J   		/
	TUNIT# 	                / 
	TTYPE# 	SRC_ID          / Source ID
	TFORM# 	2J   		/
	TUNIT# 	                / 
	TTYPE#  TYPE            / Pattern type
	TFORM#  I               / Format of column TYPE
	TUNIT#                  / Unit of column TYPE
	TTYPE#  NPIXELS         / Number of involved pixels
	TFORM#  J               / Format of column NPIXELS
	TUNIT#                  / Unit of column NPIXELS
	TTYPE#  PILEUP          / Pile-up flag
	TFORM#  I               / Format of column PILEUP
	TUNIT#                  / Unit of column PILEUP
	TTYPE#  SIGNALS         / Signals in surrounding 3x3 matrix
	TFORM#  9E              / Format of column SIGNALS
	TUNIT#                  / Unit of column SIGNALS
