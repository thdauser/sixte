XTENSION        BINTABLE        / Binary table extension
EXTNAME         EVENTS          / Extension name
ORIGIN          REMEIS          / Origin of FITS File
CREATOR         FAU             / Program that created this FITS file
MISSION         String          / Mission name
TELESCOP        String          / Telescope name
INSTRUME        String          / Instrument name
OBS_MODE        String          / Observation mode
DATAMODE        String          / Instrument data mode
FILTER          String          / CCD filter used
OBS_ID          String          / Observation Identifier
EXP_ID          String          / Exposure Identifier
OBSERVER        String          / Name of PI
OBJECT          String          / Name of observed Object
RA_OBJ          Real            / Source right ascension in degrees
DEC_OBJ         Real            / Source declination in degrees
DATE-OBS        UTC_format      / Date of the start of the observation
DATE-END        UTC_format      / Date of the end of the observation
TSTART          Real            / Start time of the observation
TSTOP           Real            / End time of the observation
TEND            Real            / End time of the observation
MJDREF          54101           / Modified Julian Date of time origin
TIMEZERO        0.0             / Time correction
TIMEUNIT        s               / Time unit
TIMESYS         TT              / Time system (Terrestial Time)
RA_PNT          Real            / actual pointing RA
DEC_PNT         Real            / actual pointing DEC
RADECSYS        FK5             / Stellar reference frame
EQUINOX         2000.0          / Coordinate system equinox
LONGSTR         OGIP 1.0        / support multi-line COMMENTs or HISTORY records
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
	TTYPE#  CHARGE          / Pixel charge
	TFORM#  E               / Format of column CHARGE
	TUNIT#  keV             / Unit of column CHARGE
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
