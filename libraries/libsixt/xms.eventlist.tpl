################################################################################
XTENSION        BINTABLE        / Binary table extension
EXTNAME         EVENTS          / Extension name
ORIGIN          REMEIS          / Origin of FITS File
CREATOR         String          / Program that created this FITS file
MISSION         String          / Mission name
TELESCOP        IXO             / Telescope name
INSTRUME        XMS             / Instrument name
OBS_MODE        String          / Observation mode
DATAMODE        String          / Instrument data mode
FRAMETIM        Real            / nominal frame time
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
TIMTSYS         TT              / Time system (Terrestial Time)
RA_PNT          Real            / actual pointing RA
DEC_PNT         Real            / actual pointing DEC
RADECSYS        FK5             / Stellar reference frame
EQUINOX         2000.0          / Coordinate system equinox
LONGSTR         OGIP 1.0        / support multi-line COMMENTs or HISTORY records
#Column definitions
	TTYPE#  TIME            / Time of event
	TFORM#  D               / Format of column TIME
	TUNIT#  s               / Unit of column TIME

	TTYPE#  PHA             / Uncorrected Event Energy
	TFORM#  I               / Format of column PHA
	TUNIT#  ADU             / Unit of column PHA

	TTYPE#	RAWX            / Event X Position
	TFORM#  I               / Format of column RAWX
	TUNIT#  pixel           / Unit of column RAWX
	TTYPE#	RAWY            / Event Y Position
	TFORM#  I               / Format of column RAWY
	TUNIT#  pixel           / Unit of column RAWY

	TTYPE#	GRADE           / Energy Accuracy (1=high, 2=intermediate, ...)
	TFORM#  I               /
	TUNIT#                  /

	TTYPE#	ARRAY           / TES Array (1=inner, 2=outer)
	TFORM#  I               /
	TUNIT#                  /
