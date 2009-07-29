################################################################################
XTENSION        BINTABLE        / Binary table extension
EXTNAME         EVENTS          / Extension name
ORIGIN          REMEIS          / Origin of FITS File
CREATOR         String          / Program that created this FITS file
MISSION         REMEIS          / Mission name
TELESCOP        IXO             / Telescope name
INSTRUME        WFI             / Instrument name
DATAMODE        FULL            / Instrument data mode
SUBMODE         ReadoutDir2     / sub-mode: 2 readout directions
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

# For the pattern recognition algorithm
COLUMN          COLUMN
ROW             ROW
COLUMNS         1024
ROWS            1024

#Column definitions
	TTYPE#  TIME            / Time of event
	TFORM#  D               / Format of column TIME
	TUNIT#  s               / Unit of column TIME

	TTYPE#  PHA             / Uncorrected Event Energy
	TFORM#  J               / Format of column PHA
	TUNIT#  ADU             / Unit of column PHA

	TTYPE#	COLUMN          / Event X Position
	TFORM#  I               / Format of column RAWX
	TUNIT#  pixel           / Unit of column RAWX
	TTYPE#	ROW             / Event Y Position
	TFORM#  I               / Format of column RAWY
	TUNIT#  pixel           / Unit of column RAWY

	TTYPE#  FRAME           / Frame number of event
	TFORM#  J               / Format of column FRAME
	
	TTYPE#  PATNUM          / Number of pixels in event
	TFORM#  I		/ Format of column PATNUM

	TTYPE#  PATID           / ID with respect to XMM
	TFORM#  I               / Format of column PATID

	TTYPE#  PILEUP          / Has energy pileup occurred? Needed for patternrecog.
	TFORM#  J               / Format of column PILEUP
