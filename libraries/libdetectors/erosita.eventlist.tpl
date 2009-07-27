################################################################################
# general keywords, common to most(all) eROSITA fits files
XTENSION        BINTABLE        / Binary table extension
EXTNAME         EVENTS          / Extension name
ORIGIN          REMEIS          / Origin of FITS File
CREATOR         String          / Program that created this FITS file
MISSION         String          / Mission name
TELESCOP        String          / Telescope name
INSTRUME        String          / Instrument name
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
# keywords specific for eROSITA event files
EXTNAME         EVENTS          / Name of FITS table extension
HDUCLASS        OGIP            / Type/source of format
HDUCLAS1        EVENTS          / data type
# Housekeeping parameter
NLLRAT          Real            / Trigger threshold
SPLTTHR         Real            / Threshold for Split events
THCODE          String          / Code name of the Threshold map
NOCODE          String          / Code name of the Noise map
OFCODE          String          / Code name of the Offset map
GAINCAME        String          / Code name of the CAMEX amplification
TIMEDEL         Real            / Length of exposure entry interval [s]
CHOPPER         String          / electronic chopper in use
#Column definitions
	TTYPE#  FRAME           / Frame number
	TFORM#  J               / Format of column FRAME
	TTYPE#  TIME            / Time of event
	TFORM#  D               / Format of column TIME
	TUNIT#  s               / Unit of column TIME

	TTYPE#  PHA             / Uncorrected Event Energy
	TFORM#  I               / Format of column PHA
	TUNIT#  ADU             / Unit of column PHA
	TTYPE#	ENERGY          / Calibrated event energy
	TFORM#  E               / Format of column RAWX
	TUNIT#  eV              / Unit of column RAWX

	TTYPE#	RAWX            / Event X Position
	TFORM#  I               / Format of column RAWX
	TUNIT#  pixel           / Unit of column RAWX
	TTYPE#	RAWY            / Event Y Position
	TFORM#  I               / Format of column RAWY
	TUNIT#  pixel           / Unit of column RAWY

	TTYPE#	RA              / right ascension (J2000)
	TFORM#  D               / Format of column RA
	TUNIT#  deg             / Unit of column RA
	TTYPE#	DEC             / declination (J2000)
	TFORM#  D               / Format of column DEC
	TUNIT#  deg             / Unit of column DEC
	
	TTYPE#  X               / sky coordinates (pixel size 0.05'')
	TFORM#  J               / Format of column X
	TUNIT#  pixel           / Unit of column X
	TTYPE#  Y               / sky coordinates (pixel size 0.05'')
	TFORM#  J               / Format of column Y
	TUNIT#  pixel           / Unit of column Y
	
	TTYPE#  PAT_TYP         / Pattern Type
	TFORM#  I		/ Format of column PAT_TYP
	TUNIT#  1234: sdtq ..   / Unit of column PAT_TYP
	TTYPE#  PAT_INF         / Pattern Information
	TFORM#  B               / Format of column PAT_INF
	TUNIT#  pattern info    / Unit of column PAT_INF
	TTYPE#  PAT_IND         / Pattern index
	TFORM#  I               / Format of column PAT_IND
	TUNIT#  pattern index   / Unit of column PAT_IND
	TTYPE#  AMP_GAIN        / Gain corrected amplitudes
	TFORM#  E               / Format of column AMP_GAIN
	TUNIT#  ADU             / Unit of column AMP_GAIN
	TTYPE#  AMP_COR         / Gain+CTI corrected amplitudes
	TFORM#  E               / Format of column AMP_GAIN
	TUNIT#  ADU             / Unit of column AMP_GAIN
