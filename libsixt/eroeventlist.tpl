XTENSION        BINTABLE        / Binary table extension
EXTNAME         EVENTS          / Extension name
HDUCLASS        OGIP            /
HDUCLAS1        EVENTS          /
ORIGIN          REMEIS          / Origin of FITS File
CREATOR         FAU             / Program that created this FITS file
MISSION         SRG             / Mission name
TELESCOP        eROSITA         / Telescope name
INSTRUME        FM4             / Instrument name
OBS_MODE        String          / Observation mode
DATAMODE        String          / Instrument data mode
FILTER          OPEN            / CCD filter used
OBS_ID          String          / Observation Identifier
EXP_ID          String          / Exposure Identifier
OBSERVER        String          / Name of PI
OBJECT          String          / Name of observed Object
RA_OBJ          Real            / Source right ascension in degrees
DEC_OBJ         Real            / Source declination in degrees
DATE-OBS        UTC_format      / Date of the start of the observation
#DATE-END        UTC_format      / Date of the end of the observation
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
NXDIM           384             /
NYDIM           384             /
PIXLEN_X        75.0            /
PIXLEN_Y        75.0            /
#Column definitions
	TTYPE#  TIME            / Time of event detection
	TFORM#  D               / Format of column TIME
	TUNIT#  s               / Unit of column TIME
	TTYPE#  FRAME           / Frame counter
	TFORM#  J               / Format of column FRAME
	TTYPE#  PHA             / Uncorrected Event Energy
	TFORM#  J               / Format of column PHA
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
	TTYPE# 	CCDNR           / 1-7 for eROSITA
	TFORM# 	I   		/ Format of column CCDNR
