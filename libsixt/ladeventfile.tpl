SIMPLE  = T
BITPIX  = 8
ORIGIN  = ECAP
CREATOR = ladsim

XTENSION= BINTABLE        / Binary table extension
EXTNAME = EVENTS          / Extension name
ORIGIN  = ECAP            / Origin of FITS File
CREATOR = ladsim          / Program that created this FITS file
OBS_MODE= POINTED         / Observation mode
DATAMODE= EVENT           / Instrument data mode
MJDREF  = 54101           / Modified Julian Date of time origin
TIMEZERO= 0.0             / Time correction
TIMEUNIT= s               / Time unit
TIMESYS = TT              / Time system (Terrestial Time)
LONGSTR = OGIP 1.0        / support multi-line COMMENTs or HISTORY records
#Column definitions
	TTYPE#  TIME            / Time of event detection
	TFORM#  D               / Format of column TIME
	TUNIT#  s               / Unit of column TIME
	TTYPE#  FRAME           / Frame counter
	TFORM#  J               / Format of column FRAME
	TTYPE#  SIGNAL          / Measured signal
	TFORM#  E               / 
	TUNIT#  keV             / 
	TTYPE#  PANEL           / Affected panel
	TFORM#  J   		/
	TTYPE#  MODULE          / Affected module
	TFORM#  J   		/
	TTYPE#  ELEMENT         / Affected element
	TFORM#  J   		/
	TTYPE#  ANODE           / Affected anode
	TFORM#  J   		/
	TTYPE# 	PH_ID           / Photon ID
	TFORM# 	2J   		/
	TTYPE# 	SRC_ID          / Source ID
	TFORM# 	2J   		/
