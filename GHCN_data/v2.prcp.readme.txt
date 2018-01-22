GHCN Version 2 Precipitation Version 2 Documentation

These files are associated with GHCN v2 precipitation:
-readme.precip.V2
-v2.prcp.Z
-v2.prcp_adj.Z
-v2.prcp.inv
-v2.prcp.failed.qc.Z
-v2.country.codes

New monthly data are added to v2.prcp a few days after the end of
the month.  Please note that sometimes these new data are later
replaced with data with different values due to, for example,
occasional corrections to the transmitted data that countries
will send over the Global Telecommunications System.

=>  readme.prcp.readme is this brief documentation file.

The files ending with .Z are compressed.  To uncompress the files use the command:  uncompress filename.  For more information,
email:  questions@ncdc.noaa.gov.

=>  v2.prcp is the raw data file.

The data are monthly total precipitation recorded
at the station in tenths of mm.  (Divide by 10 to get millimeters.)

Each line of the data file has:
 
station number which has three parts:
        country code (3 digits)
        nearest WMO station number (5 digits)
        modifier (3 digits) (this is usually 000 if it is that WMO station)
 
This file does not contain any duplicates of the stations, however, it
does have a duplicate number which indicates whether there are duplicates
in the companion file.  A 1 in the duplicate number indicates there
are no duplicates and a 0 indicates that there are duplicates in the
file v2.prcp.duplicates.  

Year:
        four digit year
 
Data:
        12 monthly values each as a 5 digit integer.  
	The data are monthly total precipitation recorded
	at the station in tenths of mm.  (Divide by 10 to get millimeters.)
	Missing monthly values are given as -9999.
	Trace precipitation is indicated by -8888.

=>  v2.prcp_adj is the adjusted data file.

This file contains far fewer station records than the raw data set,
It has the same format as the raw data file, but contains
data that has been adjusted for inhomogeneities.  Not only are these
station records free of inhomogeneities, they also
contain additional records for some stations, particularly those in
the former Soviet Union.

=>  v2.prcp.inv  Metadata file.
This metadata file contains the station id, station name,
country, latitude, longitude, and elevation
The format is as follows:
station number (i11), space (1x), station name (a20), country (a10),
latitude (f7.2), longitude (f8.2), and elevation in meters (i5)

Country codes:
        The file v2.country.codes lists the countries of the world and
GHCN's numerical country code.

=>  v2.prcp.failed.qc  Like GHCN temperature data, the 
precipitation data base has undergone a series of Quality Control
tests.  The final test removed individual data points we determined were most likely erroneous (e.g., due to errors in digitizing).
This final QC step is very important because some bad
values can make it into any data base.  However, any QC approach
that can detect most of the bad values will also flag some valid
extreme precipitation data points.  Therefore, we can say
with confidence, that some of the data points that failed
our QC are valid.  This file has the data points that failed our
QC and were removed from the data file is provided for users
who might have additional corroborating information such 
as an old news report of an unusual flood near an isolated
station.  Without such information, we do not recommend using
these data points.

=> v2.country.codes lists the countries of the world and
GHCN's numerical country code.

