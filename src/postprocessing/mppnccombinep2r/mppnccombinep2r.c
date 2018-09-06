/*
   mppnccombine - joins together netCDF data files representing a decomposed
   domain into a unified netCDF file.  It was originally
   designed to be used as a postprocessor for the parallel I/O
   programming interface "mpp_io_mod"

   AMFI version: Prajeesh A Gopinathan
	             Added option for combining and unpacking the dicomposed P-grid
			     files of AMFI model to regular lat-lon grid.
				 Added option for netcdf4 compression.

   V2.1.7:  Added option to initialize output variables with a missing_value
   from the variables of the first input file as suggested by
   Martin Schmidt (martin.schmidt@io-warnemuende.de) and
   Franz Tauber (franz.tauber@io-warnemuende.de).
   V2.1.6:  Bug fixes for greater than 4GB record sizes.  Does not contain
   V2.1.5 modifications which were a special case.
   V2.1.5:  Supports running in an MPI environment.  Reads arguments from a
   configuration file instead of from the command line; this is needed
   to work around a bug in Cray's aprun.
   V2.1.4:  Fixed a bug with file cleanup and the debugging messages.
   V2.1.3:  Fixed a run-time segmentation fault with some compilers; changed
   ncinfile allocation in process_file function.
   V2.1.2:  Fixed a bug with input files that have decomposed dimensions
   defined after the variables that use them.
   V2.1.1:  Added option (-64) for creating output files with 64-bit offset
   format requiring NetCDF 3.6.x.
   V2.1: Added an option (-h) to pad the output file's header with empty space.
   Added an option (-e #) to specify an ending number to a range of input
   filename extensions.  It no longer aborts on missing input files, but
   gives error messages at the end of all the processing.
   V2.0: Substantial rewrite; memory buffering increases speed several times.
   V1.2: Added support for specifying the start number in filename extensions.
   V1.1.1: Added a fix for dimensions that are not also variables.
   V1.1: Changed loop order for increased I/O efficiency; records are now the
   innermost loop then the variables loop.
   V1.0: Original release.

   Written by Hans Vahlenkamp (Hans.Vahlenkamp)
   Geophysical Fluid Dynamics Laboratory / NOAA
   Princeton Forrestal Campus
   Last Updated: 05/15/08
   */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <netcdf.h>
#include <time.h>
#include <utime.h>


char xgrid[1024]="p_xgrd.nc";
unsigned char xgrid_init=0, unpack=0;
int nlon=0, nlat=0, xn=0;
int *pxi, *pxj, *rxi, *rxj;
double *lonc, *latc, *xf, *rxf;

unsigned char n4=0;

#define ERR(e) {if(e){printf("Error: %s\n", nc_strerror(e)); return 2;}}

/* Information structure for a file */
struct fileinfo
{
	int ncfid;  /* ID of the input netCDF file */
	int ndims;  /* Number of dimensions */
	int nvars;  /* Number of variables */
	int ngatts;  /* Number of global attributes */
	int recdim;  /* ID of the record dimensions */
	char varname[MAX_NC_VARS][MAX_NC_NAME];  /* Names of the variables */
	nc_type datatype[MAX_NC_VARS]; /* Data types of the variables */
	int varndims[MAX_NC_VARS];  /* Number of dimensions for each variable */
	int vardim[MAX_NC_VARS][MAX_NC_DIMS];  /* Dimensions for each variable */
	int natts[MAX_NC_VARS];  /* Number of attributes for each variable */
	unsigned char vardecomp[MAX_NC_VARS];  /* Is the variable decomposed */
	unsigned char varislat[MAX_NC_VARS];
	unsigned char varislon[MAX_NC_VARS];
	unsigned char islev[MAX_NC_DIMS];
	char dimname[MAX_NC_DIMS][MAX_NC_NAME];  /* Names of the dimensions */
	long dimsize[MAX_NC_DIMS];  /* Sizes of the dimensions (decomposed) */
	long dimfullsize[MAX_NC_DIMS];  /* Full sizes of the dimensions */
	long dimregsize[MAX_NC_DIMS];  /* Full sizes of the dimensions */
	long dimstart[MAX_NC_DIMS];  /* Start positions within full dimensions */
	long dimend[MAX_NC_DIMS];  /* End positions within full dimensions */
	unsigned char varmiss[MAX_NC_VARS];  /* Does variable have missing_value */
	unsigned char varmissval[MAX_NC_VARS][8];  /* missing_value per variable */
};

/* Auxiliary function prototypes */
void usage();
int process_file(char *, unsigned char, struct fileinfo *, char *, int *,
		int *, int, int, void *[], void *[], int, unsigned char,
		unsigned char);
int process_vars(struct fileinfo *, struct fileinfo *, unsigned char, int *,
		int, int, int, void *[], void *[], unsigned char, unsigned char);
int flush_decomp(struct fileinfo *, int, int, void *[], void *[], unsigned char);
void print_debug(struct fileinfo *, unsigned char);
char *nc_type_to_str(nc_type);
int read_xgrid();
int nccp2r(int, char * []);
double modtimediff(char * []);

#ifndef lib_mppnccp2r
int main(int argc, char *argv[])
{
	return(nccp2r(argc, argv));
}
#endif


int nccp2r(int argmntc, char *argmntv[])
{
	unsigned char verbose=0;  /* Print some progress information? */
	unsigned char appendnc=0;  /* Append to an existing netCDF file? */
	unsigned char removein=0;  /* Remove the ".####" decomposed input files? */
	int nstart=0;  /* PE number of the first input netCDF file */
	int nend=(-1);  /* PE number of the last input netCDF file */
	int headerpad=16384;  /* Additional padding at the end of the header */
	int format=NC_NOCLOBBER;  /* Format of new netCDF output file */
	unsigned char missing=0;  /* Initialize output variables with */
	/* "missing_value" instead of 0 value? */
	int outputarg=(-1);  /* Argument # of the output netCDF file */
	int inputarg=(-1);  /* Argument # of first input netCDF file */
	struct stat statbuf;  /* Dummy structure for file-testing "stat" call */
	struct fileinfo *ncoutfile;  /* Information for the output file */
	char outfilename[2048], *strptr;  /* Name of the output netCDF file */
	int outlen;  /* Length of the output filename */
	char infilename[2048];  /* Name of an input file */
	unsigned char infileerror=0;  /* Errors reading an input file */
	unsigned char infileerrors=0;  /* Errors reading any input files */
	int nfiles=(-1);  /* Number of files in the decomposed domain */
	int a, f, r;  /* Loop variables */
	int status; /* Return status */
	int nrecs=1;  /* Number of records in each decomposed file */
	void *varbuf[NC_MAX_VARS];  /* Buffers for decomposed variables */
	void *vbuf[NC_MAX_VARS];  /* Buffers for upacked variables */
	size_t blksz=65536; /* netCDF block size */

	verbose=0; appendnc=0; removein=0;  
	nstart=0;  nend=(-1);  headerpad=16384;  
	format=NC_NOCLOBBER;  missing=0;  outputarg=(-1);  
	inputarg=(-1);  infileerror=0;  infileerrors=0;  
	nfiles=(-1); nrecs=1; blksz=65536; 
	n4=0;
	
	if (argmntc < 2)
	{
		usage(); return(1);
	}
	for (a=1; a < argmntc; a++)
	{
		if (!strcmp(argmntv[a],"-v")) verbose=1;
		else if (!strcmp(argmntv[a],"-vv")) verbose=2;  /* Hidden debug mode */
		else if (!strcmp(argmntv[a],"-a")) appendnc=1;
		else if (!strcmp(argmntv[a],"-r")) removein=1;
		else if (!strcmp(argmntv[a],"-n"))
		{
			a++;
			if (a < argmntc) nstart=atoi(argmntv[a]);
			else
			{
				usage(); return(1);
			}
		}
		else if (!strcmp(argmntv[a],"-e"))
		{
			a++;
			if (a < argmntc) nend=atoi(argmntv[a]);
			else
			{
				usage(); return(1);
			}
		}
		else if (!strcmp(argmntv[a],"-h"))
		{
			a++;
			if (a < argmntc) headerpad=atoi(argmntv[a]);
			else
			{
				usage(); return(1);
			}
		}
		else if (!strcmp(argmntv[a],"-64"))
			format=(NC_NOCLOBBER | NC_64BIT_OFFSET);
		else if (!strcmp(argmntv[a], "-n4")) 
		{
			format=(NC_NOCLOBBER | NC_NETCDF4 | NC_CLASSIC_MODEL);
			a++;
			if (a < argmntc) n4 = atoi(argmntv[a]);
			else
			{
				usage(); return(1);
			}
			if (n4>9) n4=9;
			if (n4<0) n4=0;
		}
		else if (!strcmp(argmntv[a], "-u")) 
		{
			unpack=1;
			a++;
			if (a < argmntc) strcpy(xgrid,argmntv[a]);
			else
			{
				usage(); return(1);
			}
		}
		else if (!strcmp(argmntv[a],"-m")) missing=1;
		else
		{
			outputarg=a; break;
		}
	}
	if (outputarg==(-1))
	{
		usage(); return(1);
	}
	if (argmntc-1 > outputarg) inputarg=outputarg+1;
	sprintf(outfilename,argmntv[outputarg]); outlen=strlen(outfilename);
	if (outlen > 4)
	{
		strptr=outfilename+outlen-5;
		if (!strcmp(strptr,".0000")) outfilename[outlen-5]='\0';
	}

	printf(" combining... %s\n",outfilename);

	/* Disable fatal returns from netCDF library functions */
	ncopts=0;

	/* Read xgrid */
	if (unpack) {
		if (!xgrid_init){
			if (!read_xgrid()) {printf("Succesfully read xgrid...\n");}
			else {fprintf(stderr,"Error reading xgrid !\n"); return(1);}
			xgrid_init=1;
		}
	}

	/* Create a new netCDF output file */
	if ((ncoutfile=(struct fileinfo *)malloc(sizeof(struct fileinfo)))==NULL)
	{
		fprintf(stderr,"Error: cannot allocate enough memory!\n"); return(1);
	}
	if (!appendnc)
	{
		if (stat(outfilename,&statbuf)==0)
		{
			fprintf(stderr,"Error: output file seems to exist already!\n");
			free(ncoutfile); return(1);
		}
		status = nc__create(outfilename, format, 0, &blksz, &ncoutfile->ncfid);
		if (status==(-1))
		{
			fprintf(stderr,"Error: cannot create the output netCDF file!\n");
			free(ncoutfile); return(1);
		}
		ncsetfill(ncoutfile->ncfid,NC_NOFILL);
	}
	/* Open an existing netCDF file for appending */
	else
	{
		if ((ncoutfile->ncfid=ncopen(outfilename,NC_WRITE))==(-1))
		{
			fprintf(stderr,"Error: cannot open the output netCDF file for appending!\n");
			free(ncoutfile); return(1);
		}
	}

	for (f=0; f < NC_MAX_VARS; f++){
		varbuf[f]=NULL;
		vbuf[f]=NULL;
	}

	/* No input files are specified on the command-line */
	if (inputarg==(-1))
	{
		if (nend > -1)
			for (r=0; r < nrecs; r++)
			{
				if (verbose) printf("record = %d\n",r);
				f=0; 
				for (a=nstart; a <= nend; a++)
				{
					sprintf(infilename,"%s.%04d",outfilename,a);
					if (verbose)
					{
						if (a==nstart && r==0) printf("  n files to go... ");
						else printf("  %d files to go... ",nend-nstart+1-f);
						printf("processing \"%s\"\n",infilename);
					}
					if (stat(infilename,&statbuf)!=0) continue;
					infileerror=process_file(infilename,appendnc,ncoutfile,
							outfilename,&nfiles,&nrecs,r,f,varbuf,vbuf,
							headerpad,verbose,missing);
					if (infileerror) infileerrors=1;
					appendnc=1; f++;
					if (f==nfiles || a==nend)
					{
						if (verbose > 1)
							printf("  Write variables from previous %d files\n",f);
						if (flush_decomp(ncoutfile,nfiles,r,varbuf,vbuf,verbose)!=0) {
							fprintf(stderr,"Error while writing to disk.."); return(1);}
						break;
					}
				}
			}
		else
		{
			nend=nstart+1;
			for (r=0; r < nrecs; r++)
			{
				if (verbose) printf("record = %d\n",r);
				f=0;
				for (a=nstart; a < nend; a++)
				{
					sprintf(infilename,"%s.%04d",outfilename,a);
					if (verbose)
					{
						if (a==nstart && r==0) printf("  n files to go... ");
						else printf("  %d files to go... ",nend-a);
						printf("processing \"%s\"\n",infilename);
					}
					infileerror=process_file(infilename,appendnc,ncoutfile,
							outfilename,&nfiles,&nrecs,r,f,varbuf,vbuf,
							headerpad,verbose,missing);
					if (infileerror) infileerrors=1;
					if (a==nstart && nfiles > 0) nend=nstart+nfiles;
					appendnc=1; f++;
					if (f==nfiles || a==(nend-1))
					{
						if (verbose > 1)
							printf("  Write variables from previous %d files\n",f);
						if (flush_decomp(ncoutfile,nfiles,r,varbuf,vbuf,verbose)!=0) {
							fprintf(stderr,"Error while writing to disk.."); return(1);}
						f=0; continue;
					}
				}
			}
		}
	}
	/* Loop over all the specified input files */
	else
		for (r=0; r < nrecs; r++)
		{
			if (verbose) printf("record = %d\n",r);
			f=0;
			for (a=inputarg; a < argmntc; a++)
			{
				if (verbose)
				{
					if ((argmntc-a)==1) printf("  1 file to go... ");
					else printf("  %d files to go... ",argmntc-a);
					printf("processing \"%s\"\n",argmntv[a]);
				}
				infileerror=process_file(argmntv[a],appendnc,ncoutfile,
						outfilename,&nfiles,&nrecs,r,f,varbuf,vbuf,
						headerpad,verbose,missing);
				if (infileerror) infileerrors=1;
				appendnc=1; f++;
				if (f==nfiles || a==(argmntc-1))
				{
					if (verbose > 1)
						printf("  Write variables from previous %d files\n",f);
					if (flush_decomp(ncoutfile,nfiles,r,varbuf,vbuf,verbose)!=0) {
						fprintf(stderr,"Error while writing to disk.."); return(1);}
					f=0; continue;
				}
			}
		}

	/* Clean up... return 1 on error, otherwise 0 */
	for (f=0; f < NC_MAX_VARS; f++)
	{
		if (varbuf[f]!=NULL){ 
			free(varbuf[f]); 
			if (verbose) printf("releasing memory of %d varbuf\n",f);
		}
		if (vbuf[f]!=NULL) { 
			free(vbuf[f]); 
			if (verbose) printf("releasing memory of %d vbuf\n",f);
		}
	}
	if (verbose) printf("released memory of ncoutfile\n");
	ncclose(ncoutfile->ncfid); free(ncoutfile);

	if (!infileerrors)
	{
		if (removein)
		{
			/* No input files are specified on the command-line */
			if (inputarg==(-1))
			{
				f=0;
				for (a=nstart; a <= nend; a++)
				{
					if (++f > nfiles) break;
					sprintf(infilename,"%s.%04d",outfilename,a);
					if (verbose) printf("Removing \"%s\"\n",infilename);
					unlink(infilename);
				}
			}
			/* Loop over all the specified input files */
			else
				for (a=inputarg; a < argmntc; a++)
				{
					if (verbose) printf("Removing \"%s\"\n",argmntv[a]);
					unlink(argmntv[a]);
				}
		}
	}
	else
		fprintf(stderr,"Warning: output file may be incomplete!\n");
	return(infileerrors);
}

int read_xgrid()
{	
	int ncfid, dimid, varid;
	char dimnm[MAX_NC_NAME];
	size_t dimsize;
	long start[1], count[1];  /* Data array sizes */
	int maxrxi, maxrxj, maxpxi, maxpxj, n;

	if ((ncfid=ncopen(xgrid,NC_NOWRITE))==(-1)) 
	{ fprintf(stderr,"Error: cannot open xgrid file!\n"); return(1); } 
	if ((dimid=ncdimid(ncfid,"xn"))==(-1)) 
	{ 
		fprintf(stderr,"Error: cannot read xgrid dimension xn!\n");
		ncclose(ncfid); return(1);
	}
	if ((ncdiminq(ncfid,dimid,dimnm,&(dimsize)))==(-1))
	{ 
		fprintf(stderr,"Error: cannot read xgrid dimension #%d's metadata!\n",dimid);
		ncclose(ncfid); return(1);
	}
	xn=dimsize;	

	if ((dimid=ncdimid(ncfid,"lon"))==(-1)) 
	{ 
		fprintf(stderr,"Error: cannot read xgrid dimension lon!\n");
		ncclose(ncfid); return(1);
	}
	if ((ncdiminq(ncfid,dimid,dimnm,&(dimsize)))==(-1))
	{ 
		fprintf(stderr,"Error: cannot read xgrid dimension #%d's metadata!\n",dimid);
		ncclose(ncfid); return(1);
	}
	nlon=dimsize;	

	if ((dimid=ncdimid(ncfid,"lat"))==(-1)) 
	{ 
		fprintf(stderr,"Error: cannot read xgrid dimension lat!\n");
		ncclose(ncfid); return(1);
	}
	if ((ncdiminq(ncfid,dimid,dimnm,&(dimsize)))==(-1))
	{ 
		fprintf(stderr,"Error: cannot read xgrid dimension #%d's metadata!\n",dimid);
		ncclose(ncfid); return(1);
	}
	nlat=dimsize;	

	printf("no of exchange grids = %d\n",xn);
	printf("nlon of regular grid = %d\n",nlon);
	printf("nlat of regular grid = %d\n",nlat);

	if ((rxi = malloc(xn * sizeof(int)))==NULL) 
	{ fprintf(stderr,"Error: cannot allocate enough memory for rxi!\n"); ncclose(ncfid); return(1); }

	if ((rxj = malloc(xn * sizeof(int)))==NULL) 
	{ fprintf(stderr,"Error: cannot allocate enough memory for rxj!\n"); ncclose(ncfid); return(1); }

	if ((pxi = malloc(xn * sizeof(int)))==NULL) 
	{ fprintf(stderr,"Error: cannot allocate enough memory for pxi!\n"); ncclose(ncfid); return(1); }

	if ((pxj = malloc(xn * sizeof(int)))==NULL) 
	{ fprintf(stderr,"Error: cannot allocate enough memory for pxj!\n"); ncclose(ncfid); return(1); }

	if ((xf = malloc(xn * sizeof(double)))==NULL) 
	{ fprintf(stderr,"Error: cannot allocate enough memory for xf!\n"); ncclose(ncfid); return(1); }

	if ((rxf = malloc(xn * sizeof(double)))==NULL) 
	{ fprintf(stderr,"Error: cannot allocate enough memory for rxf!\n"); ncclose(ncfid); return(1); }

	if ((lonc = malloc(nlon * sizeof(double)))==NULL) 
	{ fprintf(stderr,"Error: cannot allocate enough memory for lonc!\n"); ncclose(ncfid); return(1); }

	if ((latc = malloc(nlat * sizeof(double)))==NULL)
	{ fprintf(stderr,"Error: cannot allocate enough memory for latc!\n"); ncclose(ncfid); return(1); }

	if ((varid=ncvarid(ncfid,"rxi"))==(-1))	
	{ fprintf(stderr,"Error: cannot get varid for rxi!\n"); ncclose(ncfid); return(1); }

	start[0] = 0; count[0] = xn;
	if (ncvarget(ncfid,varid,start,count,rxi)==(-1)) 
	{ fprintf(stderr,"Error: cannot read variable \"%s\"'s values!\n", "rxi"); ncclose(ncfid); return(1); }

	if ((varid=ncvarid(ncfid,"rxj"))==(-1))	
	{ fprintf(stderr,"Error: cannot get varid for rxj!\n"); ncclose(ncfid); return(1); }

	start[0] = 0; count[0] = xn;
	if (ncvarget(ncfid,varid,start,count,rxj)==(-1)) 
	{ fprintf(stderr,"Error: cannot read variable \"%s\"'s values!\n", "rxj"); ncclose(ncfid); return(1); }

	if ((varid=ncvarid(ncfid,"pxi"))==(-1))	
	{ fprintf(stderr,"Error: cannot get varid for pxi!\n"); ncclose(ncfid); return(1); }

	start[0] = 0; count[0] = xn;
	if (ncvarget(ncfid,varid,start,count,pxi)==(-1)) 
	{ fprintf(stderr,"Error: cannot read variable \"%s\"'s values!\n", "pxi"); ncclose(ncfid); return(1); }

	if ((varid=ncvarid(ncfid,"pxj"))==(-1))	
	{ fprintf(stderr,"Error: cannot get varid for pxj!\n"); ncclose(ncfid); return(1); }

	start[0] = 0; count[0] = xn;
	if (ncvarget(ncfid,varid,start,count,pxj)==(-1)) 
	{ fprintf(stderr,"Error: cannot read variable \"%s\"'s values!\n", "pxj"); ncclose(ncfid); return(1); }

	if ((varid=ncvarid(ncfid,"xf"))==(-1))	
	{ fprintf(stderr,"Error: cannot get varid for xf!\n"); ncclose(ncfid); return(1); }

	start[0] = 0; count[0] = xn;
	if (ncvarget(ncfid,varid,start,count,xf)==(-1)) 
	{ fprintf(stderr,"Error: cannot read variable \"%s\"'s values!\n", "xf"); ncclose(ncfid); return(1); }

	if ((varid=ncvarid(ncfid,"rxf"))==(-1))	
	{ fprintf(stderr,"Error: cannot get varid for rxf!\n"); ncclose(ncfid); return(1); }

	start[0] = 0; count[0] = xn;
	if (ncvarget(ncfid,varid,start,count,rxf)==(-1)) 
	{ fprintf(stderr,"Error: cannot read variable \"%s\"'s values!\n", "rxf"); ncclose(ncfid); return(1); }

	if ((varid=ncvarid(ncfid,"lon"))==(-1))	
	{ fprintf(stderr,"Error: cannot get varid for lon!\n"); ncclose(ncfid); return(1); }

	start[0] = 0; count[0] = nlon;
	if (ncvarget(ncfid,varid,start,count,lonc)==(-1)) 
	{ fprintf(stderr,"Error: cannot read variable \"%s\"'s values!\n", "lon"); ncclose(ncfid); return(1); }

	if ((varid=ncvarid(ncfid,"lat"))==(-1))	
	{ fprintf(stderr,"Error: cannot get varid for lat!\n"); ncclose(ncfid); return(1); }

	start[0] = 0; count[0] = nlat;
	if (ncvarget(ncfid,varid,start,count,latc)==(-1)) 
	{ fprintf(stderr,"Error: cannot read variable \"%s\"'s values!\n", "lat"); ncclose(ncfid); return(1); }

	maxrxi=maxrxj=maxpxi=maxpxj=0;
	for (n=0; n < xn; n++) {
		if (maxrxi<rxi[n]) maxrxi=rxi[n];
		if (maxrxj<rxj[n]) maxrxj=rxj[n];
		if (maxpxj<pxj[n]) maxpxj=pxj[n];
		if (maxpxi<pxi[n]) maxpxi=pxi[n];
		// fortran indexing to c indexing
		rxi[n]=rxi[n]-1;
		pxi[n]=pxi[n]-1;
		pxj[n]=pxj[n]-1;
		rxj[n]=rxj[n]-1;
	}

	printf("max value of rxi = %d, rxj = %d, pxi = %d, pxj = %d \n", 
			maxrxi, maxrxj, maxpxi, maxpxj);

	return(0);
}


/* Print the usage message for mppnccombine */
void usage()
{
	printf("mppnccombine 2.1.7 - (written by Hans.Vahlenkamp)\n\n");
	printf("Usage:  mppnccombine [-v] [-a] [-r] [-n #] [-e #] [-h #] [-64] [-n4 #] [-m]\n");
	printf("                     output.nc [input ...]\n\n");
	printf("  -v    Print some progress information.\n");
	printf("  -a    Append to an existing netCDF file (not heavily tested...).\n");
	printf("  -r    Remove the \".####\" decomposed files after a successful run.\n");
	printf("  -n #  Input filename extensions start with number #### instead of 0000.\n");
	printf("  -e #  Ending number #### of a specified range of input filename extensions.\n");
	printf("        Files within the range do not have to be consecutively numbered.\n");
	printf("  -h #  Add a specified number of bytes of padding at the end of the header.\n");
	printf("  -64   Create netCDF output files with the 64-bit offset format.\n");
	printf("  -n4 # Create netCDF output files in NETCDF4_CLASSIC mode, # deflate level. \n");
	printf("  -u #  Unpack AMFI P-grid to latlon, # is the path to p_xgrd.nc file. \n");
	printf("  -m    Initialize output variables with a \"missing_value\" from the variables\n");
	printf("        of the first input file instead of the default 0 value.\n\n");
	printf("mppnccombine joins together an arbitrary number of netCDF input files, each\n");
	printf("containing parts of a decomposed domain, into a unified netCDF output file.\n");
	printf("An output file must be specified and it is assumed to be the first filename\n");
	printf("argument.  If the output file already exists, then it will not be modified\n");
	printf("unless the option is chosen to append to it.  If no input files are specified\n");
	printf("then their names will be based on the name of the output file plus the default\n");
	printf("numeric extension \".0000\", which will increment by 1.  There is an option for\n");
	printf("starting the filename extensions with an arbitrary number instead of 0.  There\n");
	printf("is an option for specifying an end to the range of filename extension numbers;\n");
	printf("files within the range do not have to be consecutively numbered.  If input\n");
	printf("files are specified then names will be used verbatim.\n\n");
	printf("A value of 0 is returned if execution completed successfully; a value of 1\n");
	printf("otherwise.\n");
}


/* Open an input file and get some information about it, define the   */
/* structure of the output file if necessary, prepare to copy all the */
/* variables at the current record to memory                          */
int process_file(char *ncname, unsigned char appendnc,
		struct fileinfo *ncoutfile, char *outncname, int *nfiles,
		int *nrecs, int r, int f, void *varbuf[], void *vbuf[], int headerpad,
		unsigned char verbose, unsigned char missing)
{
	struct fileinfo *ncinfile;  /* Information about an input netCDF file */
	int nfiles2;  /* Number of files in the decomposed domain */
	int d, v, n;  /* Loop variables */
	int dimid;  /* ID of a dimension */
	int decomp[4];  /* "domain_decomposition = #0, #1, #2, #3" attribute */
	/*  #0 starting position of original dimension   */
	/*  #1 ending position of original dimension     */
	/*  #2 starting position of decomposed dimension */
	/*  #3 ending position of decomposed dimension   */
	char attname[MAX_NC_NAME];  /* Name of a global or variable attribute */
	unsigned char ncinfileerror=0;  /* Were there any file errors? */
	char cart[]="N";

	strcpy(cart,"N");

	/* Information for netCDF input file */
	if ((ncinfile=(struct fileinfo *)malloc(sizeof(struct fileinfo)))==NULL)
	{
		fprintf(stderr,"Error: cannot allocate enough memory!\n"); return(1);
	}

	/* Open an input netCDF file */
	if ((ncinfile->ncfid=ncopen(ncname,NC_NOWRITE))==(-1))
	{
		fprintf(stderr,"Error: cannot open input file \"%s\"\n",ncname);
		free(ncinfile); return(1);
	}

	/* Determine the number of files in the decomposed domain */
	if (ncattget(ncinfile->ncfid,NC_GLOBAL,"NumFilesInSet",
				(void *)&nfiles2)==(-1))
	{
		if (*nfiles==1)
		{
			fprintf(stderr,"Error: missing the \"NumFilesInSet\" global attribute!\n");
			return(1);
		}
		else if (*nfiles==(-1))
		{
			fprintf(stderr,"Warning: missing the \"NumFilesInSet\" global attribute.\n");
		}
	}
	*nfiles=nfiles2;

	/* Get some general information about the input netCDF file */
	if (ncinquire(ncinfile->ncfid,&(ncinfile->ndims),&(ncinfile->nvars),
				&(ncinfile->ngatts),&(ncinfile->recdim))==(-1))
	{
		fprintf(stderr,"Error: cannot read the file's metadata!\n");
		ncclose(ncinfile->ncfid); free(ncinfile); return(1);
	}

	/* Get some information about the dimensions */
	for (d=0; d < ncinfile->ndims; d++)
	{
		if ((ncdiminq(ncinfile->ncfid,d,ncinfile->dimname[d],
						&(ncinfile->dimsize[d])))==(-1))
		{
			fprintf(stderr,"Error: cannot read dimension #%d's metadata!\n",d);
			ncclose(ncinfile->ncfid); free(ncinfile); return(1);
		}
		ncinfile->dimfullsize[d]=ncinfile->dimsize[d];
		ncinfile->dimregsize[d]=ncinfile->dimsize[d];
		ncinfile->dimstart[d]=1; ncinfile->dimend[d]=(-1);
	}

	/* Save some information for the output file */
	if (r==0)
	{
		ncoutfile->nvars=ncinfile->nvars; ncoutfile->recdim=ncinfile->recdim;
	}

	/* Get some information about the variables */
	for (v=0; v < ncinfile->nvars; v++)
	{
		if ((ncvarinq(ncinfile->ncfid,v,ncinfile->varname[v],
						&(ncinfile->datatype[v]),&(ncinfile->varndims[v]),
						ncinfile->vardim[v],&(ncinfile->natts[v])))==(-1))
		{
			fprintf(stderr,"Error: cannot read variable #%d's metadata!\n",v);
			ncclose(ncinfile->ncfid); free(ncinfile); return(1);
		}

		/* If the variable is also a dimension then get decomposition info */
		ncinfile->varislat[v]=ncinfile->varislon[v]=0;
		if ((dimid=ncdimid(ncinfile->ncfid,ncinfile->varname[v]))!=(-1))
		{
			if (ncattget(ncinfile->ncfid,v,"domain_decomposition",
						(void *)decomp)!=(-1))
			{
				ncinfile->dimfullsize[dimid]=decomp[1]-decomp[0]+1;
				ncinfile->dimstart[dimid]=decomp[2]-(decomp[0]-1);
				ncinfile->dimend[dimid]=decomp[3]-(decomp[0]-1);
			}
			else
			{
				ncinfile->dimfullsize[dimid]=ncinfile->dimsize[dimid];
				ncinfile->dimregsize[dimid]=ncinfile->dimsize[dimid];
				ncinfile->dimstart[dimid]=1; ncinfile->dimend[dimid]=(-1);
			}

			if (ncattget(ncinfile->ncfid,v,"cartesian_axis", cart)!=(-1))
			{   
				if(verbose>1) printf("%s %s\n", ncinfile->varname[v], cart);
				if (!strcmp(cart,"Z")) { ncinfile->islev[dimid]=1; }
				else { ncinfile->islev[dimid]=0; }
				if (!strcmp(cart,"X")) { 
					if(verbose>1) printf("%s %s\n", ncinfile->varname[v], cart);
					ncinfile->dimregsize[dimid]=nlon; 
					ncinfile->varislon[v]=1;
				}
				if (!strcmp(cart,"Y")) {
					if(verbose>1) printf("%s %s\n", ncinfile->varname[v], cart);
					ncinfile->dimregsize[dimid]=nlat; 
					ncinfile->varislat[v]=1;
				}
			}
			else
			{
				ncinfile->islev[dimid]=0;
			}
		}
	}

	/* Get some additional information about the variables */
	for (v=0; v < ncinfile->nvars; v++)
	{
		/* Does the variable have a decomposed dimension? */
		ncinfile->vardecomp[v]=0;
		for (d=0; d < ncinfile->varndims[v]; d++)
		{
			if (ncinfile->dimend[ncinfile->vardim[v][d]]!=(-1))
			{
				ncinfile->vardecomp[v]=1; 
			}
		}

		/* Save some information for the output file */
		if (r==0)
		{
			ncoutfile->varndims[v]=ncinfile->varndims[v];
			for (d=0; d < ncinfile->ndims; d++)
				ncoutfile->dimfullsize[d]=ncinfile->dimfullsize[d];
			for (d=0; d < ncinfile->ndims; d++)
				ncoutfile->dimregsize[d]=ncinfile->dimregsize[d];
			for (d=0; d < ncinfile->varndims[v]; d++)
				ncoutfile->vardim[v][d]=ncinfile->vardim[v][d];
			ncoutfile->vardecomp[v]=ncinfile->vardecomp[v];
			ncoutfile->varislat[v]=ncinfile->varislat[v];
			ncoutfile->varislon[v]=ncinfile->varislon[v];
			ncoutfile->datatype[v]=ncinfile->datatype[v];
			strcpy(ncoutfile->varname[v],ncinfile->varname[v]);
			ncoutfile->varmiss[v]=0;
		}
	}

	/* If the output netCDF file was just created then define its structure */
	if (!appendnc)
	{
		if (verbose) printf("    Creating output \"%s\"\n",outncname);

		/* Define the dimensions */
		for (d=0; d < ncinfile->ndims; d++)
		{
			if (d==ncinfile->recdim)
				ncdimdef(ncoutfile->ncfid,ncinfile->dimname[d],NC_UNLIMITED);
			else ncdimdef(ncoutfile->ncfid,ncinfile->dimname[d],
					ncinfile->dimregsize[d]);
		}

		/* Define the variables and copy their attributes */
		for (v=0; v < ncinfile->nvars; v++)
		{
			ncvardef(ncoutfile->ncfid,ncinfile->varname[v],ncinfile->datatype[v],
					ncinfile->varndims[v],ncinfile->vardim[v]);
			if (n4) { ERR(nc_def_var_deflate(ncoutfile->ncfid,v,1,1,n4)) }
			for (n=0; n < ncinfile->natts[v]; n++)
			{
				ncattname(ncinfile->ncfid,v,n,attname);
				if (missing)
				{
					if (!strcmp(attname,"missing_value"))
					{
						ncoutfile->varmiss[v]=1;
						ncattget(ncinfile->ncfid,v,"missing_value",
								(void *)(ncoutfile->varmissval[v]));
					}
				}
				if (!strcmp(attname,"domain_decomposition")) continue;
				else
				{
					if (ncattcopy(ncinfile->ncfid,v,attname,ncoutfile->ncfid,v)==(-1))
					{
						fprintf(stderr,"Error: cannot copy variable \"%s\"'s attributes!\n",
								ncinfile->varname[v]);
						free(ncinfile); return(1);
					}
				}
			}
		}

		/* Copy the global attributes */
		for (n=0; n < ncinfile->ngatts; n++)
		{
			ncattname(ncinfile->ncfid,NC_GLOBAL,n,attname);
			if (!strcmp(attname,"NumFilesInSet")) continue;
			else if (!strcmp(attname,"filename"))
				ncattput(ncoutfile->ncfid,NC_GLOBAL,attname,NC_CHAR,
						strlen(outncname),(void *)outncname);
			else
			{
				if (ncattcopy(ncinfile->ncfid,NC_GLOBAL,attname,ncoutfile->ncfid,
							NC_GLOBAL)==(-1))
				{
					fprintf(stderr,"Error: cannot copy the file's global attributes!\n");
					return(1);
				}
			}
		}

		/* Definitions done */
		nc__enddef(ncoutfile->ncfid,headerpad,4,0,4);
	}

	/* Copy all data values of the dimensions and variables to memory */
	ncinfileerror=process_vars(ncinfile,ncoutfile,appendnc,nrecs,r,*nfiles,
			f,varbuf,vbuf,verbose,missing);

	/* Done */
	ncclose(ncinfile->ncfid); free(ncinfile); return(ncinfileerror);
}


/* Copy all data values in an input file at the current record to memory */
int process_vars(struct fileinfo *ncinfile, struct fileinfo *ncoutfile,
		unsigned char appendnc, int *nrecs, int r, int nfiles,
		int f, void *varbuf[], void *vbuf[], unsigned char verbose,
		unsigned char missing)
{
	int v, d, i, j, k, l, b, s;  /* Loop variables */
	int dimid;  /* ID of a dimension */
	void *values;  /* Current data values */
	long instart[MAX_NC_DIMS], outstart[MAX_NC_DIMS];  /* Data array sizes */
	long count[MAX_NC_DIMS];                           /*        "         */
	long long recsize;  /* Decomposed size of one record of a variable */
	long long recfullsize;  /* Non-decomposed size of one record of a variable */
	long long recregsize;  /* Non-decomposed regular size of one record of a variable */
	int varrecdim;  /* Variable's record dimension */
	static unsigned char first=1;  /* First time reading variables? */
	int imax, jmax, kmax, lmax;
	int imaxfull, jmaxfull, kmaxfull, lmaxfull;
	int imaxjmaxfull, imaxjmaxkmaxfull;
	int offset, ioffset, joffset, koffset, loffset;
	long long varbufsize, vbufsize;

	/* Check the number of records */

	first=1;

	if (*nrecs==1) *nrecs=ncinfile->dimsize[ncinfile->recdim];
	else
		if (ncinfile->dimsize[ncinfile->recdim] != *nrecs)
		{
			fprintf(stderr,"Error: different number of records than the first input file!\n");
			return(1);
		}

	/* Loop over all the variables */
	for (v=0; v < ncinfile->nvars; v++)
	{
		if (verbose > 1) printf("    variable = %s\n",ncinfile->varname[v]);

		/* Get read/write dimension sizes for the variable */
		recsize=1; recfullsize=1; varrecdim=(-1); recregsize=1;
		outstart[0]=0; outstart[1]=0; outstart[2]=0; outstart[3]=0;
		for (d=0; d < ncinfile->varndims[v]; d++)
		{
			if (ncinfile->vardim[v][d]==ncinfile->recdim)
			{
				count[d]=1; varrecdim=d;
			}
			else
			{
				count[d]=ncinfile->dimsize[ncinfile->vardim[v][d]];
				recsize*=count[d]; instart[d]=0;
				outstart[d]=ncinfile->dimstart[ncinfile->vardim[v][d]]-1;
				recfullsize*=ncinfile->dimfullsize[ncinfile->vardim[v][d]];
				recregsize*=ncinfile->dimregsize[ncinfile->vardim[v][d]];
			}
			if (verbose > 1)
				printf("      dim %d:  instart=%ld  outstart=%ld  count=%ld\n",d,
						instart[d],outstart[d],count[d]);
		}

		/* Prevent unnecessary reads/writes */
		if (r > 0)
		{
			/* Prevent unnecessary reads/writes of the dimensions */
			if ((dimid=ncdimid(ncinfile->ncfid,ncinfile->varname[v]))!=(-1))
			{
				if (ncinfile->recdim==dimid)
				{
					if (f!=0) continue;
				}
				else continue;
			}
			/* Prevent unnecessary reads/writes of the variables */
			else
			{
				/* Prevent unnecessary reads/writes of non-decomposed variables
				   if (ncinfile->vardecomp[v]!=1 && appendnc) continue; */

				/* Non-record variables */
				if (varrecdim==(-1)) continue;

				/* Non-decomposed record variables */
				if (ncinfile->vardecomp[v]!=1 && f > 0) continue;
			}
		}
		else
		{
			if (ncinfile->vardecomp[v]!=1 && appendnc) continue;
		}

		/* Allocate a buffer for the variable's record */
		if ((values=malloc(nctypelen(ncinfile->datatype[v])*recsize))==NULL)
		{
			fprintf(stderr,"Error: cannot allocate %lld bytes for decomposed variable \"%s\"'s values!\n",
					nctypelen(ncinfile->datatype[v])*recsize,ncinfile->varname[v]);
			return(1);
		}

		/* Read the variable */
		if (varrecdim!=(-1)) instart[varrecdim]=outstart[varrecdim]=r;
		if (ncvarget(ncinfile->ncfid,v,instart,count,values)==(-1))
		{
			fprintf(stderr,"Error: cannot read variable \"%s\"'s values!\n",
					ncinfile->varname[v]);
			return(1);
		}

		/* Write the buffered variable immediately if it's not decomposed */
		if (ncinfile->vardecomp[v]!=1)
		{
			if (verbose > 1)
				printf("      writing %lld bytes to file\n",
						nctypelen(ncinfile->datatype[v])*recsize);
			if (ncvarput(ncoutfile->ncfid,v,outstart,count,values)==(-1))
			{
				fprintf(stderr,"Error: cannot write variable \"%s\"'s values!\n",
						ncinfile->varname[v]);
				return(1);
			}
		}
		else if (unpack && (ncinfile->varislat[v]==1))
		{
			if (first) {
				if(verbose) printf(" writing lat in place of %s to file\n", ncinfile->varname[v]);
				outstart[0]=0; count[0]=nlat;
				if (ncvarput(ncoutfile->ncfid,v,outstart,count,latc)==(-1))
				{
					fprintf(stderr,"Error: cannot write variable \"%s\"'s values!\n",
							ncinfile->varname[v]);
					return(1);
				}
			}
		}
		else if (unpack && (ncinfile->varislon[v]==1))
		{
			if (first) {
				if(verbose) printf(" writing lon in place of %s to file\n", ncinfile->varname[v]);
				outstart[0]=0; count[0]=nlon;
				if (ncvarput(ncoutfile->ncfid,v,outstart,count,lonc)==(-1))
				{
					fprintf(stderr,"Error: cannot write variable \"%s\"'s values!\n",
							ncinfile->varname[v]);
					return(1);
				}
			}
		}
		/* Save the buffer */
		else
		{
			/* Allocate a buffer for the variable's non-decomposed record size */
			if (first)
			{
				varbufsize=nctypelen(ncinfile->datatype[v])*recfullsize;
				vbufsize=nctypelen(ncinfile->datatype[v])*recregsize;

				if (verbose > 0) printf(" allocating %lld bytes for full domain var %s\n",
						varbufsize, ncinfile->varname[v]);

				if ((varbuf[v]=calloc(varbufsize,1))==NULL)
				{
					fprintf(stderr,"Error: cannot allocate %lld bytes for entire variable \"%s\"'s values!\n",
							varbufsize,ncinfile->varname[v]); return(1);
				}

				if (unpack){ 
					if (verbose > 0) printf(" allocating %lld bytes for unpacked var %s\n",
							vbufsize, ncinfile->varname[v]);
					if ((vbuf[v]=calloc(vbufsize,1))==NULL)
					{
						fprintf(stderr,"Error: cannot allocate %lld bytes for unpacked variable \"%s\"'s values!\n",
								vbufsize,ncinfile->varname[v]); return(1);
					}
				}
				if (missing && ncoutfile->varmiss[v])
					switch (ncinfile->datatype[v])
					{
						case NC_BYTE:
						case NC_CHAR:
							for (s=0; s < recfullsize; s++)
								*((unsigned char *)(varbuf[v])+s)=
									*((unsigned char *)(ncoutfile->varmissval[v]));
							break;
						case NC_SHORT:
							for (s=0; s < recfullsize; s++)
								*((short *)(varbuf[v])+s)=
									*((short *)(ncoutfile->varmissval[v]));
							break;
						case NC_INT:
							for (s=0; s < recfullsize; s++)
								*((int *)(varbuf[v])+s)=
									*((int *)(ncoutfile->varmissval[v]));
							break;
						case NC_FLOAT:
							for (s=0; s < recfullsize; s++)
								*((float *)(varbuf[v])+s)=
									*((float *)(ncoutfile->varmissval[v]));
							break;
						case NC_DOUBLE:
							for (s=0; s < recfullsize; s++)
								*((double *)(varbuf[v])+s)=
									*((double *)(ncoutfile->varmissval[v]));
							break;
					}
			}
			if (varbuf[v]==NULL)
			{
				fprintf(stderr,"Internal memory usage error!\n"); return(1);
			}
			if (vbuf[v]==NULL)
			{
				fprintf(stderr,"Internal memory usage error!\n"); return(1);
			}
			if (verbose > 1)
				printf("      writing %lld bytes to memory\n",
						nctypelen(ncinfile->datatype[v])*recsize);

			imax=ncinfile->dimsize[ncinfile->vardim[v][ncinfile->varndims[v]-1]];
			if (ncinfile->varndims[v] > 1)
			{
				dimid=ncinfile->vardim[v][ncinfile->varndims[v]-2];
				if (dimid==ncinfile->recdim) jmax=1;
				else jmax=ncinfile->dimsize[dimid];
			}
			else jmax=1;
			if (ncinfile->varndims[v] > 2)
			{
				dimid=ncinfile->vardim[v][ncinfile->varndims[v]-3];
				if (dimid==ncinfile->recdim) kmax=1;
				else kmax=ncinfile->dimsize[dimid];
			}
			else kmax=1;
			if (ncinfile->varndims[v] > 3)
			{
				dimid=ncinfile->vardim[v][ncinfile->varndims[v]-4];
				if (dimid==ncinfile->recdim) lmax=1;
				else lmax=ncinfile->dimsize[dimid];
			}
			else lmax=1;
			if (verbose > 1)
				printf("      imax=%d  jmax=%d  kmax=%d  lmax=%d\n",imax,jmax,
						kmax,lmax);

			imaxfull=ncinfile->dimfullsize[ncinfile->vardim[v][ncinfile->varndims[v]-1]];
			if (ncinfile->varndims[v] > 1)
				jmaxfull=ncinfile->dimfullsize[ncinfile->vardim[v][ncinfile->varndims[v]-2]];
			else jmaxfull=1;
			if (ncinfile->varndims[v] > 2)
				kmaxfull=ncinfile->dimfullsize[ncinfile->vardim[v][ncinfile->varndims[v]-3]];
			else kmaxfull=1;
			if (ncinfile->varndims[v] > 3)
			{
				if (ncinfile->vardim[v][ncinfile->varndims[v]-4]!=ncinfile->recdim)
					lmaxfull=ncinfile->dimfullsize[ncinfile->vardim[v][ncinfile->varndims[v]-4]];
				else lmaxfull=1;
			}
			else lmaxfull=1;
			if (verbose > 1)
				printf("      imaxfull=%d  jmaxfull=%d  kmaxfull=%d  lmaxfull=%d\n",
						imaxfull,jmaxfull,kmaxfull,lmaxfull);
			imaxjmaxfull=imaxfull*jmaxfull;
			imaxjmaxkmaxfull=imaxfull*jmaxfull*kmaxfull;

			ioffset=outstart[ncinfile->varndims[v]-0-1];
			if (ncinfile->varndims[v] > 1)
				joffset=outstart[ncinfile->varndims[v]-1-1];
			else joffset=0;
			if (ncinfile->varndims[v] > 2)
				koffset=outstart[ncinfile->varndims[v]-2-1];
			else koffset=0;
			if (ncinfile->varndims[v] > 3)
				loffset=outstart[ncinfile->varndims[v]-3-1];
			else loffset=0;
			if (varrecdim!=(-1))
			{
				switch (ncinfile->varndims[v])
				{
					case 1:
						ioffset=0;
						break;
					case 2:
						joffset=0;
						break;
					case 3:
						koffset=0;
						break;
					case 4:
						loffset=0;
						break;
				}
			}
			if (verbose > 1)
				printf("      ioffset=%d  joffset=%d  koffset=%d  loffset=%d\n",
						ioffset,joffset,koffset,loffset);
			switch (ncinfile->datatype[v])
			{
				case NC_BYTE:
				case NC_CHAR:
					if (verbose > 1) printf("      start copying byte/char\n");
					b=0;
					for (l=0; l < lmax; l++)
						for (k=0; k < kmax; k++)
							for (j=0; j < jmax; j++)
								for (i=0; i < imax; i++)
								{
									offset=(i+ioffset)+
										(j+joffset)*imaxfull+
										(k+koffset)*imaxjmaxfull+
										(l+loffset)*imaxjmaxkmaxfull;
									*((unsigned char *)(varbuf[v])+offset)=
										*((unsigned char *)values+(b++));
								}
					if (verbose > 1) printf("      end copying byte/char\n");
					break;
				case NC_SHORT:
					if (verbose > 1) printf("      start copying short\n");
					b=0;
					for (l=0; l < lmax; l++)
						for (k=0; k < kmax; k++)
							for (j=0; j < jmax; j++)
								for (i=0; i < imax; i++)
								{
									offset=(i+ioffset)+
										(j+joffset)*imaxfull+
										(k+koffset)*imaxjmaxfull+
										(l+loffset)*imaxjmaxkmaxfull;
									*((short *)(varbuf[v])+offset)=
										*((short *)values+(b++));
								}
					if (verbose > 1) printf("      end copying short\n");
					break;
				case NC_INT:
					if (verbose > 1) printf("      start copying int\n");
					b=0;
					for (l=0; l < lmax; l++)
						for (k=0; k < kmax; k++)
							for (j=0; j < jmax; j++)
								for (i=0; i < imax; i++)
								{
									offset=(i+ioffset)+
										(j+joffset)*imaxfull+
										(k+koffset)*imaxjmaxfull+
										(l+loffset)*imaxjmaxkmaxfull;
									*((int *)(varbuf[v])+offset)=
										*((int *)values+(b++));
								}
					if (verbose > 1) printf("      end copying int\n");
					break;
				case NC_FLOAT:
					if (verbose > 1) printf("      start copying float\n");
					b=0;
					for (l=0; l < lmax; l++)
						for (k=0; k < kmax; k++)
							for (j=0; j < jmax; j++)
								for (i=0; i < imax; i++)
								{
									offset=(i+ioffset)+
										(j+joffset)*imaxfull+
										(k+koffset)*imaxjmaxfull+
										(l+loffset)*imaxjmaxkmaxfull;
									*((float *)(varbuf[v])+offset)=
										*((float *)values+(b++));
								}
					if (verbose > 1) printf("      end copying float\n");
					break;
				case NC_DOUBLE:
					if (verbose > 1) printf("      start copying double\n");
					b=0;
					for (l=0; l < lmax; l++)
						for (k=0; k < kmax; k++)
							for (j=0; j < jmax; j++)
								for (i=0; i < imax; i++)
								{
									offset=(i+ioffset)+
										(j+joffset)*imaxfull+
										(k+koffset)*imaxjmaxfull+
										(l+loffset)*imaxjmaxkmaxfull;
									*((double *)(varbuf[v])+offset)=
										*((double *)values+(b++));
								}
					if (verbose > 1) printf("      end copying double\n");
					break;
			}
		}

		/* Deallocate the decomposed variable's buffer */
		free(values);
	}
	first=0; return(0);
}


int flush_decomp(struct fileinfo *ncoutfile, int nfiles, int r,
		void *varbuf[], void *vbuf[], unsigned char verbose)
{
	int v, d;  /* Loop variable */
	long outstart[MAX_NC_DIMS];  /* Data array sizes */
	long count[MAX_NC_DIMS];     /*        "         */
	int varrecdim;  /* Position of a variable's record dimension */
	int imaxp, jmaxp, kmaxp, lmaxp;
	int imaxjmaxp, imaxjmaxkmaxp;
	int imaxr, jmaxr, kmaxr, lmaxr;
	int imaxjmaxr, imaxjmaxkmaxr;
	int i, l, n, k, poffset, roffset, varhaslev;

	if (verbose > 0)
	{
		printf("    nvars=%d, r=%d\n",ncoutfile->nvars,r);
	}

	/* Write out all the decomposed variables */
	for (v=0; v < ncoutfile->nvars; v++)
	{
		if (ncoutfile->vardecomp[v]==0) continue;
		if (ncoutfile->varislat[v]==1) continue;
		if (ncoutfile->varislon[v]==1) continue;
		if (verbose > 0) printf("    v=%d (%s)\n",v,ncoutfile->varname[v]);
		varrecdim=(-1);
		varhaslev=0;
		for (d=0; d < ncoutfile->varndims[v]; d++)
		{
			outstart[d]=0;
			if (ncoutfile->vardim[v][d]==ncoutfile->recdim)
			{
				count[d]=1; varrecdim=d;
			}
			else
			{
				if(unpack)count[d]=ncoutfile->dimregsize[ncoutfile->vardim[v][d]];
				else count[d]=ncoutfile->dimfullsize[ncoutfile->vardim[v][d]];
			}
			if (verbose > 0)
				printf("      d=%d:  outstart=%ld  count=%ld\n",d,outstart[d],
						count[d]);
		}

		if (ncoutfile->varndims[v]>3) varhaslev=1;
		else if ((ncoutfile->varndims[v]==3) && (varrecdim==(-1))) varhaslev=1;

		if (varrecdim!=(-1)) outstart[varrecdim]=r;
		if (varrecdim==(-1) && r > 0) continue;

		//combined -> unpacked (varbuf -> vbuf)
		if (unpack) {
			imaxp=ncoutfile->dimfullsize[ncoutfile->vardim[v][ncoutfile->varndims[v]-1]];

			if (ncoutfile->varndims[v] > 1)
				jmaxp=ncoutfile->dimfullsize[ncoutfile->vardim[v][ncoutfile->varndims[v]-2]];
			else jmaxp=1;

			if (ncoutfile->varndims[v] > 2)
			{
				if (ncoutfile->vardim[v][ncoutfile->varndims[v]-3]!=ncoutfile->recdim)
					kmaxp=ncoutfile->dimfullsize[ncoutfile->vardim[v][ncoutfile->varndims[v]-3]];
				else kmaxp=1;
			}
			else kmaxp=1;

			if (ncoutfile->varndims[v] > 3)
			{
				if (ncoutfile->vardim[v][ncoutfile->varndims[v]-4]!=ncoutfile->recdim)
					lmaxp=ncoutfile->dimfullsize[ncoutfile->vardim[v][ncoutfile->varndims[v]-4]];
				else lmaxp=1;
			}
			else lmaxp=1;

			if (verbose)
				printf(" %s  imaxp=%d  jmaxp=%d  kmaxp=%d  lmaxp=%d\n",
						ncoutfile->varname[v],imaxp,jmaxp,kmaxp,lmaxp);

			imaxjmaxp=imaxp*jmaxp;
			imaxjmaxkmaxp=imaxp*jmaxp*kmaxp;

			imaxr=ncoutfile->dimregsize[ncoutfile->vardim[v][ncoutfile->varndims[v]-1]];

			if (ncoutfile->varndims[v] > 1)
				jmaxr=ncoutfile->dimregsize[ncoutfile->vardim[v][ncoutfile->varndims[v]-2]];
			else jmaxr=1;

			if (ncoutfile->varndims[v] > 2)
			{
				if (ncoutfile->vardim[v][ncoutfile->varndims[v]-3]!=ncoutfile->recdim)
					kmaxr=ncoutfile->dimregsize[ncoutfile->vardim[v][ncoutfile->varndims[v]-3]];
				else kmaxr=1;
			}
			else kmaxr=1;

			if (ncoutfile->varndims[v] > 3)
			{
				if (ncoutfile->vardim[v][ncoutfile->varndims[v]-4]!=ncoutfile->recdim)
					lmaxr=ncoutfile->dimregsize[ncoutfile->vardim[v][ncoutfile->varndims[v]-4]];
				else lmaxr=1;
			}
			else lmaxr=1;

			if (verbose)
				printf(" %s  imaxr=%d  jmaxr=%d  kmaxr=%d  lmaxr=%d\n",
						ncoutfile->varname[v],imaxr,jmaxr,kmaxr,lmaxr);

			imaxjmaxr=imaxr*jmaxr;
			imaxjmaxkmaxr=imaxr*jmaxr*kmaxr;

			if (varhaslev!=0)
			{
				if (!((lmaxr==lmaxp) && (imaxr=imaxp))) {
					fprintf(stderr,"dimension mismatch lmax = %d, %d, imax = %d, %d\n", 
							lmaxr,lmaxp,imaxr,imaxp);
					return(1);
				}

				switch (ncoutfile->datatype[v])
				{
					case NC_SHORT:
						if (verbose > 1) printf("      start unpacking short\n");
						for (n=0; n < imaxjmaxkmaxr; n++) { *((short *)(vbuf[v])+n) = 0; } 
						for (l=0; l < lmaxp; l++)
							for (n=0; n < xn; n++)
								for (i=0; i < imaxp; i++)
								{
									poffset=i+
										pxj[n]*imaxp+
										pxi[n]*imaxjmaxp+
										l*imaxjmaxkmaxp;

									roffset=i+
										rxj[n]*imaxr+
										rxi[n]*imaxjmaxr+
										l*imaxjmaxkmaxr;

									*((short *)(vbuf[v])+roffset)=
										*((short *)(vbuf[v])+roffset) +
										*((short *)(varbuf[v])+poffset) *
										rxf[n];
								}
						if (verbose > 1) printf("      end unpacking short\n");
						break;
					case NC_INT:
						if (verbose > 1) printf("      start unpacking int\n");
						for (n=0; n < imaxjmaxkmaxr; n++) { *((int *)(vbuf[v])+n) = 0; } 
						for (l=0; l < lmaxp; l++)
							for (n=0; n < xn; n++)
								for (i=0; i < imaxp; i++)
								{
									poffset=i+
										pxj[n]*imaxp+
										pxi[n]*imaxjmaxp+
										l*imaxjmaxkmaxp;

									roffset=i+
										rxj[n]*imaxr+
										rxi[n]*imaxjmaxr+
										l*imaxjmaxkmaxr;

									*((int *)(vbuf[v])+roffset)=
										*((int *)(vbuf[v])+roffset) +
										*((int *)(varbuf[v])+poffset) *
										rxf[n];
								}
						if (verbose > 1) printf("      end unpacking int\n");
						break;
					case NC_FLOAT:
						if (verbose > 1) printf("      start unpacking float\n");
						for (n=0; n < imaxjmaxkmaxr; n++) { *((float *)(vbuf[v])+n) = 0; } 
						for (l=0; l < lmaxp; l++)
							for (n=0; n < xn; n++)
								for (i=0; i < imaxp; i++)
								{
									poffset=i+
										pxj[n]*imaxp+
										pxi[n]*imaxjmaxp+
										l*imaxjmaxkmaxp;

									roffset=i+
										rxj[n]*imaxr+
										rxi[n]*imaxjmaxr+
										l*imaxjmaxkmaxr;

									*((float *)(vbuf[v])+roffset)=
										*((float *)(vbuf[v])+roffset) +
										*((float *)(varbuf[v])+poffset) *
										rxf[n];
								}
						if (verbose > 1) printf("      end unpacking float\n");
						break;
					case NC_DOUBLE:
						if (verbose > 1) printf("      start unpacking double\n");
						for (n=0; n < imaxjmaxkmaxr; n++) { *((double *)(vbuf[v])+n) = 0; } 
						for (l=0; l < lmaxp; l++)
							for (n=0; n < xn; n++)
								for (i=0; i < imaxp; i++)
								{
									poffset=i+
										pxj[n]*imaxp+
										pxi[n]*imaxjmaxp+
										l*imaxjmaxkmaxp;

									roffset=i+
										rxj[n]*imaxr+
										rxi[n]*imaxjmaxr+
										l*imaxjmaxkmaxr;

									*((double *)(vbuf[v])+roffset)=
										*((double *)(vbuf[v])+roffset) +
										*((double *)(varbuf[v])+poffset) *
										rxf[n];
								}
						if (verbose > 1) printf("      end unpacking double\n");
						break;
				} //switch case
			} // if varhaslev
			else 
			{
				if (!((lmaxr==lmaxp) && (kmaxr=kmaxp))) {
					fprintf(stderr,"dimension mismatch lmax = %d, %d, kmax = %d, %d\n", 
							lmaxr,lmaxp,kmaxr,kmaxp);
					return(1);
				}

				switch (ncoutfile->datatype[v])
				{
					case NC_SHORT:
						if (verbose > 1) printf("      start unpacking short\n");
						for (n=0; n < imaxjmaxkmaxr; n++) { *((short *)(vbuf[v])+n) = 0; } 
						for (l=0; l < lmaxp; l++)
							for (k=0; k < kmaxp; k++)
								for (n=0; n < xn; n++)
								{
									poffset=pxj[n]+
										pxi[n]*imaxp+
										k*imaxjmaxp+
										l*imaxjmaxkmaxp;

									roffset=rxj[n]+
										rxi[n]*imaxr+
										k*imaxjmaxr+
										l*imaxjmaxkmaxr;

									*((short *)(vbuf[v])+roffset)=
										*((short *)(vbuf[v])+roffset) +
										*((short *)(varbuf[v])+poffset) *
										rxf[n];
								}
						if (verbose > 1) printf("      end unpacking short\n");
						break;
					case NC_INT:
						if (verbose > 1) printf("      start unpacking int\n");
						for (n=0; n < imaxjmaxkmaxr; n++) { *((int *)(vbuf[v])+n) = 0; } 
						for (l=0; l < lmaxp; l++)
							for (k=0; k < kmaxp; k++)
								for (n=0; n < xn; n++)
								{
									poffset=pxj[n]+
										pxi[n]*imaxp+
										k*imaxjmaxp+
										l*imaxjmaxkmaxp;

									roffset=rxj[n]+
										rxi[n]*imaxr+
										k*imaxjmaxr+
										l*imaxjmaxkmaxr;

									*((int *)(vbuf[v])+roffset)=
										*((int *)(vbuf[v])+roffset) +
										*((int *)(varbuf[v])+poffset) *
										rxf[n];
								}
						if (verbose > 1) printf("      end unpacking int\n");
						break;
					case NC_FLOAT:
						if (verbose > 1) printf("      start unpacking float\n");
						for (n=0; n < imaxjmaxkmaxr; n++) { *((float *)(vbuf[v])+n) = 0; } 
						for (l=0; l < lmaxp; l++)
							for (k=0; k < kmaxp; k++)
								for (n=0; n < xn; n++)
								{
									poffset=pxj[n]+
										pxi[n]*imaxp+
										k*imaxjmaxp+
										l*imaxjmaxkmaxp;

									roffset=rxj[n]+
										rxi[n]*imaxr+
										k*imaxjmaxr+
										l*imaxjmaxkmaxr;

									*((float *)(vbuf[v])+roffset)=
										*((float *)(vbuf[v])+roffset) +
										*((float *)(varbuf[v])+poffset) *
										rxf[n];
								}
						if (verbose > 1) printf("      end unpacking float\n");
						break;
					case NC_DOUBLE:
						if (verbose > 1) printf("      start unpacking double\n");
						for (n=0; n < imaxjmaxkmaxr; n++) { *((double *)(vbuf[v])+n) = 0; } 
						for (l=0; l < lmaxp; l++)
							for (k=0; k < kmaxp; k++)
								for (n=0; n < xn; n++)
								{
									poffset=pxj[n]+
										pxi[n]*imaxp+
										k*imaxjmaxp+
										l*imaxjmaxkmaxp;

									roffset=rxj[n]+
										rxi[n]*imaxr+
										k*imaxjmaxr+
										l*imaxjmaxkmaxr;

									*((double *)(vbuf[v])+roffset)=
										*((double *)(vbuf[v])+roffset) +
										*((double *)(varbuf[v])+poffset) *
										rxf[n];
								}
						if (verbose > 1) printf("      end unpacking double\n");
						break;
				} //switch case
			} // if varhaslev
		}//if unpack

		//done unpacking
		if (verbose) printf("      writing to disk\n");
		if (unpack) {
			if (ncvarput(ncoutfile->ncfid,v,outstart,count,vbuf[v])==(-1))
			{
				fprintf(stderr,"Error: cannot write variable \"%d\"'s values!\n",
						v); return(1);
			}
		}
		else
		{
			if (ncvarput(ncoutfile->ncfid,v,outstart,count,varbuf[v])==(-1))
			{
				fprintf(stderr,"Error: cannot write variable \"%d\"'s values!\n",
						v); return(1);
			}
		}

		if (verbose) printf("      done writing to disk\n");
	}
	return(0);
}


double modtimediff (char *filename[]) {
	const int narg=2;
	struct stat foo1, foo2;
	time_t mtime1, mtime2;
	double diff_t;

	if (stat(filename[0], &foo1) < 0) { perror(filename[0]); return (-1); }
	mtime1 = foo1.st_mtime; /* seconds since the epoch */

	if (stat(filename[1], &foo2) < 0) { perror(filename[1]); return (-1); }
	mtime2 = foo2.st_mtime; 

	diff_t = difftime(mtime2,mtime1);
	return(diff_t);
}

int show_status (double percent)
{
    int x;
	double percent1;

	percent1 = percent;
	//if (percent1>100.) percent1=100.;

    //for(x = 0; x < percent1; x++)
    //{
    //   printf("|");

    //}
       printf("%.2f%%\r", percent1);
       fflush(stdout);

    return(EXIT_SUCCESS);
}

double modtime (char *filename[]) {
	struct stat foo1, foo2;
	time_t mtime1, mtime2;
	struct tm str_time;
	double diff_t;

	str_time.tm_year = 1987-1900;
	str_time.tm_mon = 3;
	str_time.tm_mday = 26;
	str_time.tm_hour = 0;
	str_time.tm_min = 0;
	str_time.tm_sec = 0;
	str_time.tm_isdst = 0;

	mtime1 = mktime(&str_time);

	if (stat(filename[0], &foo2) < 0) { perror(filename[0]); return (-1); }
	mtime2 = foo2.st_mtime; 

	diff_t = difftime(mtime2,mtime1);
	return(diff_t);
}

int rmfile(char *filename[]) {
	return(remove(filename[0]));
}

