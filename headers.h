//headers.h
#ifndef HEADERFILE_H
#define HEADERFILE_H

//define structures
typedef struct
{
	char dump1[64]; 
	char dump2[64];
	char log[64];
	char guess[64];
} Filenames;

typedef struct
{
	int id, mol, type;
	double q, mass, x, y, z;
	double* observables;
} ATOM;

typedef struct
{
	double x, y, z, mass;
	int type;
} COORD;

typedef struct
{
	int id, mol, type;
	double q, mass, x, y, z;
	int num_in_site;
	double* observables;
	COORD* coord;
	int* matches;
} SITE;

typedef struct
{
	double xmin, xmax, ymin, ymax, zmin, zmax;
	int timestep;
	int num_atoms;
	int num_mol;
	int num_observables;
	ATOM* atoms;
	SITE* sites;
	int* type;
	int* type_num;
	int** type_list;
	int type_count;
	
} Frame;

typedef struct
{
	int num_cg_sites, num_cg_types, num_frames;
	PROTO* prototype;
	
	int num_observables;
	int output_flag;
	int sensitivity_flag;
	void (*read_function)(Frame*, FILE*);
	void (*output_function)(Frame*, FILE*);
	void (*header_function)(FILE*);
	
	int frame;
	int log_type;
	double scaleU;
	int guess_type;
	double log_value;	//should be able to write this variable out
	double guess;		//should be able to write this variable out
	Filenames files;
	double timestep, volume;	//only used for formatting output logfile
	
	double* log_values;
	double* guesses;
	
	double* charge;
	int num_charges, num_files, num_outfile;
	FILE** file_point;
	FILE** outfile;
	char* name;
	
} Controller;
#endif