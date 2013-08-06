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
	int num;
	int* num_list;
} PROTO;

typedef struct
{
	int num_cg_sites, num_cg_types, num_fg_sites, num_fg_types;
	int max_to_map, num_frames;
	//int threshold;
	PROTO* prototype;
	
	int geometry_map_flag, observable_map_flag, num_observables;
	int output_flag;
	int sensitivity_flag;
	int map_style_flag;
	int debug_flag;
	int sign_flag;
	
	int num_map;
	int frame;
	int* map;

	int log_type;
	double scaleU1;
	double scaleU2;
	double scaleF;
	int guess_type;
	int sens_map_flag;
	double log_value;	//should be able to write this variable out
	double guess;		//should be able to write this variable out
	Filenames files;
	double timestep, volume;	//only used for formatting output logfile
	
	double* log_values;
	double* guesses;
	
	int* charge;
	int num_charges, num_files, num_outfile;
	FILE** file_point;
	FILE** outfile;
	char* name;
} Controller;

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
	int type_count;
	
} Frame;

#endif