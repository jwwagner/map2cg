//MAP2CG
//By Jacob Wagner
// Created July 8, 2013
//Creates a data file and a topology file where a specified coarse-grained mapping
//has been applied to the frames of a custom LAMMPS dump.
//Useful in mapping an observable(s) from a fine-grained (FG) to a coarse-grained (CG) system.

//Acknowledgements: Some intput handling functions are taken from James Dama's MS-CG Force Matching Code
// Additionally, one of the outputs is simillar to the A2CG program by Zhen that uses GROMAC's libraries.

//INPUT: Custom LAMMPS DUMP (see next line) and a TOPOLOGY mapping file
//1: dump tj2 all custom 250 MeOH_sample.dat	id	mol type q mass x y z fx fy fz

//2: topology file should be of the form:
// int	#number of CG sites
// int  #number CG types
// int	#number of FG sites
// int  #number FG types
// int  #max number of FG sites in single CG site
//		#this line intentionally left blank
// int	#geometry mapping flag (0 = CoM geometry mapping, 1 = center of geometry mapping)
// int  #observable mapping flag (0 = additive mappign)
// int	#number observables to map
//	 	#this line intentionally left blank
// int  #mapping style flag (0 = 1:1 molecule entirely, 1 = implicit solvent)
//IF mapping style flag == 1
// int	#number of types to map
// int #type ID's to be mapped
// ...
//ENDIF
//		#end of input file, please left blank line at end

//call with one of 4 options
// exec.x #assumes input files dump.dat and map.top
// exec.x -f name.dat #assumes map.top
// exec.x -f1 name.top #assumes dump.dat
// exec.x -f name.dat -f1 name.top

//includes
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <signal.h>
#include <limits.h>
#include <time.h>

//define structures
typedef struct
{
	int num_cg_sites, num_cg_types, num_fg_sites, num_fg_types;
	int max_to_map;
	int geometry_map_flag, observable_map_flag, num_observables;
	int map_style_flag;
	int num_map;
	int frame;
	int* map;
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
} COORD;

typedef struct
{
	int id, mol, type;
	double q, mass, x, y, z;
	int num_in_site;
	double* observables;
	COORD* coord;
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

//function prototype definitions
void parse_command_line_arguments(int, char**, char*, char*, char*);
void read_topology_file(Controller*, char*);
void read_frame(Controller*, Frame*, FILE*, int*);
void process_frame(Controller*, Frame*, Frame*);
void output_frame(Controller*, Frame*, char*);
void output_topology(Controller*, Frame*);
void check_file_extension(char*, char*);
void report_traj_input_suffix_error(char*);
void report_usage_error(char*);
void free_allocation(Controller*, Frame*, Frame*);

//begin main function
int main(int argc,char *argv[])
{

	//declare variables
   	double start_cputime = clock();
   	Controller controls;
   	Frame frame;
   	Frame outframe;
   	frame.num_atoms = 0;
	outframe.num_atoms = 0;

	char* datfile = malloc(64 * sizeof(char));
	char* topfile = malloc(64 * sizeof(char));
	char* outfile = malloc(64 * sizeof(char));
	FILE* df;
	FILE* of;
	int frame_count = 0;
	int cont_flag = 1;
	controls.frame = 0;
	
	//read topology/control file
	parse_command_line_arguments(argc, argv, datfile, topfile, outfile);
 
	printf("assigned file names are %s and %s\n", datfile, topfile);
	
	//read topology/control file
	read_topology_file(&controls, topfile);
	
	//open dump file and toggle output file to reset
	df = fopen(datfile, "rt");
	of = fopen(outfile, "w+");
	fclose(of);
	
	//read first frame
	frame_count++;
	controls.frame++;
	read_frame(&controls, &frame, df, &cont_flag);
	
	printf("main says num_observables is %d\n", controls.num_observables);
	
	while(cont_flag == 1)
		{
		//process/map frame
		process_frame(&controls, &frame, &outframe);
		printf("finished processing frame\n");
		
		//output mapped frame and observables
		output_frame(&controls, &outframe, outfile);
		printf("finished output for frame %d\n", frame_count);
		
		//read next frame or set flag if done
		frame_count++;
		controls.frame++;
		read_frame(&controls, &frame, df, &cont_flag);
		printf("cont_flag is %d\n", cont_flag);
		}
		
	//close dump file
	fclose(df);
	
	//create top.in file for FM
	output_topology(&controls, &outframe);
	
	//free allocated variables
	free(datfile);
	free(topfile);
	free(outfile);
	if( (controls.map_style_flag == 0) || (controls.map_style_flag == 1) ) free(controls.map);
	free_allocation(&controls, &frame, &outframe);
}

//should add support for -o flag to specify output
void parse_command_line_arguments(int num_arg, char** arg, char* filename1, char* filename2, char* outfile)
{
    if (num_arg != 1 && num_arg != 3 && num_arg != 5 && num_arg != 7) report_usage_error(arg[0]);
    else if (num_arg == 3) 
    	{
        if (strcmp(arg[1], "-f") == 0)
        	{
        	sscanf(arg[2], "%s", filename1);
        	check_file_extension(arg[2], "dat");
        	strcpy(filename2, "map.top");
        	strcpy(outfile, "out.dat");
        	}
        else if  (strcmp(arg[1], "-f1") == 0)
        	{
        	sscanf(arg[2], "%s", filename2);
        	check_file_extension(arg[2], "top");
        	strcpy(filename1, "dump.dat");
        	strcpy(outfile, "out.dat");
        	}
        else if	(strcmp(arg[1], "-o") == 0)
        	{
        	sscanf(arg[2], "%s", outfile);
        	//check_file_extension(arg[2], "dat");
        	strcpy(filename1, "dump.dat");
        	strcpy(filename2, "map.top");
        	}
        else  report_usage_error(arg[0]);
    	} 
    else if (num_arg == 5) 
    	{
        if ( ((strcmp(arg[1], "-f") != 0) || (strcmp(arg[3], "-f1") != 0)) && ((strcmp(arg[1], "-f1") != 0) || (strcmp(arg[3], "-o") != 0)) && ((strcmp(arg[1], "-f") != 0) || (strcmp(arg[3], "-o") != 0)) ) report_usage_error(arg[0]);
        else if ( strcmp(arg[1], "-f") == 0 && strcmp(arg[3], "-f1") == 0 )
        	{
        	sscanf(arg[2], "%s", filename1);
        	sscanf(arg[4], "%s", filename2);
        	check_file_extension(arg[2], "dat");
        	check_file_extension(arg[4], "top");
        	strcpy(outfile, "out.dat");
        	}
        else if ( strcmp(arg[1], "-f") == 0 && strcmp(arg[3], "-o") == 0 )
        	{
        	sscanf(arg[2], "%s", filename1);
        	sscanf(arg[4], "%s", outfile);
        	check_file_extension(arg[2], "dat");
        	//check_file_extension(arg[4], "dat");
        	strcpy(filename2, "map.top");
        	}
         else if (strcmp(arg[1], "-f1") == 0 && strcmp(arg[3], "-o") == 0 )
     		{
     		sscanf(arg[2], "%s", filename1);
        	sscanf(arg[4], "%s", outfile);
        	check_file_extension(arg[2], "top");
        	//check_file_extension(arg[4], "dat");
     		strcpy(filename2, "map.top");
     		}
     	}
	else if (num_arg == 7) 
    	{
        if (strcmp(arg[1], "-f") != 0 || strcmp(arg[3], "-f1") != 0 || strcmp(arg[5], "-o") != 0) report_usage_error(arg[0]);
		else
			{
			sscanf(arg[2], "%s", filename1);
        	sscanf(arg[4], "%s", filename2);
        	sscanf(arg[6], "%s", outfile);
        	check_file_extension(arg[2], "dat");
        	check_file_extension(arg[4], "top");
        	//check_file_extension(arg[6], "dat");
			}
		}
     else if (num_arg == 1) 
     	{
        strcpy(filename1, "dump.dat");
        strcpy(filename2, "map.top");
        strcpy(outfile, "out.dat");
     	}
}

void read_topology_file(Controller *control, char* topfile)
{
	int i;
	FILE* fr = fopen(topfile, "rt");
    char line[100];
    
    fgets(line,100,fr);//1
    sscanf(line, "%d", &control->num_cg_sites);
    
    fgets(line,100,fr);//2
    sscanf(line, "%d", &control->num_cg_types);
    
    fgets(line,100,fr);//3
    sscanf(line, "%d", &control->num_fg_sites);
    
    fgets(line,100,fr);//4
    sscanf(line, "%d", &control->num_fg_types);
	
    fgets(line,100,fr);//5
    sscanf(line, "%d", &control->max_to_map);
	
    fgets(line,100,fr);//6
    
    fgets(line,100,fr);//7
    sscanf(line, "%d", &control->geometry_map_flag);
    
    fgets(line,100,fr);//8
    sscanf(line, "%d", &control->observable_map_flag);
    
    fgets(line,100,fr);//9
    sscanf(line, "%d", &control->num_observables);
    
    fgets(line,100,fr);//10
    
    fgets(line,100,fr);//11
    sscanf(line, "%d", &control->map_style_flag);
    
    //print information read in file
    printf("num_cg_sites %d\n", control->num_cg_sites);
    printf("num_cg_types %d\n", control->num_cg_types);
    printf("num_fg_sites %d\n", control->num_fg_sites);
    printf("num_fg_types %d\n", control->num_fg_types);
    printf("max_to_map %d\n", control->max_to_map);
    printf("\n");
    printf("geometry_map_flag %d\n", control->geometry_map_flag);
    printf("observable_map_flag %d\n", control->observable_map_flag);
    printf("num_observables %d\n", control->num_observables);
    printf("\n");
    printf("map_style_flag %d\n", control->map_style_flag);
    
    if(control->map_style_flag == 0)
    	{
    	control->num_map = control->num_fg_types;
    	control->map = malloc(control->num_fg_types * sizeof(int));
    	for(i=0; i < control->num_fg_types; i++)
    		{
    		control->map[i] = i+1;
    		}
    	}
    if(control->map_style_flag == 1)
    	{
    	fgets(line,100,fr);//12
    	sscanf(line, "%d", &control->num_map);
    	control->map = malloc(control->num_map * sizeof(int));
    	for(i=0; i < control->num_map; i++)
    		{
    		fgets(line,100,fr);
    		sscanf(line, "%d", &control->map[i]);
    		}
    	}
    	
    printf("finished reading top file\n\n");
}


void read_frame(Controller* control, Frame* frame, FILE* df, int* flag)
{
	int i, j;
	char line[100];
	char test[15];
	frame->num_observables = control->num_observables;
	
	printf("begin reading frame\n");
	//check if EOF
	if( fgets (line, 100, df) == NULL )
		{
		*flag = 0;
		printf("end of file reached\n");
		return;
		}
	
	//check if content matches expected 1st line of frame
	memcpy( test, &line[0], 14);
	test[14] = '\0';
	if( strcmp(test, "ITEM: TIMESTEP") != 0)
		{
			printf("error in frame format\n");
			printf("line found was %s\n", test);
			printf("expected line is ITEM: TIMESTEP\n");
			*flag = 0;
			return;
		}
	
	printf("passed EOF and format test\n");
	//so, we assume content is correct now for next frame and begin reading
	//read in timestep
	fgets(line,100,df);
    sscanf(line, "%d", &frame->timestep);
    
    //read in number of atoms
    fgets(line,100,df);
	fgets(line,100,df);
    sscanf(line, "%d", &i);
    
    printf("timestep is %d\n", frame->timestep);
    printf("number of atoms is %d\n", i);
    
    //check if num_atoms changed
    if(i != frame->num_atoms)
    	{
    	if(frame->num_atoms == 0)	//do initial allocation
    		{
    		printf("intial allocation of atoms \n");
    		frame->num_atoms = i;
			frame->atoms = malloc(frame->num_atoms * sizeof(ATOM));
			frame->num_observables = control->num_observables;
			printf("set frame->num_observables to %d\n", frame->num_observables);
			}
		else //number of atoms has changed
			{
			printf("number of atoms has  changed");
			//free all existing data for atoms
			for(j = 0; j < frame->num_atoms; j++)
				{
				free(frame->atoms[j].observables);	
				}
			free(frame->atoms);
			
			//reallocate atoms space
			frame->num_atoms = i;
			frame->atoms = malloc(frame->num_atoms * sizeof(ATOM));
			frame->num_mol = -1;
			}
			
		//allocate observable information
		for(i=0; i < frame->num_atoms; i++)
			{
			frame->atoms[i].observables = malloc(frame->num_observables * sizeof(double));
			}
		}
	
	//read in box information
	fgets(line,100,df);
	fgets(line,100,df);
    sscanf(line, "%lf %lf", &frame->xmin, &frame->xmax);
	fgets(line,100,df);
    sscanf(line, "%lf %lf", &frame->ymin, &frame->ymax);
	fgets(line,100,df);
    sscanf(line, "%lf %lf", &frame->zmin, &frame->zmax);
    
    printf("box is %lf %lf by %lf %lf by %lf %lf\n", frame->xmin, frame->xmax, frame->ymin, frame->ymax, frame->zmin, frame->zmax);
    //read in atoms and observables for frame
	fgets(line,100,df);
	for(i=0; i < frame->num_atoms; i++)
		{
		fgets(line,100,df);
		//determine reading based on how many arguements were supplied (expected number is frame->num_obserables + 8)
		switch(frame->num_observables)
		{
			case 0:
    			sscanf(line, "%d %d %d %lf %lf %lf %lf %lf", &frame->atoms[i].id, \
    			&frame->atoms[i].mol, &frame->atoms[i].type, &frame->atoms[i].q, &frame->atoms[i].mass, \
    			&frame->atoms[i].x, &frame->atoms[i].y, &frame->atoms[i].z);
				break;
			
			case 1:
				sscanf(line, "%d %d %d %lf %lf %lf %lf %lf %lf", &frame->atoms[i].id, \
				&frame->atoms[i].mol, &frame->atoms[i].type, &frame->atoms[i].q, &frame->atoms[i].mass, \
				&frame->atoms[i].x, &frame->atoms[i].y, &frame->atoms[i].z, &frame->atoms[i].observables[0]);
				break;
			case 2:
				sscanf(line, "%d %d %d %lf %lf %lf %lf %lf %lf %lf", &frame->atoms[i].id, \ 
				&frame->atoms[i].mol, &frame->atoms[i].type, &frame->atoms[i].q, &frame->atoms[i].mass, \
				&frame->atoms[i].x, &frame->atoms[i].y, &frame->atoms[i].z, &frame->atoms[i].observables[0], \
				&frame->atoms[i].observables[1]);
				break;
				
			case 3:
				sscanf(line, "%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf", &frame->atoms[i].id, \
				&frame->atoms[i].mol, &frame->atoms[i].type, &frame->atoms[i].q, &frame->atoms[i].mass, \
				&frame->atoms[i].x, &frame->atoms[i].y, &frame->atoms[i].z, &frame->atoms[i].observables[0], \ 
				&frame->atoms[i].observables[1], &frame->atoms[i].observables[2]);
				break;
				
			case 4:
				sscanf(line, "%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf", &frame->atoms[i].id, \
				&frame->atoms[i].mol, &frame->atoms[i].type, &frame->atoms[i].q, &frame->atoms[i].mass, \ 
				&frame->atoms[i].x, &frame->atoms[i].y, &frame->atoms[i].z, &frame->atoms[i].observables[0], \ 
				&frame->atoms[i].observables[1], &frame->atoms[i].observables[2], &frame->atoms[i].observables[3]);
				break;
				
			case 5:
				sscanf(line, "%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &frame->atoms[i].id, \
				&frame->atoms[i].mol, &frame->atoms[i].type, &frame->atoms[i].q, &frame->atoms[i].mass, \ 
				&frame->atoms[i].x, &frame->atoms[i].y, &frame->atoms[i].z, &frame->atoms[i].observables[0], \ 
				&frame->atoms[i].observables[1], &frame->atoms[i].observables[2], &frame->atoms[i].observables[3], \
				&frame->atoms[i].observables[4]);
				break;
				
			case 6:
				sscanf(line, "%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &frame->atoms[i].id, \
				&frame->atoms[i].mol, &frame->atoms[i].type, &frame->atoms[i].q, &frame->atoms[i].mass, \ 
				&frame->atoms[i].x, &frame->atoms[i].y, &frame->atoms[i].z, &frame->atoms[i].observables[0], \ 
				&frame->atoms[i].observables[1], &frame->atoms[i].observables[2], &frame->atoms[i].observables[3], \
				&frame->atoms[i].observables[4], &frame->atoms[i].observables[5]);
				break;
			
			case 7:
				sscanf(line, "%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &frame->atoms[i].id, \
				&frame->atoms[i].mol, &frame->atoms[i].type, &frame->atoms[i].q, &frame->atoms[i].mass, \ 
				&frame->atoms[i].x, &frame->atoms[i].y, &frame->atoms[i].z, &frame->atoms[i].observables[0], \ 
				&frame->atoms[i].observables[1], &frame->atoms[i].observables[2], &frame->atoms[i].observables[3], \
				&frame->atoms[i].observables[4], &frame->atoms[i].observables[5], &frame->atoms[i].observables[6]);
				break;
			
			case 8:
				sscanf(line, "%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &frame->atoms[i].id, \
				&frame->atoms[i].mol, &frame->atoms[i].type, &frame->atoms[i].q, &frame->atoms[i].mass, \ 
				&frame->atoms[i].x, &frame->atoms[i].y, &frame->atoms[i].z, &frame->atoms[i].observables[0], \ 
				&frame->atoms[i].observables[1], &frame->atoms[i].observables[2], &frame->atoms[i].observables[3], \
				&frame->atoms[i].observables[4], &frame->atoms[i].observables[5], &frame->atoms[i].observables[6], \
				&frame->atoms[i].observables[7]);
				break;
				
			case 9:
				sscanf(line, "%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &frame->atoms[i].id, \
				&frame->atoms[i].mol, &frame->atoms[i].type, &frame->atoms[i].q, &frame->atoms[i].mass, \ 
				&frame->atoms[i].x, &frame->atoms[i].y, &frame->atoms[i].z, &frame->atoms[i].observables[0], \ 
				&frame->atoms[i].observables[1], &frame->atoms[i].observables[2], &frame->atoms[i].observables[3], \
				&frame->atoms[i].observables[4], &frame->atoms[i].observables[5], &frame->atoms[i].observables[6], \
				&frame->atoms[i].observables[7], &frame->atoms[i].observables[8]);
				break;
			
			default:
				printf("number of observables requested %d is not supported \n", frame->num_observables);
				break;
		}
		printf("frame->atoms[%d].mol %d and frame->num_mol %d\n", i, frame->atoms[i].mol, frame->num_mol);
		if(frame->atoms[i].mol > frame->num_mol)	frame->num_mol = frame->atoms[i].mol;
	}
	//finished reading frame
}

void process_frame(Controller* control, Frame* inframe, Frame* outframe)
{
	int i, j;
	int mol_count = 0;
	int mol_val;
	int key[inframe->num_mol];
	int site_count;
	double outx, outy, outz, tot_mass;
	double dist, box;
	outframe->type_count = 0;
	
	printf("processing frame %d\n", control->frame);
	
	if(control->frame == 1)
		{
		outframe->type = malloc(control->num_cg_types * sizeof(int));
		outframe->type_num = malloc(control->num_cg_types * sizeof(int));
		}
		
	for(i = 0; i< control->num_cg_types; i++)
		{
		outframe->type[i] = -1;
		outframe->type_num[i] = 0;
		}
		
	printf("start process frame\n");
	
	//copy basic information from inframe to outframe
	outframe->xmin = inframe->xmin;
	outframe->xmax = inframe->xmax;
	outframe->ymin = inframe->ymin;
	outframe->ymax = inframe->ymax;
	outframe->zmin = inframe->zmin;
	outframe->zmax = inframe->zmax;
	
	outframe->timestep = inframe->timestep;
	outframe->num_observables = inframe->num_observables;
	
	//determine if we need to allocate sites
	if(control->num_cg_sites != outframe->num_atoms) //assumes unitialized (0) although flexibility could be added to change number of sites (e.g. Grand Canonical Ensemble)
		{
		printf("allocate sites for outframe to size %d\n", control->num_cg_sites);
		outframe->num_atoms = control->num_cg_sites;
		outframe->sites = malloc(outframe->num_atoms * sizeof(SITE));
		
		//also need to allocate observables
		for(i=0; i < outframe->num_atoms; i++)
			{
			outframe->sites[i].observables = malloc(outframe->num_observables * sizeof(double));
			outframe->sites[i].coord = malloc(control->max_to_map *sizeof(COORD));
			}
		printf("observables and coord space allocated\n");
		}
	
	//intialize key
	for(i = 0; i < inframe->num_mol; i++)
		{
		key[i] = -1;
		}
	for(i = 0; i < control->num_cg_types; i++)
		{
		outframe->type[i] = -1;
		}

	//reset frame info
	for(i = 0; i < outframe->num_atoms; i++)
		{
		outframe->sites[i].q = 0.0;
		outframe->sites[i].num_in_site = 0;
		
		for(j = 0; j < outframe->num_observables; j++)
			{
			outframe->sites[i].observables[j] = 0.0;
			}
		
		for(j = 0; j < control->max_to_map; j++)
			{
			outframe->sites[i].coord[j].x = 0.0;
			outframe->sites[i].coord[j].y = 0.0;
			outframe->sites[i].coord[j].z = 0.0;
			}
		}
	
	printf("accumulate information from fine_grained sites as we go for frame at timestep %d\n", inframe->timestep);
	//accumulate information from fine_grained sites as we go (assuming molecule to molecule mapping)
	if(control->map_style_flag == 0) //map all atoms
		{
		
		if(control->observable_map_flag == 0) //sum obesrvables
			{
			//process and sort all FG atoms
			for(i = 0; i < inframe->num_atoms; i++)
				{
				//check to see if mol key is set
				mol_val = inframe->atoms[i].mol - 1;
				//printf("key[%d] = %d\n", mol_val, key[mol_val]);
				if(key[mol_val] == -1)
					{
					printf("set information on key[%d] to %d\n", mol_val, mol_count);
					key[mol_val] = mol_count;
					outframe->sites[key[mol_val]].id = key[mol_val] + 1;
					outframe->sites[key[mol_val]].mol = key[mol_val] + 1;
					
					//also set molecule type by looking to match first id
					outframe->sites[key[mol_val]].type = -1;
					for(j = 0; j < outframe->type_count; j++)
						{
						if( outframe->type[j] == inframe->atoms[i].type) outframe->sites[key[mol_val]].type = j + 1;
						}
					//assign type if not set before
					if(outframe->sites[key[mol_val]].type == -1)
						{
						outframe->sites[key[mol_val]].type = outframe->type_count + 1;
						outframe->type[outframe->type_count] = inframe->atoms[i].type;
						outframe->type_count++;
						printf("TYPE_COUNT IS %d ON TYPE %d\n", outframe->type_count, inframe->atoms[i].type);
						}
					mol_count++;
					outframe->type_num[ outframe->sites[key[mol_val]].type - 1]++;
					}					
				
				//printf("transfer information to CG site\n");
				//transfer information to CG site
				site_count = outframe->sites[ key[mol_val] ].num_in_site;
				//printf("site_count is %d\n", site_count);
				 
				outframe->sites[key[mol_val]].coord[site_count].x = inframe->atoms[i].x;
				outframe->sites[key[mol_val]].coord[site_count].y = inframe->atoms[i].y;
				outframe->sites[key[mol_val]].coord[site_count].z = inframe->atoms[i].z;
				outframe->sites[key[mol_val]].coord[site_count].mass = inframe->atoms[i].mass;
				outframe->sites[key[mol_val]].q += inframe->atoms[i].q;
				//printf("outframe->sites[key[%d]].q = %lf\n", mol_val, outframe->sites[key[mol_val]].q);
				
				for(j = 0; j < outframe->num_observables; j++)
					{
					outframe->sites[key[mol_val]].observables[j] += inframe->atoms[i].observables[j];
					//printf("observable[%d] is %lf\n", j, outframe->sites[key[mol_val]].observables[j]);
					}
				printf("finished observable transfer at num_in_site %d\n", outframe->sites[ key[mol_val] ].num_in_site);
				outframe->sites[ key[mol_val] ].num_in_site++;
				}
					
			//apply averaging and processing
			if(control->geometry_map_flag == 0) //map to CoM
				{
				for(i = 0; i < outframe->num_atoms; i++)
					{
					outx = 0.0;
					outy = 0.0;
					outz = 0.0;
					tot_mass = 0.0;
					
					for(j = 0; j < outframe->sites[i].num_in_site; j++)
						{
						//check to see if coordinates are wrapped (and reset out* values)
						dist = outframe->sites[i].coord[j].x - outframe->sites[i].coord[0].x;
						box = outframe->xmax - outframe->xmin;
						if( abs(dist) >= (box/2.0) )
							{
							if(dist > 0) outframe->sites[i].coord[j].x -= box;
							else		 outframe->sites[i].coord[j].x += box;
							}
						
						
						dist = outframe->sites[i].coord[j].y - outframe->sites[i].coord[0].y;
						box = outframe->ymax - outframe->ymin;
						if( abs(dist) >= (box/2.0) )
							{
							if(dist > 0) outframe->sites[i].coord[j].y -= box;
							else		 outframe->sites[i].coord[j].y += box;
							}
						
						
						dist = outframe->sites[i].coord[j].z - outframe->sites[i].coord[0].z;
						box = outframe->zmax - outframe->zmin;
						if( abs(dist) >= (box/2.0) )
							{
							if(dist > 0) outframe->sites[i].coord[j].z -= box;
							else		 outframe->sites[i].coord[j].z += box;
							}
							
						//add value to weighted sum
						outx += outframe->sites[i].coord[j].x * outframe->sites[i].coord[j].mass;
						outy += outframe->sites[i].coord[j].y * outframe->sites[i].coord[j].mass;
						outz += outframe->sites[i].coord[j].z * outframe->sites[i].coord[j].mass;
						tot_mass += outframe->sites[i].coord[j].mass;
						}
					
					//calculate average positions
					outframe->sites[i].x = outx / tot_mass;
					outframe->sites[i].y = outy / tot_mass;
					outframe->sites[i].z = outz / tot_mass;
					printf("done finding center of mass for molecules %d \n", i);

					//check if final coordinate is in box (wrap)
					if( outframe->sites[i].x > outframe->xmax ) outframe->sites[i].x -= outframe->xmax - outframe->xmin;
					if( outframe->sites[i].x < outframe->xmin ) outframe->sites[i].x += outframe->xmax - outframe->xmin;
					if( outframe->sites[i].y > outframe->ymax ) outframe->sites[i].y -= outframe->ymax - outframe->ymin;
					if( outframe->sites[i].y < outframe->ymin ) outframe->sites[i].y += outframe->ymax - outframe->ymin;
					if( outframe->sites[i].z > outframe->zmax ) outframe->sites[i].z -= outframe->zmax - outframe->zmin;
					if( outframe->sites[i].z < outframe->zmin ) outframe->sites[i].z += outframe->zmax - outframe->zmin;
						
					//set total mass
					outframe->sites[i].mass = tot_mass;
					}
								
				}
			else if(control->geometry_map_flag == 1) //map to CoG
				{
				for(i = 0; i < outframe->num_atoms; i++)
					{
					outx = 0.0;
					outy = 0.0;
					outz = 0.0;
					tot_mass = 0.0;
					
					for(j = 0; j < outframe->sites[i].num_in_site; j++)
						{
						//check to see if coordinates are wrapped (and reset out* values)
						dist = outframe->sites[i].coord[j].x - outframe->sites[i].coord[0].x;
						box = outframe->xmax - outframe->xmin;
						if( abs(dist) >= (box/2.0) )
							{
							if(dist > 0) outframe->sites[i].coord[j].x -= box;
							else		 outframe->sites[i].coord[j].x += box;
							}
						
						
						dist = outframe->sites[i].coord[j].y - outframe->sites[i].coord[0].y;
						box = outframe->ymax - outframe->ymin;
						if( abs(dist) >= (box/2.0) )
							{
							if(dist > 0) outframe->sites[i].coord[j].y -= box;
							else		 outframe->sites[i].coord[j].y += box;
							}
						
						dist = outframe->sites[i].coord[j].z - outframe->sites[i].coord[0].z;
						box = outframe->zmax - outframe->zmin;
						if( abs(dist) >= (box/2.0) )
							{
							if(dist > 0) outframe->sites[i].coord[j].z -= box;
							else		 outframe->sites[i].coord[j].z += box;
							}
							
						//add value to weighted sum
						outx += outframe->sites[i].coord[j].x;
						outy += outframe->sites[i].coord[j].y;
						outz += outframe->sites[i].coord[j].z;
						tot_mass += outframe->sites[i].coord[j].mass;
						}
					
					//calculate average positions
					outframe->sites[i].x = outx / outframe->sites[i].num_in_site;
					outframe->sites[i].y = outy / outframe->sites[i].num_in_site;
					outframe->sites[i].z = outz / outframe->sites[i].num_in_site;
					
					//check if final coordinate is in box (wrap)
					if( outframe->sites[i].x > outframe->xmax ) outframe->sites[i].x -= outframe->xmax - outframe->xmin;
					if( outframe->sites[i].x < outframe->xmin ) outframe->sites[i].x += outframe->xmax - outframe->xmin;
					if( outframe->sites[i].y > outframe->ymax ) outframe->sites[i].y -= outframe->ymax - outframe->ymin;
					if( outframe->sites[i].y < outframe->ymin ) outframe->sites[i].y += outframe->ymax - outframe->ymin;
					if( outframe->sites[i].z > outframe->zmax ) outframe->sites[i].z -= outframe->zmax - outframe->zmin;
					if( outframe->sites[i].z < outframe->zmin ) outframe->sites[i].z += outframe->zmax - outframe->zmin;
						
					//set total mass
					outframe->sites[i].mass = tot_mass;
					}
				}		
			}
		}
		
	else if(control->map_style_flag == 1) //map only select atoms
		{
		int go_flag;
		
		if(control->observable_map_flag == 0) //sum obesrvables
			{
			//process and sort all FG atoms
			for(i = 0; i < inframe->num_atoms; i++)
				{
				
				//see if molecule is in include list
				go_flag = 0;
				for(j = 0; j < control->num_map; j++)
					{
					if(control->map[j] == inframe->atoms[i].id)
						{
						go_flag = 1;
						break;
						}
					}
				
				if(go_flag == 0) continue;
				
				//check to see if mol key is set
				if(key[i] == -1)
					{
					key[i] = mol_count;
					outframe->sites[key[i]].id = key[i];
					outframe->sites[key[i]].mol = key[i];
					mol_count++;
					}					
					
				//transfer information to CG site
				site_count = outframe->sites[ key[i] ].num_in_site;
				
				//how to set type??
				outframe->sites[key[i]].coord[site_count].x = inframe->atoms[i].x;
				outframe->sites[key[i]].coord[site_count].y = inframe->atoms[i].y;
				outframe->sites[key[i]].coord[site_count].z = inframe->atoms[i].z;
				outframe->sites[key[i]].coord[site_count].mass = inframe->atoms[i].mass;
				outframe->sites[key[i]].q += inframe->atoms[i].q;
				
				for(j = 0; j < outframe->num_observables; j++)
					{
					outframe->sites[key[i]].observables[j] += inframe->atoms[i].observables[j];
					}
				
				outframe->sites[ key[i] ].num_in_site++;
				}
					
			//apply averaging and processing
			if(control->geometry_map_flag == 0) //map to CoM
				{
				for(i = 0; i < outframe->num_atoms; i++)
					{
					outx = 0.0;
					outy = 0.0;
					outz = 0.0;
					tot_mass = 0.0;
					
					for(j = 0; j < outframe->sites[i].num_in_site; j++)
						{
						//check to see if coordinates are wrapped (and reset out* values)
						dist = outframe->sites[i].coord[j].x - outframe->sites[i].coord[0].x;
						box = outframe->xmax - outframe->xmin;
						if( abs(dist) >= (box/2.0) )
							{
							if(dist > 0) outframe->sites[i].coord[j].x -= box;
							else		 outframe->sites[i].coord[j].x += box;
							}
						
						
						dist = outframe->sites[i].coord[j].y - outframe->sites[i].coord[0].y;
						box = outframe->ymax - outframe->ymin;
						if( abs(dist) >= (box/2.0) )
							{
							if(dist > 0) outframe->sites[i].coord[j].y -= box;
							else		 outframe->sites[i].coord[j].y += box;
							}
						
						
						dist = outframe->sites[i].coord[j].z - outframe->sites[i].coord[0].z;
						box = outframe->zmax - outframe->zmin;
						if( abs(dist) >= (box/2.0) )
							{
							if(dist > 0) outframe->sites[i].coord[j].z -= box;
							else		 outframe->sites[i].coord[j].z += box;
							}
							
						//add value to weighted sum
						outx += outframe->sites[i].coord[j].x * outframe->sites[i].coord[j].mass;
						outy += outframe->sites[i].coord[j].y * outframe->sites[i].coord[j].mass;
						outz += outframe->sites[i].coord[j].z * outframe->sites[i].coord[j].mass;
						tot_mass += outframe->sites[i].coord[j].mass;
						}
					
					//calculate average positions
					outframe->sites[i].x = outx / tot_mass;
					outframe->sites[i].y = outy / tot_mass;
					outframe->sites[i].z = outz / tot_mass;
					
					//check if final coordinate is in box (wrap)
					if( outframe->sites[i].x > outframe->xmax ) outframe->sites[i].x -= outframe->xmax - outframe->xmin;
					if( outframe->sites[i].x < outframe->xmin ) outframe->sites[i].x += outframe->xmax - outframe->xmin;
					if( outframe->sites[i].y > outframe->ymax ) outframe->sites[i].y -= outframe->ymax - outframe->ymin;
					if( outframe->sites[i].y < outframe->ymin ) outframe->sites[i].y += outframe->ymax - outframe->ymin;
					if( outframe->sites[i].z > outframe->zmax ) outframe->sites[i].z -= outframe->zmax - outframe->zmin;
					if( outframe->sites[i].z < outframe->zmin ) outframe->sites[i].z += outframe->zmax - outframe->zmin;
						
					//set total mass
					outframe->sites[i].mass = tot_mass;
					}
				
				}
			else if(control->geometry_map_flag == 1) //map to CoG
				{
				for(i = 0; i < outframe->num_atoms; i++)
					{
					outx = 0.0;
					outy = 0.0;
					outz = 0.0;
					tot_mass = 0.0;
					
					for(j = 0; j < outframe->sites[i].num_in_site; j++)
						{
						//check to see if coordinates are wrapped (and reset out* values)
						dist = outframe->sites[i].coord[j].x - outframe->sites[i].coord[0].x;
						box = outframe->xmax - outframe->xmin;
						if( abs(dist) >= (box/2.0) )
							{
							if(dist > 0) outframe->sites[i].coord[j].x -= box;
							else		 outframe->sites[i].coord[j].x += box;
							}
						
						
						dist = outframe->sites[i].coord[j].y - outframe->sites[i].coord[0].y;
						box = outframe->ymax - outframe->ymin;
						if( abs(dist) >= (box/2.0) )
							{
							if(dist > 0) outframe->sites[i].coord[j].y -= box;
							else		 outframe->sites[i].coord[j].y += box;
							}
						
						
						dist = outframe->sites[i].coord[j].z - outframe->sites[i].coord[0].z;
						box = outframe->zmax - outframe->zmin;
						if( abs(dist) >= (box/2.0) )
							{
							if(dist > 0) outframe->sites[i].coord[j].z -= box;
							else		 outframe->sites[i].coord[j].z += box;
							}
							
						//add value to weighted sum
						outx += outframe->sites[i].coord[j].x;
						outy += outframe->sites[i].coord[j].y;
						outz += outframe->sites[i].coord[j].z;
						tot_mass += outframe->sites[i].coord[j].mass;
						}
					
					//calculate average positions
					outframe->sites[i].x = outx / outframe->sites[i].num_in_site;
					outframe->sites[i].y = outy / outframe->sites[i].num_in_site;
					outframe->sites[i].z = outz / outframe->sites[i].num_in_site;
					
					//check if final coordinate is in box (wrap)
					if( outframe->sites[i].x > outframe->xmax ) outframe->sites[i].x -= outframe->xmax - outframe->xmin;
					if( outframe->sites[i].x < outframe->xmin ) outframe->sites[i].x += outframe->xmax - outframe->xmin;
					if( outframe->sites[i].y > outframe->ymax ) outframe->sites[i].y -= outframe->ymax - outframe->ymin;
					if( outframe->sites[i].y < outframe->ymin ) outframe->sites[i].y += outframe->ymax - outframe->ymin;
					if( outframe->sites[i].z > outframe->zmax ) outframe->sites[i].z -= outframe->zmax - outframe->zmin;
					if( outframe->sites[i].z < outframe->zmin ) outframe->sites[i].z += outframe->zmax - outframe->zmin;
						
					//set total mass
					outframe->sites[i].mass = tot_mass;
					}
				}		
			}
		}
	outframe->num_mol = mol_count;		
}

void output_frame(Controller* control, Frame* outframe, char* outfile)
{
	int i, j;
	FILE* of = fopen(outfile, "a+");
	
	//output frame header
	fprintf(of, "ITEM: TIMESTEP\n");
	fprintf(of, "%d\n", outframe->timestep);
	fprintf(of, "ITEM: NUMBER OF ATOMS\n");
	fprintf(of, "%d\n", outframe->num_atoms); //could also consider using outframe->num_mol
	fprintf(of, "ITEM: BOX BOUNDS pp pp pp\n");
	fprintf(of, "%lf %lf\n%lf %lf\n %lf, %lf\n", outframe->xmin, outframe->xmax, outframe->ymin, outframe->ymax, outframe->zmin, outframe->zmax);

	switch(outframe->num_observables)
		{	
		case 0:
			fprintf(of, "ITEM: ATOMS id mol type q mass x y z \n");
			break;
		
		case 1:
			fprintf(of, "ITEM: ATOMS id mol type q mass x y z U \n");
			break;
			
		case 3:
			fprintf(of, "ITEM: ATOMS id mol type q mass x y z fx fy fz \n");
			break;
		
		case 4:
			fprintf(of, "ITEM: ATOMS id mol type q mass x y z dfx dfy dfz dU \n");
			break;
			
		case 6:
			fprintf(of, "ITEM: ATOMS id mol type q mass x y z fx fy fz dfx dfy dfz \n");
			break;
			
		case 7:
			fprintf(of, "ITEM: ATOMS id mol type q mass x y z fx fy fz dfx dfy dfz dU \n");
			break;
			
		default:
			printf("output number of observables not supported so giving generic 3 observable header\n");
			fprintf(of, "ITEM: ATOMS id mol type q mass x y z fx fy fz dfx dfy dfz dU \n");
			//for(i = 0; i < outframe->num_atoms; i++)
		}
		
	for( i = 0; i < outframe->num_atoms; i++)
		{
		switch(outframe->num_observables)
			{
			case 0:
				fprintf(of, "%d %d %d %lf %lf %lf %lf %lf\n", outframe->sites[i].id, outframe->sites[i].mol, \
				outframe->sites[i].type, outframe->sites[i].q, outframe->sites[i].mass, outframe->sites[i].x, \
				outframe->sites[i].y, outframe->sites[i].z);
				break;
			
			case 1:
				fprintf(of, "%d %d %d %lf %lf %lf %lf %lf %lf\n", outframe->sites[i].id, outframe->sites[i].mol, \
				outframe->sites[i].type, outframe->sites[i].q, outframe->sites[i].mass, outframe->sites[i].x, \
				outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0]);
				break;
				
			case 2:
				fprintf(of, "%d %d %d %lf %lf %lf %lf %lf %lf %lf\n", outframe->sites[i].id, outframe->sites[i].mol, \
				outframe->sites[i].type, outframe->sites[i].q, outframe->sites[i].mass, outframe->sites[i].x, \
				outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0], \
				outframe->sites[i].observables[1]);
				break;

			case 3:
				fprintf(of, "%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf\n", outframe->sites[i].id, outframe->sites[i].mol, \
				outframe->sites[i].type, outframe->sites[i].q, outframe->sites[i].mass, outframe->sites[i].x, \
				outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0], \
				outframe->sites[i].observables[1], outframe->sites[i].observables[2]);
				break;
				
			case 4:
				fprintf(of, "%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", outframe->sites[i].id, outframe->sites[i].mol, \
				outframe->sites[i].type, outframe->sites[i].q, outframe->sites[i].mass, outframe->sites[i].x, \
				outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0], \
				outframe->sites[i].observables[1], outframe->sites[i].observables[2], \
				outframe->sites[i].observables[3]);
				break;
			
			case 5:
				fprintf(of, "%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", outframe->sites[i].id, outframe->sites[i].mol, \
				outframe->sites[i].type, outframe->sites[i].q, outframe->sites[i].mass, outframe->sites[i].x, \
				outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0], \
				outframe->sites[i].observables[1], outframe->sites[i].observables[2], \
				outframe->sites[i].observables[3], outframe->sites[i].observables[4]);
				break;
				
			case 6:
				fprintf(of, "%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", outframe->sites[i].id, outframe->sites[i].mol, \
				outframe->sites[i].type, outframe->sites[i].q, outframe->sites[i].mass, outframe->sites[i].x, \
				outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0], \
				outframe->sites[i].observables[1], outframe->sites[i].observables[2], \
				outframe->sites[i].observables[3], outframe->sites[i].observables[4], \
				outframe->sites[i].observables[5]);
				break;
				
			case 7:
				fprintf(of, "%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", outframe->sites[i].id, outframe->sites[i].mol, \
				outframe->sites[i].type, outframe->sites[i].q, outframe->sites[i].mass, outframe->sites[i].x, \
				outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0], \
				outframe->sites[i].observables[1], outframe->sites[i].observables[2], \
				outframe->sites[i].observables[3], outframe->sites[i].observables[4], \
				outframe->sites[i].observables[5], outframe->sites[i].observables[6]);
				break;
				
			case 8:
				fprintf(of, "%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", outframe->sites[i].id, outframe->sites[i].mol, \
				outframe->sites[i].type, outframe->sites[i].q, outframe->sites[i].mass, outframe->sites[i].x, \
				outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0], \
				outframe->sites[i].observables[1], outframe->sites[i].observables[2], \
				outframe->sites[i].observables[3], outframe->sites[i].observables[4], \
				outframe->sites[i].observables[5], outframe->sites[i].observables[6], \
				outframe->sites[i].observables[7]);
				break;
				
			case 9:
				fprintf(of, "%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", outframe->sites[i].id, outframe->sites[i].mol, \
				outframe->sites[i].type, outframe->sites[i].q, outframe->sites[i].mass, outframe->sites[i].x, \
				outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0], \
				outframe->sites[i].observables[1], outframe->sites[i].observables[2], \
				outframe->sites[i].observables[3], outframe->sites[i].observables[4], \
				outframe->sites[i].observables[5], outframe->sites[i].observables[6], \
				outframe->sites[i].observables[7], outframe->sites[i].observables[8]);
				break;
			
			default:
				printf("number of observables to output is not supported in file writing\n");
				break;
			}
		}
		
	fclose(of);
}

void output_topology(Controller* control, Frame* outframe)
{
	FILE* fp = fopen("top.in", "w+");
	int i, j;
	
	fprintf(fp, "cgsites %d\n", outframe->num_atoms);
	
	fprintf(fp, "cgtypes %d\n", outframe->type_count);
	for(i = 0; i < outframe->type_count; i++)
		{
		fprintf(fp, "%d\n", i+1);
		}
	
	fprintf(fp, "moltypes %d\n", outframe->type_count);
	for(i = 0; i < outframe->type_count; i++)
		{
		fprintf(fp, "mol %d %d\n", 1, 3); 
		}
	
	fprintf(fp, "sitetypes\n");
	for(i = 0; i < 1; i++)
		{
		fprintf(fp, "%d\n", i+1);
		}
	
	fprintf(fp, "bonds %d\n", 0);
	//for(i = 0; i < 0; i++)
	//	{
	//	fprintf(fp, "%d %d\n", 0, 0); 
	//	}
	
	fprintf(fp, "system %d\n", outframe->type_count); 
	for(i = 0; i < outframe->type_count; i++)
		{
		fprintf(fp, "%d %d\n", i+1, outframe->type_num[i]);
		}
	fprintf(fp, "\n");
	
	fclose(fp);
}


void check_file_extension(char* name, char* suffix)
{
    char temp;
    int len, pos, i;
    len = strlen(name);
    pos = -1;
    for (i = 0; i < len; i++) {
        if (name[i] == '.') pos = i;
    }
    if (pos < 0) report_traj_input_suffix_error(suffix);
    if (strcmp(&name[pos + 1], suffix) != 0) report_traj_input_suffix_error(suffix);
}

void report_traj_input_suffix_error(char *suffix)
{
    printf("Failed to find a valid file extension (.%s) in the supplied trajectory filename.\n", suffix);
    exit(EXIT_SUCCESS);
}

void report_usage_error(char *exe_name)
{
    printf("Usage: %s -f file.dat -f1 file1.top -o out.dat with missing file flag-name pairs allowed\n", exe_name);
    exit(EXIT_SUCCESS);
}

void free_allocation(Controller* control, Frame* inframe, Frame* outframe)
{
	int i;
	
	//free all data for atoms
	for(i = 0; i < inframe->num_atoms; i++)
		{
		free(inframe->atoms[i].observables);	
		}
	free(inframe->atoms);

	for(i =	0; i < outframe->num_atoms; i++)
		{
		free(outframe->sites[i].observables);
		free(outframe->sites[i].coord);
		}
	free(outframe->sites);
	free(outframe->type);
	free(outframe->type_num);
}