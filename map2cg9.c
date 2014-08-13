//MAP2CG
//By Jacob Wagner
// Created July 8, 2013
// REV9: November 13, 2013
//Creates a data file and a topology file where a specified coarse-grained mapping
//has been applied to the frames of a custom LAMMPS dump.
//Useful in mapping an observable(s) from a fine-grained (FG) to a coarse-grained (CG) system.

//Acknowledgements: One of the outputs is simillar to the A2CG program by Zhen that uses GROMAC's libraries.

//INPUT: Custom LAMMPS DUMP (see next line) and a TOPOLOGY mapping file
//1: dump tj2 all custom 250 MeOH_sample.dat	id	mol type q mass x y z fx fy fz

//2: topology file should be of the form:
// int	#number of CG sites
// int  #number CG types
// int	#number of FG sites
// int  #number FG types
// int  #max number of FG sites in single CG site
//		#this line intentionally left blank
// int *** #number of atoms in site  for #cg types listed above
// int *** #list types for each atom in cg molecule (# should match above line's value)
// ...     #continue until #CG types is fulfilled
//		#this line intentionally left blank
// int	#geometry mapping flag (0 = CoM geometry mapping, 1 = center of geometry mapping)
// int  #observable mapping flag (0 = additive mappign)
// int	#number observables to map
// int  #output flag (0 = all, 1 = minimal, 2 = 1 value is #4 observable, 3 = 2nd sets of 3 values)
// int  #sensitivity flag (0 = only map 1 file, 1 = map 2 files along with log data and guess, 2 = convert "minimal" output to filler "all" output, 3 = convert MSCGFM (1_1.dat) to (tab_ff.dat), 4 = sort CG frame by type, 5 = bootstrapping rearrangement, 6 = II CG charge force derivative, 7 = IJ CG charge force derivative)
//	 	#this line intentionally left blank
// int  #mapping style flag (0 = 1:1 molecule entirely, 1 = implicit solvent, 2 = truncate trajectory)
//IF mapping style flag == 1
// int	#number of types to map
// int #type ID's to be mapped
// ...
//IF mapping style flag == 2
// int #number of atoms before truncate
//ENDIF
//		#this line intentionally left blank
//IF sensitivity flag == 1 (Composite FM derivative processing)
//RECOMMEND: 4-1 debug flag (or 4+1, 8-1, 8+1)
// %s filename for 1st dump (original data, regular parameters)
// %s filename for 2nd dump (rerun data, force derivative wrt parameter, e.g., epsilon_11)
// %s filename for log.lammps (rerun data, potential derivative wrt parameter, e.g., epsilon_11)
// %s filename for guess of dU/d\lambda (needs as many rows as frames in other inputs)
// %s filename for composite output (lammps dump file)
// int	#log type (0 = vdW, 1 = Columb)
// int  #guess type (0 = 1 column, 1 = log file similar log.lammps file)
// int  #sensitivity mapping flag (0 = dump files are FG, 1 = dump files are CG)
// int	#debug_flag (0 = normal, 1 = U/N_cg, 2 = U/temp, 3 = all/temp, 4=1+2, 5=1+3, 6 = N_cg and N_fg, 7 = 2 + 6, 8 = U/temp + U/3Ncg, 9 = U/temp + U/3N_fg and U/3N_cg)
// int	#sign_flag (-1 = negative, +1 = positive)
//ENDIF
//IF sensitivity flag == 2 (convert "minimal" format to "all" adding filler id, mol, type, q, and mass fields)
// #blank
// %lf 	#charge
// %lf 	#mass
//IF sensitivity flag == 3 (convert FM output to LAMMPS TABULATED potential)
// #blank
// %s name of ID for tab potential (e.g. SENS_MEOH)
// int #units flag (0 convert from nm and kj/mol to A and kcal/mol, 1 no conversion)
//ENDIF
//IF sensitivity flag == 4 (sort by type within frame -- no inputs needed)
//ENDIF
//IF sensitivity flag == 5 (bootstrapping creation)
// #blank
// int int #start number of file and end number of file
// int     #frames in trajectory
// %s  #base name of file with ending "#.dat" -- output will be "base#map.dat"
//IF sensitivity flag == 6 (II charge frames)
// #blank
// int #charge value
// %s input filename 
// %s filename for output
//		#end of input file, please left blank line at end
//IF sensitivity flag == 7 (IJ charge frames)
// #blank
// int int #charge values for each site in order of input files (II, JJ)
// %s filename for mixed charge (IJ)
// %s filename for 1st self interaction (II) 
// %s filename for 2nd self interaction (JJ)
// %s filename for output *I
// %s filename for output *J
//		#end of input file, please left blank line at end
//ENDIF

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

//structure defintions
#include "headers.h"

//external files
#include "output_lmp.c"
#include "process_frame.c"
#include "reading.c"

//TBA
void handle_cg_frame_read(Controller*, Frame*);
//function prototype definitions
void do_simple_map(Controller*, Frame*, Frame*, char*, char*);
void do_sensitivity_map(Controller*, Frame*, Frame*, char*);
void do_force_conversion(Controller*, char*, char*);
void do_bootstrapping(Controller*, char*);
void do_ii_charge_derivative(Controller*, Frame*, Frame*);
void do_ij_charge_derivative(Controller*);
void do_map_to_full(Controller*, Frame*, Frame*, char*, char*);
void do_sort_frame_by_type(Controller*, Frame*, Frame*, char*, char*);
void sensitivity_no_mapping(Controller*, FILE*, FILE*, FILE*, FILE*, FILE*, Frame*, Frame*, Frame*, char*);
void sensitivity_with_mapping(Controller*, FILE*, FILE*, FILE*, FILE*, FILE*, Frame*, Frame*, Frame*, char*);
void handle_cg_frame_read(Controller*, Frame*);
//free 
void free_allocation(Controller*);
void free_sensitivity_allocation(Controller*, Frame*, Frame*, Frame*); 
//void free_charge_dump_allocation(Controller*, Frame**, Frame**);
void free_inframes(Controller*, Frame*);
void free_outframes(Controller*, Frame*);

//////////////////////////
///      MAIN   	  ///
/////////////////////////

int main(int argc,char *argv[])
{
	//declare variables
   	double start_cputime = clock();
   	Controller controls;
   	Frame frame;
   	Frame outframe;
   	frame.num_atoms = 0;
   	frame.num_mol = 0;
   	frame.num_observables = 0;
	outframe.num_atoms = 0;
	outframe.num_mol = 0;
	outframe.num_observables = 0;
	controls.frame = 0;
	controls.sensitivity_flag = 0;
	controls.num_cg_sites = 0;
	controls.guess_type = 0;
	controls.log_type = 0;
	controls.sens_map_flag = 0;
	
	char* datfile = malloc(64 * sizeof(char));
	char* topfile = malloc(64 * sizeof(char));
	char* outfile = malloc(64 * sizeof(char));
	
	//read topology/control file
	parse_command_line_arguments(argc, argv, datfile, topfile, outfile);
 
	printf("assigned file names are %s and %s and %s\n", datfile, topfile, outfile);
	
	//read topology/control file
	printf("to read_topology_file\n");
	read_topology_file(&controls, topfile);
	printf("finished reading topology_file\n");
	
	//determine appropriate controller function
	if(controls.sensitivity_flag == 0) 
		{
		do_simple_map(&controls, &frame, &outframe, datfile, outfile);
		}
	else if(controls.sensitivity_flag == 1) 
		{
		do_sensitivity_map(&controls, &frame, &outframe, outfile);	
		}
	else if(controls.sensitivity_flag == 2) 
		{
		do_map_to_full(&controls, &frame, &outframe, datfile, outfile);
		}
	else if(controls.sensitivity_flag == 3) 
		{
		do_force_conversion(&controls, datfile, outfile);
		}
	else if(controls.sensitivity_flag == 4) 
		{
		do_sort_frame_by_type(&controls, &frame, &outframe, datfile, outfile);
		}
	else if(controls.sensitivity_flag == 5) 
		{
		do_bootstrapping(&controls, datfile);
		}
	else if(controls.sensitivity_flag == 6) 
		{
		do_ii_charge_derivative(&controls, &frame, &outframe);
		}
	else if(controls.sensitivity_flag == 7) 
		{
		do_ij_charge_derivative(&controls);
		}
	//free allocated variables
	printf("free basic files\n");
	free(datfile);
	free(topfile);
	free(outfile);
	
	printf("free map information\n");
	free(controls.map);
				
	//print out run statistics
	double final_cputime = clock();
	printf("total run time was %lf seconds\n", (final_cputime - start_cputime)/( (double) CLOCKS_PER_SEC ) );
}

//////////////////////////
///  do_simple_map 	  ///
/////////////////////////

void do_simple_map(Controller* controls, Frame* inframe, Frame* outframe, char* datfile, char* outfile)
{
	int cont_flag = 1;
	int frame_count = 0;
	FILE* df;
	FILE* of;

	//open dump file and toggle output file to reset
	df = fopen(datfile, "rt");
	of = fopen(outfile, "w+");
	fclose(of);
	
	//read first frame
	frame_count++;
	controls->frame++;
	read_frame(controls, inframe, df, &cont_flag);
	//printf("main says num_observables is %d\n", controls->num_observables);
	
	while(cont_flag == 1)
		{
		//process/map frame
		process_frame(controls, inframe, outframe);		
		//output mapped frame and observables
		output_frame(controls, outframe, outfile);
		
		//read next frame or set flag if done
		frame_count++;
		controls->frame++;
		read_frame(controls, inframe, df, &cont_flag);
		}
		
	//close dump file
	fclose(df);
	
	//create top.in file for FM
	output_topology(controls, outframe);
	
	free_inframes(controls, inframe);
	free_outframes(controls, outframe);
	free_allocation(controls);
}

//////////////////////////////////
///   do_sensitivity_map	  ///
////////////////////////////////

void do_sensitivity_map(Controller* controls, Frame* inframe1, Frame* outframe, char* outfile)
{
	int cont_flag = 1;
	int frame_count = 0;

   	Frame inframe2;
   	inframe2.num_atoms = 0;
   	inframe2.num_mol = 0;
   	inframe2.num_observables = 0;
 
	FILE* df1;
	FILE* df2;
	FILE* lf;
	FILE* gf;
	FILE* of;
	
	//open dump file and toggle output file to reset
	df1 = fopen(controls->files.dump1, "rt");
	df2 = fopen(controls->files.dump2, "rt");
	lf  = fopen(controls->files.log,   "rt");
	gf  = fopen(controls->files.guess, "rt");
	of  = fopen(outfile, "w+");
	fclose(of);
	
	//check that all files opened correctly
	if(df1 == NULL)
		{
		printf("dump file #1 specified does not exist\n");
		}
	if(df2 == NULL) 
		{
		printf("dump file #2 specified does not exist\n");
		}
	if(lf == NULL) 
		{
		printf("log file specified does not exist\n");
		}
	if(gf == NULL) 
		{
		printf("guess file specified does not exist\n");
		}

	//exit if there is an error
	if( (df1 == NULL) || (df2 == NULL) || (lf == NULL) || (gf == NULL) ) 
		{
		exit(EXIT_SUCCESS);
		}
		
	//read first frames
	frame_count++;
	controls->frame++;
	
	if( controls->sens_map_flag == 0) 
		{
		sensitivity_no_mapping(controls, df1, df2, lf, gf, of, inframe1, &inframe2, outframe, outfile);
		}
	else if( controls->sens_map_flag == 1) 
		{
		sensitivity_with_mapping(controls, df1, df2, lf, gf, of, inframe1, &inframe2, outframe, outfile);
		}
	
	//close files for frame, log, and guess reading
	fclose(df1);
	fclose(df2);
	fclose(lf);
	fclose(gf);
	
	//free allocated varaibles
	free_sensitivity_allocation(controls, inframe1, &inframe2, outframe); 
}

/////////////////////////////
///   do_map_to_full	  ///
/////////////////////////////

void do_map_to_full(Controller* control, Frame* inframe, Frame* outframe, char* infile, char* outfile)
{
	int cont_flag = 1;
	int frame_count = 0;
	int i = 0;	
	FILE* IF = fopen(infile, "rt");
	
	if(IF == NULL) 
		{
		printf("file number %d in list does not exist\n", (i+1) );
		cont_flag = 0;
		}
	
	//exit if there is an error
	if( cont_flag == 0 )  
		{
		//close files openend and free allocated space
		fclose(IF);
		exit(EXIT_SUCCESS);
		}
	
	//read first frames
	frame_count++;
	control->frame++;
	cont_flag = 1;
	read_frame_minimal(control, inframe, IF, &cont_flag);
	
	while(cont_flag == 1)
		{
		//process/map frame
		process_minimal_frame(control, inframe, outframe);
		//output mapped frame and observables
		output_frame(control, outframe, outfile);
		
		//read next frame or set flag if done
		frame_count++;
		control->frame++;
		read_frame_minimal(control, inframe, IF, &cont_flag);
		printf("cont_flag is %d\n", cont_flag);
		}

	//create top.in file for FM
	//output_topology(control, outframe);
	
	//free allocated varaibles	
	//close read files and output files
	//also, free frame content 
	fclose(IF);
	free_inframes(control, inframe);
	free_outframes(control, outframe );
	
	//free file pointer holders	
	//free prototypes
	for(i = 0; i < control->num_cg_types; i++)
		{
		free(control->prototype[i].num_list);
		}
	free(control->prototype);
}

//////////////////////////////////
///   do_force_conversion	  ///
////////////////////////////////

void do_force_conversion(Controller* control, char* infile, char* outfile)
{
	//declare variables
	double* distance = NULL;
	double* force = NULL;
	double* potential = NULL;
	int num_entries, i;

	//read file to determine number of lines and process input by allocating space
	//declare variables
	//int i;
	int num_lines = 0;
	int flag = 1;
	char line[100];
	FILE* fp = fopen(infile, "rt");
	
	//determine number of lines in file
	while(flag == 1)
		{
		//check if EOF
		if( fgets (line, 100, fp) == NULL )
			{
			flag = 0;
			printf("end of file reached in force file!\n");
			}
		else
			{
			num_lines++;
			}
		}
		
	//allocate space
	distance = malloc( num_lines * sizeof(double));
	force = malloc( num_lines * sizeof(double));
	
	//rewind and read in data
	printf("do rewind\n");
	rewind(fp);
	
	for(i = 0; i < num_lines; i++)
		{
		fgets(line, 100, fp);
		sscanf(line, "%lf %lf", &distance[i], &force[i]);
		printf("read at site %d values of distance %lf and force %lf\n", i, distance[i], force[i]);
		}
	
	//close input file
	fclose(fp);
	
	num_entries = num_lines;
	printf("number of lines %d = %d num lines\n", num_entries, num_lines);
	//read_force_file(control, infile, distance, force, &num_entries);	
	//integrate to get potential
	potential = malloc(num_entries * sizeof(double));
	
	printf("integrate\n");
	potential[num_entries - 1] = 0.0;
	for(i = num_entries - 2; i >= 0; i--)
		{		
		potential[i] = potential[i+1] + 0.5 * (distance[i + 1] - distance[i]) * (force[i] + force[i + 1]);
		}
	
	//convert ALL units
	if(control->units_flag == 0) 
		{
		printf("convert units for %d entries\n", num_entries);
		for(i = 0; i < num_entries; i++)
			{
			//distance (nm -> A)
			distance[i] *= 10.0;
		
			//potential (kj/mol -> kcal/mol)
			potential[i] /= 4.184;
			
			//force (kj / mol nm -> kcal/mol A)
			force[i] /= 41.84;
			}
		} 
	else if(control->units_flag == 1) 
		{
		printf("no unit conversion applied\n");
		}
	 
	printf("write_output\n");
	//write output
	output_force_file(control, num_entries, distance, potential, force, outfile);
	
	//free variables
	free(distance);
	free(force);
	free(potential);
	//free controller data alloacted
	for(i = 0; i < control->num_cg_types; i++)
		{
		free(control->prototype[i].num_list);
		}
	free(control->prototype);
}

///////////////////////////////
///		do_bootstrapping	///
///////////////////////////////

void do_bootstrapping(Controller* control, char* datfile)
{
	//declare variables
	int cont_flag = 1;
	int file_count = control->charge[0];
	int i = 0;
	printf("allocating space\n");
	Frame* mapped_frames = malloc(control->num_frames * sizeof(Frame));
	int* frame_order = malloc(control->num_frames * sizeof(int));
	FILE* fp = fopen(datfile, "rt");
	
	
	printf("read_and_process_bootstrapping_trajectory\n");
	//read in and process/map set of frames
	read_and_process_bootstrapping_trajectory(control, mapped_frames, fp, &cont_flag);	
	fclose(fp);
	
	//read in first bootstrapping file
	printf("read first bootstrapping file\n");
	read_bootstrapping_file(control, frame_order, &file_count);
		
	while(cont_flag == 1)
		{
		//create frame set and output results
		create_and_output_bootstrapping_file(control, mapped_frames, frame_order, &file_count);		
		file_count++;
		if(file_count > control->charge[1]) 
			{
			cont_flag = 0;
			}
		else
			{
			//read next bootstrapping file
			printf("read boostrap file #%d\n", file_count);
			read_bootstrapping_file(control, frame_order, &file_count);
			}
		}
		
	//clean-up and free variables
	printf("clean-up and free variables\n");
	free(frame_order);
	control->frame = 2;
	for(i = 0; i < control->num_frames; i++)
		{
		//free information in mapped_frames
		free_outframes(control, &mapped_frames[i] );
		}
		
	free(mapped_frames);
	free(control->charge);
    free(control->name);
    
    //free prototypes 
    for(i = 0; i < control->num_cg_types; i++)
		{
		free(control->prototype[i].num_list);
	}
	free(control->prototype);
}

///////////////////////////////////////////
///		do_ii_charge_derivative			///
///////////////////////////////////////////

void do_ii_charge_derivative(Controller* control, Frame* inframe, Frame* outframe)
{
	int cont_flag = 1;
	int frame_count = 0;
	int i = 0;	
	
	if(control->file_point[i] == NULL) 
		{
		printf("file number %d in list does not exist\n", (i+1) );
		cont_flag = 0;
		}
	
	//exit if there is an error
	if( cont_flag == 0 )  
		{
		//close files openend and free allocated space
		fclose(control->file_point[i]);
		free(control->file_point);
		fclose(control->outfile[i]);
		free(control->outfile);
		exit(EXIT_SUCCESS);
		}
	
	//read first frames
	frame_count++;
	control->frame++;
	cont_flag = 1;		
	read_frame(control, inframe, control->file_point[0], &cont_flag);

	while(cont_flag == 1)
		{
		//process/map frame
		process_ii_charge_frames(control, inframe, outframe);		
		//output mapped frame and observables
		output_charge_frames(control, outframe);
		
		//read next frame or set flag if done
		frame_count++;
		control->frame++;
		read_frame(control, inframe, control->file_point[0], &cont_flag);
		printf("cont_flag is %d\n", cont_flag);
		}

	printf("to output_topology\n");
	//create top.in file for FM
	output_topology(control, outframe);
	
	//close read files and output files
	//also, free frame content 
	for(i = 0; i < control->num_files; i++)
		{
		fclose(control->file_point[i]);
		free_inframes(control, inframe);
		}
	for(i = 0; i < control->num_charges; i++)
		{
		fclose(control->outfile[i]);
		free_outframes(control, outframe );
		}

	//free file pointer holders
	free(control->file_point);
	free(control->outfile);
	free(control->charge);
	
	//free prototypes
	for(i = 0; i < control->num_cg_types; i++)
		{
		free(control->prototype[i].num_list);
		}
	free(control->prototype);
}

///////////////////////////////////////////
///		do_ij_charge_derivative			///
///////////////////////////////////////////

void do_ij_charge_derivative(Controller* control)
{
	int cont_flag = 1;
	int frame_count = 0;
	int i = 0;

   	Frame* inframes;
   	inframes = malloc(control->num_files * sizeof(Frame));
	Frame* outframes;
	outframes = malloc(control->num_outfile * sizeof(Frame));
	
	//initalize num_mol for in/out-frames and check that all files opened correctly
	for(i = 0; i < control->num_files; i++)
		{
		inframes[i].num_mol = 0;
		inframes[i].num_atoms = 0;
		if(control->file_point[i] == NULL) 
			{
			printf("file number %d in list does not exist\n", (i+1) );
			cont_flag = 0;
			}
		}
	for(i = 0; i < control->num_outfile; i++)
		{
		outframes[i].num_mol = 0;
		outframes[i].num_atoms = 0;
		if(control->file_point[i] == NULL) 
			{
			printf("file number %d in list does not exist\n", (i+1) );
			cont_flag = 0;
			}
		}
	//exit if there is an error
	if( cont_flag == 0 )  
		{
		//close files openend and free allocated space
		for(i = 0; i < control->num_files; i++) 
			{
			fclose(control->file_point[i]);
			}
		free(control->file_point);
		for(i = 0; i < control->num_outfile; i++) 
			{
			fclose(control->outfile[i]);
			}
		free(control->outfile);
		free(inframes);
		free(outframes);
		exit(EXIT_SUCCESS);
		}
	
	//read first frames
	frame_count++;
	control->frame++;
	cont_flag = 1;
	read_charge_frames(control, inframes, &cont_flag);
	
	while(cont_flag == 1)
		{
		//process/map frame
		process_ij_charge_frames(control, inframes, outframes);		
		//output mapped frame and observables
		output_charge_frames(control, outframes);
		
		//read next frame or set flag if done
		frame_count++;
		control->frame++;
		read_charge_frames(control, inframes, &cont_flag);
		}
	
	printf("to output_topology\n");
	//create top.in file for FM
	output_topology(control, &outframes[0] );
	
	//free allocated variables	
	//close read files and output files
	//also, free frame content 
	for(i = 0; i < control->num_files; i++)
		{
		fclose(control->file_point[i]);
		free_inframes(control, &inframes[i] );
		}
	for(i = 0; i < control->num_charges; i++)
		{
		fclose(control->outfile[i]);
		free_outframes(control, &outframes[i] );
		}
	//free frame holders
	free(inframes);
	free(outframes);
	
	//free file pointer holders
	free(control->file_point);
	free(control->outfile);
	free(control->charge);
	
	//free prototypes
	for(i = 0; i < control->num_cg_types; i++)
		{
		free(control->prototype[i].num_list);
		}
	free(control->prototype);
}

///////////////////////////////////////
///		do_sort_frame_by_type		///
///////////////////////////////////////

void do_sort_frame_by_type(Controller* control, Frame* inframe, Frame* outframe, char* infile, char* outfile)
{
	int cont_flag = 1;
	int frame_count = 0;
	int i = 0;	
	FILE* IF = fopen(infile, "rt");
	
	if(IF == NULL) 
		{
		printf("file number %d in list does not exist\n", (i+1) );
		cont_flag = 0;
		}
	
	//exit if there is an error
	if( cont_flag == 0 )  
		{
		fclose(IF);
		exit(EXIT_SUCCESS);
		}
	
	//read first frames
	frame_count++;
	control->frame++;
	cont_flag = 1;
	read_frame(control, inframe, IF, &cont_flag);
	
	while(cont_flag == 1)
		{
		//process/map frame
		process_frame_order(control, inframe, outframe);		
		//output mapped frame and observables
		output_frame(control, outframe, outfile);
		
		//read next frame or set flag if done
		frame_count++;
		control->frame++;
		read_frame(control, inframe, IF, &cont_flag);
		printf("cont_flag is %d\n", cont_flag);
		}

	//create top.in file for FM
	output_topology(control, outframe);
	
	//free allocated variables
	//close read files and output files
	fclose(IF);
	for(i = 0; i < control->num_cg_types; i++)
		{
		free(inframe->type_list[i]);
		}
	free(inframe->type_list);
	free_inframes(control, inframe);
	free_outframes(control, outframe );
	
	//free prototypes
	for(i = 0; i < control->num_cg_types; i++)
		{
		free(control->prototype[i].num_list);
		}
	free(control->prototype);
	free(control->order);
}

///////////////////////////////////////////
///		Sensitivity_Slave_Functions		///
///////////////////////////////////////////

///////////////////////////////////////
///		sensitivity_no_mapping		///
///////////////////////////////////////

void sensitivity_no_mapping(Controller* control, FILE* df1, FILE* df2, FILE* lf, FILE* gf, FILE* of, Frame* inframe1, Frame* inframe2,  Frame* outframe, char* outfile)
{
	//declare variables 
	int cont_flag = 1;
	int frame_count = control->frame;
	
	printf("reading first frame data for sensitivity map\n");
	read_frames_and_log(control, inframe1, inframe2, df1, df2, lf, gf, &cont_flag);
	printf("1st reading reports log value is %lf and guess is %lf\n", control->log_value, control->guess);
	
	while(cont_flag == 1)
		{
		//process/map frame
		process_no_map_frames_and_log(control, inframe1, inframe2, outframe);
	
		if(control->frame == 1) 
			{
			printf("to output_topology\n");
			//create top.in file for FM
			output_topology(control, outframe);
			}
			
		//output mapped frame and observables
		output_frame(control, outframe, outfile);
		
		//read next frame or set flag if done
		frame_count++;
		control->frame++;
		read_frames_and_log(control, inframe1, inframe2, df1, df2, lf, gf, &cont_flag);
		}
}

///////////////////////////////////////////
///		sensitivity_with_mapping		///
///////////////////////////////////////////

void sensitivity_with_mapping(Controller* control, FILE* df1, FILE* df2, FILE* lf, FILE* gf, FILE* of, Frame* inframe1, Frame* inframe2, Frame* outframe, char* outfile)
{
	//declare variables 
	int cont_flag = 1;
	int frame_count = control->frame;

	printf("reading first frame data for sensitivity map\n");
	read_frames_and_log(control, inframe1, inframe2, df1, df2, lf, gf, &cont_flag);
	
	printf("main says num_observables is %d\n", control->num_observables);
	printf("1st reading reports log value is %lf and guess is %lf\n", control->log_value, control->guess);
	
	while(cont_flag == 1)
		{
		//process/map frame
		process_frames_and_log(control, inframe1, inframe2, outframe);
	
		if(control->frame == 1) 
			{
			printf("to output_topology\n");
			//create top.in file for FM
			output_topology(control, outframe);
			}
			
		//output mapped frame and observables
		output_frame(control, outframe, outfile);
		
		//read next frame or set flag if done
		frame_count++;
		control->frame++;
		read_frames_and_log(control, inframe1, inframe2, df1, df2, lf, gf, &cont_flag);
		}
}

////////////////////////////////////////
///		handle_cg_frame_read		///
//////////////////////////////////////

void handle_cg_frame_read(Controller* control, Frame* inframe)
{

}

///////////////////////////////
///		free_allocation		///
///////////////////////////////

void free_allocation(Controller* control)
{
	int i;
	//free controller data alloacted
	for(i = 0; i < control->num_cg_types; i++)
		{
		free(control->prototype[i].num_list);
		}
	free(control->prototype);
}

void free_sensitivity_allocation(Controller* control, Frame* inframe1, Frame* inframe2, Frame* outframe)
{
	int i;
	
	//free controller data alloacted
	for(i = 0; i < control->num_cg_types; i++)
		{
		free(control->prototype[i].num_list);
		}
	free(control->prototype);

	//free inframe information
	for(i = 0; i < inframe1->num_atoms; i++)
		{
		free(inframe1->atoms[i].observables);
		}
	for(i = 0; i < inframe2->num_atoms; i++)
		{
		free(inframe2->atoms[i].observables);
		}
	free(inframe1->atoms);
	free(inframe2->atoms);
	
	//free all data for atoms
	for(i = 0; i < outframe->num_atoms; i++)
		{
		free(outframe->sites[i].observables);
		}
	
	printf("inframe1->num_atoms = %d and inframe2->num_atoms = %d\n", inframe1->num_atoms, inframe2->num_atoms);

	free(outframe->sites);
	free(outframe->type);
	free(outframe->type_num);
}

void free_inframes(Controller* control, Frame* inframe)
{
	int i;
	for(i = 0; i < inframe->num_atoms; i++)
		{
		free(inframe->atoms[i].observables);	
		}
	free(inframe->atoms);
}
	
void free_outframes(Controller* control, Frame* outframe)
{
	int i;
	for(i =	0; i < outframe->num_atoms; i++)
		{
		free(outframe->sites[i].observables);
		free(outframe->sites[i].coord);
		free(outframe->sites[i].matches);
		}	
	free(outframe->type);
	free(outframe->type_num);
	
	if(control->frame > 1) 
		{
		free(outframe->sites);
		}
}