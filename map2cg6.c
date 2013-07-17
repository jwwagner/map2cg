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
// int *** #number of atoms in site  for #cg types listed above
// int *** #list types for each atom in cg molecule (# should match above line's value)
// ...     #continue until #CG types is fulfilled
//		#this line intentionally left blank
// int	#geometry mapping flag (0 = CoM geometry mapping, 1 = center of geometry mapping)
// int  #observable mapping flag (0 = additive mappign)
// int	#number observables to map
// int  #output flag (0 = all, 1 = minimal, 2 = 1 value is #4 observable, 3 = 2nd sets of 3 values)
// int  #sensitivity flag (0 = only map 1 file, 1 = map 2 files along with log data and guess, 2 unused, 3 = convert MSCGFM (1_1.dat) to (tab_ff.dat))
//	 	#this line intentionally left blank
// int  #mapping style flag (0 = 1:1 molecule entirely, 1 = implicit solvent)
//IF mapping style flag == 1
// int	#number of types to map
// int #type ID's to be mapped
// ...
//ENDIF
//		#this line intentionally left blank
//IF sensitivity flag == 1
// %s filename for 1st dump (original data, regular parameters)
// %s filename for 2nd dump (rerun data, force derivative wrt parameter, e.g., epsilon_11)
// %s filename for log.lammps (rerun data, potential derivative wrt parameter, e.g., epsilon_11)
// %s filename for guess of dU/d\lambda (needs as many rows as frames in other inputs)
// %s filename for composite output (lammps dump file)
// int	#log type (0 = vdW, 1 = Columb)
// int  #guess type (0 = 1 column, 1 = log file similar log.lammps file)
//ENDIF
//IF sensitivity flag == 3
// %s name of ID for tab potential (e.g. SENS_MEOH)
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

//structure defintions
#include "headers.h"

//external files
#include "output_lmp.c"
#include "process_frame.c"
#include "reading.c"

//function prototype definitions
void do_simple_map(Controller*, Frame*, Frame*, char*, char*);
void do_sensitivity_map(Controller*, Frame*, Frame*, char*);
void do_force_conversion(Controller*, char*, char*);
void free_allocation(Controller*, Frame*, Frame*);
void free_sensitivity_allocation(Controller*, Frame*, Frame*, Frame*); 

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
	outframe.num_atoms = 0;
	controls.frame = 0;
	
	char* datfile = malloc(64 * sizeof(char));
	char* topfile = malloc(64 * sizeof(char));
	char* outfile = malloc(64 * sizeof(char));
	
	//read topology/control file
	parse_command_line_arguments(argc, argv, datfile, topfile, outfile);
 
	printf("assigned file names are %s and %s and %s\n", datfile, topfile, outfile);
	
	//read topology/control file
	read_topology_file(&controls, topfile);
	
	//determine appropriate controller function
	if(controls.sensitivity_flag == 0) do_simple_map(&controls, &frame, &outframe, datfile, outfile);
	else if(controls.sensitivity_flag == 1) do_sensitivity_map(&controls, &frame, &outframe, outfile);	
	else if(controls.sensitivity_flag == 3) do_force_conversion(&controls, datfile, outfile);
	
	//free allocated variables
	printf("free basic files\n");
	free(datfile);
	free(topfile);
	free(outfile);
	
	printf("free map information\n");
	if( (controls.map_style_flag == 0) || (controls.map_style_flag == 1) ) free(controls.map);
		
	if(controls.sensitivity_flag == 0)
		{
		printf("to free_allocation\n");
		free_allocation(&controls, &frame, &outframe);
		}
		
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
	
	printf("main says num_observables is %d\n", controls->num_observables);
	
	while(cont_flag == 1)
		{
		//process/map frame
		process_frame(controls, inframe, outframe);
		//printf("finished processing frame\n");
		
		//output mapped frame and observables
		output_frame(controls, outframe, outfile);
		//printf("finished output for frame %d\n", frame_count);
		
		//read next frame or set flag if done
		frame_count++;
		controls->frame++;
		read_frame(controls, inframe, df, &cont_flag);
		//printf("cont_flag is %d\n", cont_flag);
		}
		
	//close dump file
	fclose(df);
	
	//create top.in file for FM
	output_topology(controls, outframe);
	
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
	if(df1 == NULL) printf("dump file #1 specified does not exist\n");
	if(df2 == NULL) printf("dump file #2 specified does not exist\n");
	if(lf == NULL) printf("log file specified does not exist\n");
	if(gf == NULL) printf("guess file specified does not exist\n");

	//exit if there is an error
	if( (df1 == NULL) || (df2 == NULL) || (lf == NULL) || (gf == NULL) )  exit(EXIT_SUCCESS);
		
	//read first frames
	frame_count++;
	controls->frame++;
	
	printf("reading first frame data for sensitivity map\n");
	read_frames_and_log(controls, inframe1, &inframe2, df1, df2, lf, gf, &cont_flag);
	
	printf("main says num_observables is %d\n", controls->num_observables);
	printf("1st reading reports log value is %lf and guess is %lf\n", controls->log_value, controls->guess);
	
	while(cont_flag == 1)
		{
		//process/map frame
		process_frames_and_log(controls, inframe1, &inframe2, outframe);
		//printf("finished processing frame\n");
		
		//output mapped frame and observables
		output_frame(controls, outframe, outfile);
		//printf("finished output for frame %d\n", frame_count);
		
		//read next frame or set flag if done
		frame_count++;
		controls->frame++;
		read_frames_and_log(controls, inframe1, &inframe2, df1, df2, lf, gf, &cont_flag);
		//printf("cont_flag is %d\n", cont_flag);
		}

	//close files for frame, log, and guess reading
	printf("close df1\n");
	fclose(df1);
	printf("close df2\n");
	fclose(df2);
	printf("close lf\n");
	fclose(lf);
	printf("closes gf\n");
	fclose(gf);
	
	printf("to output_topology\n");
	//create top.in file for FM
	output_topology(controls, outframe);
	
	printf("to free_sensitivity_allocation");
	//free allocated varaibles
	free_sensitivity_allocation(controls, inframe1, &inframe2, outframe); 
	printf("finished with sensitivity function\n");
}

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
			//printf("num_lines is %d\n", num_lines);
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
		//printf("%d: %s\n", i, line);
		sscanf(line, "%lf %lf", &distance[i], &force[i]);
		printf("read at site %d values of distance %lf and force %lf\n", i, distance[i], force[i]);
		}
	
	//close input file
	fclose(fp);
	
	num_entries = num_lines;
	printf("number of lines %d = %d num lines\n", num_entries, num_lines);
	
	//read_force_file(control, infile, distance, force, &num_entries);
	
	//printf("returned values include num_entries %d\n", num_entries);	
	//printf("pointer values are distance %d force %d and potential %d\n", &distance, &force, &potential);
	//printf("first values are distance[0] = %lf\n", distance[0]);
	//printf("first values are force[0] = %lf\n", force[0]);
	
	//integrate to get potential
	potential = malloc(num_entries * sizeof(double));
	
	printf("integrate\n");
	potential[num_entries - 1] = 0.0;
	for(i = num_entries - 2; i >= 0; i--)
		{
		//printf("set value at site i = %d\n", i);
		//printf("use value of distance 2 %lf\n", distance[i]);
		//printf("use value of distance 1 %lf\n", distance[i+1]);
		//printf("use value of force 1 %lf\n", force[i+1]);
		//printf("use value of force 2 %lf\n", force[i]);
		//printf("set value of potential %lf\n", potential[i]);
		
		potential[i] = potential[i+1] + 0.5 * (distance[i + 1] - distance[i]) * (force[i] + force[i + 1]);
		}
	
	//convert ALL units
	printf("convert units\n");
	for(i = 0; i > num_entries; i++)
		{
		//distance (nm -> A)
		distance[i] *= 10.0;
		
		//potential (kj/mol -> kcal/mol)
		potential[i] /= 4.184;
		
		//force (kj / mol nm -> kcal/mol A)
		force[i] /= 41.84;
		}
	
	printf("write_output\n");
	//write output
	output_force_file(control, num_entries, distance, potential, force, outfile);
	
	printf("free variables\n");
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

//////////////////////////////
///   free_allocation	  ///
////////////////////////////

void free_allocation(Controller* control, Frame* inframe, Frame* outframe)
{
	int i;
	
	//free controller data alloacted
	for(i = 0; i < control->num_cg_types; i++)
		{
		free(control->prototype[i].num_list);
		}
	free(control->prototype);
		
	
	//free all data for atoms
	printf("free 1 inframe\n");
	for(i = 0; i < inframe->num_atoms; i++)
		{
		free(inframe->atoms[i].observables);	
		}
	free(inframe->atoms);

	printf("free 1 outframe\n");
	for(i =	0; i < outframe->num_atoms; i++)
		{
		free(outframe->sites[i].observables);
		free(outframe->sites[i].coord);
		}
	free(outframe->sites);
	free(outframe->type);
	free(outframe->type_num);
	printf("finished with free_allocation\n");
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
	printf("free inframe information\n");
	for(i = 0; i < inframe1->num_atoms; i++)
		{
		free(inframe1->atoms[i].observables);	
		free(inframe2->atoms[i].observables);
		}
	free(inframe1->atoms);
	free(inframe2->atoms);

	
	//free all data for atoms
	printf("free outframe information\n");
	for(i = 0; i < outframe->num_atoms; i++)
		{
		free(outframe->sites[i].observables);
		}
	
	printf("free sites\n");
	free(outframe->sites);
	
	printf("free type\n");
	free(outframe->type);
	
	printf("free type_num");
	free(outframe->type_num);
}
