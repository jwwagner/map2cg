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
// int  #output flag (0 = all, 1 = minimal, 2 = 1 value is #4 observable, 3 = 2nd sets of 3 values)
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

//structure defintions
#include "headers.h"

//external files
#include "output_lmp.c"
#include "process_frame.c"
#include "reading.c"

//function prototype definitions
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
 
	printf("assigned file names are %s and %s and %s\n", datfile, topfile, outfile);
	
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
		//printf("finished processing frame\n");
		
		//output mapped frame and observables
		output_frame(&controls, &outframe, outfile);
		//printf("finished output for frame %d\n", frame_count);
		
		//read next frame or set flag if done
		frame_count++;
		controls.frame++;
		read_frame(&controls, &frame, df, &cont_flag);
		//printf("cont_flag is %d\n", cont_flag);
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