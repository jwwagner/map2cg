//reading.c

#include "headers.h"

void parse_command_line_arguments(int, char**, char*, char*, char*);
void read_topology_file(Controller*, char*);
void read_frame(Controller*, Frame*, FILE*, int*);
void read_frames_and_log(Controller*, Frame*, Frame*, FILE*, FILE*, FILE*, FILE*, int*);
void read_logfile(Controller*, FILE*, int*);
void read_guess(Controller*, FILE*, int*);
void check_file_extension(char*, char*);
void report_traj_input_suffix_error(char*);
void report_usage_error(char*);

//////////////////////////////////////////
///   parse_command_line_arguments	  ///
////////////////////////////////////////

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

//////////////////////////////////
///   read_topology_file	  ///
////////////////////////////////

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
    sscanf(line, "%d", &control->output_flag);
    //printf("O10:%s\n", line);
    
    fgets(line,100,fr);//11
    sscanf(line, "%d", &control->sensitivity_flag);
    
    fgets(line,100,fr);//12
    //printf("O11:%s\n", line);
    
    fgets(line,100,fr);//13
    //printf("O12:%s\n", line);
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
    printf("output_flag %d \n", control->output_flag);
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
    	
    if(control->sensitivity_flag == 1)
    	{
    	fgets(line,100,fr);//blank line    	
    	
    	fgets(line,100,fr);//1st dump file
    	sscanf(line, "%s", control->files.dump1);
    	
    	fgets(line,100,fr);//2nd dump file
    	sscanf(line, "%s", control->files.dump2);
    	
    	fgets(line,100,fr);//log file
    	sscanf(line, "%s", control->files.log);
    	
    	fgets(line,100,fr);//guess file
    	sscanf(line, "%s", control->files.guess);
    	
    	fgets(line,100,fr);//guess file
    	sscanf(line, "%d", &control->log_type);
    	}
    	
    printf("finished reading top file\n\n");
}


//////////////////////////
///   read_frame	  ///
/////////////////////////

void read_frame(Controller* control, Frame* frame, FILE* df, int* flag)
{
	int i, j;
	char line[100];
	char test[15];
	frame->num_observables = control->num_observables;
	
	//printf("begin reading frame\n");
	//check if EOF
	if( fgets (line, 100, df) == NULL ) //1
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
	fgets(line,100,df); //2
    sscanf(line, "%d", &frame->timestep);
    
    //read in number of atoms
    fgets(line,100,df); //3 
	fgets(line,100,df); //4
    sscanf(line, "%d", &i);
    
    printf("timestep is %d\n", frame->timestep);
    //printf("number of atoms is %d\n", i);
    
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
			frame->num_mol = 0;
			}
			
		//allocate observable information
		for(i=0; i < frame->num_atoms; i++)
			{
			frame->atoms[i].observables = malloc(frame->num_observables * sizeof(double));
			}
		printf("observable space allocated\n");
		}
	
	//read in box information
	fgets(line,100,df); //5
	fgets(line,100,df); //6
    sscanf(line, "%lf %lf", &frame->xmin, &frame->xmax);
	fgets(line,100,df); //7
    sscanf(line, "%lf %lf", &frame->ymin, &frame->ymax);
	fgets(line,100,df); //8
    sscanf(line, "%lf %lf", &frame->zmin, &frame->zmax);
    
    //printf("box is %lf %lf by %lf %lf by %lf %lf\n", frame->xmin, frame->xmax, frame->ymin, frame->ymax, frame->zmin, frame->zmax);
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
		if( (frame->atoms[i].mol > frame->num_mol) || (frame->num_mol > frame->num_atoms) ) frame->num_mol = frame->atoms[i].mol;
		//printf("frame->atoms[%d].mol %d and frame->num_mol %d\n", i, frame->atoms[i].mol, frame->num_mol);
		//if(frame->atoms[i].mol > frame->num_mol)	frame->num_mol = frame->atoms[i].mol;
	}
	//finished reading frame
}

//////////////////////////////////
///   read_frames_and_log	  ///
/////////////////////////////////

void read_frames_and_log(Controller* control, Frame* inframe1, Frame* inframe2, FILE* df1, FILE* df2, FILE* lf, FILE* gf, int* flag)
{
	//flag variables
	int test1 = 1;
	int test2 = 1;
	int test_log = 1;
	int guess_read = 1;
	
	//read frames as usual
	//printf("read frame 1\n");
	read_frame(control, inframe1, df1, &test1);
	//printf("read frame 2\n");
	read_frame(control, inframe2, df2, &test2);
	//printf("moving to logfile read\n");
	read_logfile(control, lf, &test_log);
	//printf("moving to guess read\n");
	read_guess(control, gf, &guess_read);
	//do test on flags to determine composite (0 in any is a fail)
	*flag = test1 && test2 && test_log && guess_read;
	//printf("flag result is %d\n", *flag);
}

  ////////////////////////
 //   read_logfile	  ///
////////////////////////

void read_logfile(Controller* control, FILE* lf, int* flag)
{
	//declare varaibles
	int i;
	double junk, junk1, junk2;
	double val1, val2;
	char line[100];
	char test[8];
	
	//check if we are at EOF
	if( fgets(line, 100, lf) == NULL )
		{
		*flag = 0;
		printf("end of file reached in logfile!\n");
		return;
		}

	//check if we need to skip header info
	if(control->frame == 1)
	{
		i = 1;
		char test2[5];
		
		while(i == 1)
		{
			//check line
			//printf("intial log file line skipped is %s\n", line);
			memcpy(test2, &line[0], 4);
			test2[4] = '\0';
			if( strcmp(test2, "Step") == 0) i = 0;			
			//read next line
			fgets(line, 100, lf);
			
			//check if we are at EOF
			if( fgets (line, 100, lf) == NULL )
			{
				*flag = 0;
				printf("end of file reached in logfile!\n");
				return;
			}

		}
	}
	
	//see if this is a WARNING line or actual content
	//printf("\nwarningtest  logfile line reads %s\n", line);
	memcpy(test, &line[0], 7);
	test[7] = '\0';
	if( strcmp(test, "WARNING") == 0) fgets(line, 100, lf);
	
	//read actual content for energy value
	//printf("acutal logfile line reads %s\n", line);
	sscanf(line, "%lf %lf %lf %lf %lf", &junk, &junk1, &val1, &val2, &junk2);
	//printf("values read are %lf %lf %lf %lf %lf\n", junk, junk1, val1, val2, junk2); 
	if(control->log_type == 0) control->log_value = val1;
	else if(control->log_type == 1) control->log_value = val2;
}

/////////////////////////
///   read_guess	  ///
////////////////////////

void read_guess(Controller* control, FILE* gf, int* flag)
{
	char line[100];
	
	//check if we are at EOF
	if( fgets (line, 100, gf) == NULL )
		{
		*flag = 0;
		printf("end of file reached in GUESS file!\n");
		return;
		}

	//read guess value
	sscanf(line, "%lf", &control->guess);
}

/////////////////////////////////////
///   err_checking_functions	  ///
////////////////////////////////////

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
