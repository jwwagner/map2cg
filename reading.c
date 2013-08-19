//reading.c

#include "headers.h"

//main functions called from other files
void parse_command_line_arguments(int, char**, char*, char*, char*);
void read_topology_file(Controller*, char*);
void read_frame(Controller*, Frame*, FILE*, int*);
void read_frames_and_log(Controller*, Frame*, Frame*, FILE*, FILE*, FILE*, FILE*, int*);
void read_charge_frames(Controller*, Frame*, int*);
void read_charge_log(Controller*, int*);
void read_and_process_bootstrapping_trajectory(Controller*, Frame*, FILE*, int*);
void read_bootstrapping_file(Controller*, int*, int*);
//slave functions
void read_logfile(Controller*, FILE*, double*, int*);
void read_guess(Controller*, FILE*, int*);
void read_force_file(Controller*, char*, double*, double*,  int*);
void read_number_in_line(int, char*, int*);
void read_number_in_line_float(int, char*, double*);
void check_file_extension(char*, char*);
void report_traj_input_suffix_error(char*);
void report_usage_error(char*);
void eof_exit(Controller*, Frame*);
void eof_exit2(Controller*, FILE*, int*);

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
	int i, j;
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
    
    //read in prototype information
    fgets(line,100,fr);//7
    int* temp_mol = malloc(control->num_cg_types * sizeof(int));
    int* temp_type = malloc(control->max_to_map * sizeof(int));
    
    read_number_in_line(control->num_cg_types, line, temp_mol); 
	control->prototype = malloc(control->num_cg_types * sizeof(PROTO));
    for(i = 0; i < control->num_cg_types; i++)
    	{
    	control->prototype[i].num = temp_mol[i];
    	control->prototype[i].num_list = malloc(control->prototype[i].num * sizeof(int));
    	
    	fgets(line,100,fr);
    	read_number_in_line(control->prototype[i].num, line, temp_type);
    	
    	for(j = 0; j < control->prototype[i].num; j++)
    		{
    		control->prototype[i].num_list[j] = temp_type[j];
    		}
    	}
    free(temp_mol);
	free(temp_type);
	
    fgets(line,100,fr);//8+
    
    fgets(line,100,fr);//9+
    sscanf(line, "%d", &control->geometry_map_flag);
    
    fgets(line,100,fr);//10+
    sscanf(line, "%d", &control->observable_map_flag);
    
    fgets(line,100,fr);//11+
    sscanf(line, "%d", &control->num_observables);
    
    fgets(line,100,fr);//12+
    sscanf(line, "%d", &control->output_flag);
    //printf("O10:%s\n", line);
    
    fgets(line,100,fr);//13+
    sscanf(line, "%d", &control->sensitivity_flag);
    
    fgets(line,100,fr);//14+
    //printf("O11:%s\n", line);
    
    fgets(line,100,fr);//15+
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
    printf("sensitivity_flag %d\n", control->sensitivity_flag);
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
    else if(control->map_style_flag == 1)
    	{
    	printf("reading num_map\n");
    	fgets(line,100,fr);//12
    	fgets(line,100,fr);
    	sscanf(line, "%d", &control->num_map);
    	printf("line is %s\n", line);
    	
    	printf("num_map is %d\n", control->num_map);
    	control->map = malloc(control->num_map * sizeof(int));
    	for(i=0; i < control->num_map; i++)
    		{
    		fgets(line,100,fr);
    		sscanf(line, "%d", &control->map[i]);
    		//printf("reading %d as %d\n", i, control->map[i]);
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
    	
    	fgets(line,100,fr);//log type
    	sscanf(line, "%d", &control->log_type);
    	
    	fgets(line,100,fr);//guess type
    	sscanf(line, "%d", &control->guess_type);
    	
    	fgets(line,100,fr);//mapping type
    	sscanf(line, "%d", &control->sens_map_flag);
    	
    	fgets(line,100,fr);//debug flag
    	sscanf(line, "%d", &control->debug_flag);
    	
    	//set scaleF flag
    	double temp = 300.0 * 0.00198720414; //kcal/(mol K)
    	if( (control->debug_flag == 3) || (control->debug_flag == 5) ) control->scaleF = temp;
    	else control->scaleF = 1.0;
    	
    	//set scaleU flag
    	control->scaleU1 = 1.0;
    	control->scaleU2 = 1.0;
    	if( ((control->debug_flag >= 2) && (control->debug_flag <= 5)) || (control->debug_flag >= 7) )
    		{
    		control->scaleU1 *= temp;
    		control->scaleU2 *= temp;
    		}
    		
    	if( (control->debug_flag == 1) || (control->debug_flag == 4) || (control->debug_flag == 5) ) 
    		{
    		control->scaleU1 *= (double) control->num_cg_sites;
    		control->scaleU2 *= (double) control->num_cg_sites;
    		}
    	if( (control->debug_flag == 6) || (control->debug_flag == 7) )
    		{
    		control->scaleU1 *= (double) control->num_fg_sites; 
    		control->scaleU2 *= (double) control->num_cg_sites;
    		}
    		
    	if( (control->debug_flag == 8) )
    		{
    		control->scaleU1 *= ((double) control->num_cg_sites) * 3.0;
    		control->scaleU2 *= ((double) control->num_cg_sites) * 3.0;
    		}
    	if( (control->debug_flag == 9) )
    		{
    		control->scaleU1 *= ((double) control->num_fg_sites) * 3.0;
    		control->scaleU2 *= ((double) control->num_cg_sites) * 3.0;
    		}
    	    	
    	printf("control->debug_flag = %d\n", control->debug_flag);
    	printf("control->num_fg_sites = %d\n", control->num_fg_sites);
    	printf("control->num_cg_sites = %d\n", control->num_cg_sites);
    	printf("SCALE F = %lf\n", control->scaleF);
    	printf("SCALE U1 = %lf \n", control->scaleU1);
    	printf("SCALE U2 = %lf \n", control->scaleU2);
    	fgets(line,100,fr);//sign flag
    	sscanf(line, "%d", &control->sign_flag);
    	
    	if( abs(control->sign_flag) != 1 ) printf("ERROR: INVALID SIGN FLAG = %d\n", control->sign_flag);
    	else printf("control->sign_flag = %d\n", control->sign_flag);
    	}
    	
    else if(control->sensitivity_flag == 3)
    	{
    	fgets(line,100,fr);//blank line
    	
    	fgets(line, 100,fr);//1st dump file
    	sscanf(line,"%s", control->files.guess); //name of tabulated output
    	}
    
    else if( (control->sensitivity_flag == 2) || (control->sensitivity_flag == 4) )
    	{
		//routine specific variables
		char name[64];
		double* temp;
	
    	fgets(line,100,fr);//blank line
    	
    	fgets(line, 100,fr);//number of charges
    	sscanf(line,"%d", &control->num_charges);
    	control->num_files = 0;
    	for(i = 1; i <= control->num_charges; i++) control->num_files += i;
    	printf("num charges is %d and num_charges %d vs %d\n", control->num_charges, control->num_files, (int) round( control->num_files * ( (double) control->num_files * 0.5 + 0.5) ) );
    	
    	//read in charge values
    	temp = malloc( control->num_charges * sizeof(double) );
    	for(i = 0; i < control->num_charges; i++) 
    		{
    		temp[i] = 0;
    		}
    	control->charge = malloc( control->num_charges * sizeof(int) );
    	fgets(line, 100,fr);
    	read_number_in_line_float(control->num_charges, line, temp);
    	for(i = 0; i < control->num_charges; i++) 
    		{
    		printf("numbers read are %lf\n", temp[i] );
    		control->charge[i] = temp[i]; 
    		}
    	free(temp);

    	control->file_point = malloc( control->num_files * sizeof(FILE*) );
    	//read in the necessary number of filenames
    	printf("number of files expected is %d\n", control->num_files);
    	for(i = 0; i < control->num_files; i++)
    		{
    		fgets(line, 100,fr);
    		sscanf(line,"%s", name);
    		control->file_point[i] = fopen(name, "rt");
    		}
    		
    	//read # output files (as separator)
    	fgets(line, 100, fr);
    	sscanf(line,"%d", &i);
    	if(i != control->num_charges) printf("ERROR: number of CHARGES does not agree with number of OUTPUT files\n");
    	
    	//allocate space for output file pointers and read in (and open files)
    	control->outfile = malloc( control->num_charges * sizeof(FILE*) );
    	for(i = 0; i < control->num_charges; i++)
    		{
    		fgets(line, 100, fr);
    		sscanf(line, "%s", name);
    		control->outfile[i] = fopen( name, "w+");
    		}	
    	
    	}
    	
    else if(control->sensitivity_flag == 5)
    	{
    	printf("control->sensitivity_flag = %d\n", control->sensitivity_flag);
    	control->charge = malloc( 2 * sizeof(int) );
    	control->name = malloc( 32 * sizeof(char) );
    	
    	fgets(line,100,fr);//blank line
    	
    	fgets(line,100,fr);//start and end numbers for bootstrap files
    	sscanf(line,"%d %d", &control->charge[0], &control->charge[1]); //name of tabulated output
    	printf("bootstrap:%d %d\n", control->charge[0], control->charge[1]);
    	
    	fgets(line, 100,fr); //number of frames
    	sscanf(line,"%d", &control->num_frames);
    	
    	fgets(line,100,fr);//base for bootstrap filenames
    	sscanf(line,"%s", control->name);
    	
    	}

	fclose(fr);
    printf("finished reading top file\n\n");
}


//////////////////////////
///   read_frame	  ///
/////////////////////////

void read_frame(Controller* control, Frame* frame, FILE* df, int* flag)
{
	//printf("in read\n");
	int i, j;
	char line[100];
	char test[15];
	
	//printf("control->num_observables is %d\n", control->num_observables);
	//printf("frame->num_observables is %d\n", frame->num_observables);
	frame->num_observables = control->num_observables;
	
	//printf("begin reading frame\n");
	//check if EOF
	if( fgets (line, 100, df) == NULL ) //1
		{
		*flag = 0;
		printf("end of file reached in read_frame\n");
		eof_exit(control, frame);
		printf("finished eof_exit\n");
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
			eof_exit(control, frame);
			return;
		}
	
	//printf("passed EOF and format test\n");
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
    		//printf("intial allocation of atoms \n");
    		if(i != control->num_fg_sites) 
    			{
    				if(control->sens_map_flag == 1) 
    					{
    					if(i != control->num_cg_sites) printf("num cg sites = %d does not match num atoms %d\n", control->num_cg_sites, i);
    					}
    				else
    					{
    					printf("num fg sites = %d does not match num atoms %d\n", control->num_fg_sites, i);
    					}
    			}
    		frame->num_atoms = i;
			frame->atoms = malloc(frame->num_atoms * sizeof(ATOM));
			frame->num_observables = control->num_observables;
			for(j = 0; j < frame->num_atoms; j++) frame->atoms[j].mol = 0;
			//printf("set frame->num_observables to %d\n", frame->num_observables);
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
			for(j = 0; j < frame->num_atoms; j++) frame->atoms[j].mol = 0;
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
		//if( i < 10  && control->frame == 1) printf("read frame->atoms[%d].x = %lf\n", i, frame->atoms[i].x);
		
		if( (frame->atoms[i].mol > frame->num_mol) || (frame->num_mol > frame->num_atoms) ) frame->num_mol = frame->atoms[i].mol;
		//printf("frame->atoms[%d].mol %d and frame->num_mol %d\n", i, frame->atoms[i].mol, frame->num_mol);

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
	double double_val = 0.0;
	
	//read frames as usual
	//printf("read frame 1\n");
	read_frame(control, inframe1, df1, &test1);
	if(test1 == 0) 
		{
		*flag = 0;
		return;	
		}
	
	//printf("read frame 2\n");
	read_frame(control, inframe2, df2, &test2);
	if(test2 == 0)
		{
		*flag = 0;
		return;
		}
		
	printf("1:inframe1->num_atoms = %d and inframe2->num_atoms = %d\n", inframe1->num_atoms, inframe2->num_atoms);
	
	//printf("moving to logfile read\n");
	read_logfile(control, lf, &double_val, &test_log);
	control->log_value = double_val;
	
	if(test_log == 0) 
		{
		*flag = 0;
		return;
		}
	
	//printf("moving to guess read\n");
	read_guess(control, gf, &guess_read);
	if(guess_read == 0)
		{
		*flag = 0;
		return;
		}
	
	//do test on flags to determine composite (0 in any is a fail)
	*flag = test1 && test2 && test_log && guess_read;
	//printf("flag result is %d\n", *flag);
	
	printf("2:inframe1->num_atoms = %d and inframe2->num_atoms = %d\n", inframe1->num_atoms, inframe2->num_atoms);
}

  ////////////////////////////////
 //   read_charge_frames	  ///
////////////////////////////////
void read_charge_frames(Controller* control, Frame* inframes, int* flag)
{
	//declare variables
	int i;
	int temp = 1;
	
	//read each frame in turn
	for(i = 0; i < control->num_files; i++)
		{
		//printf("read frame %d\n", i);
		read_frame(control, &inframes[i], control->file_point[i], &temp);
		//printf("result is %d\n", temp);
		if(temp == 0)
			{
			*flag = 0;
			return;
			}
		}
}

  //////////////////////////
 //   read_charge_log	  ///
//////////////////////////////

void read_charge_log(Controller* control, int* flag)
{
	//declare variables
	int i;
	int temp = 1;
	double double_temp = 0.0;
	
	//read in value for each of the log files and check that the process completed correctly
	for(i = 0; i < control->num_files; i++)
		{
		read_logfile(control, control->file_point[i], &double_temp, &temp);
		control->log_values[i] = double_temp;
		
		if(temp == 0)
			{
			*flag = 0;
			return;
			}
		}
}

  ////////////////////////////////////////////////////
 //   read_and_process_bootstrapping_trajectory	  ///
/////////////////////////////////////////////////////

void read_and_process_bootstrapping_trajectory(Controller* control, Frame* outframes, FILE* fp, int* flag)
{
	//declare variables
	int i;
	int temp = 1;
	Frame inframe;
	inframe.num_atoms = 0;
	inframe.num_mol = 0;
	control->frame = 1;
	
	for(i=0; i < control->num_frames; i++)
		{
		//read each frame into the appropriate structure
		outframes[i].num_atoms = 0;
		//printf("read_frame #%d\n", i);
		read_frame(control, &inframe, fp, &temp);
	
		if(temp == 0)
			{
			*flag = 0;
			break;
			}
		
		//map that fine-grained frame to the appropriate coarse-grained frame
		printf("process_frame #%d\n", i);
		process_frame(control, &inframe, &(outframes[i]) );
		}
	
	//free temp inframe
	for(i = 0; i < inframe.num_atoms; i++)
		{
		free(inframe.atoms[i].observables);	
		}
	free(inframe.atoms);
}

  ///////////////////////////////////
 //   read_bootstrapping_file	  ///
/////////////////////////////////////

void read_bootstrapping_file(Controller* control, int* frame_order, int* count)
{
	//declare variables
	int i;
	char filename[32] = "";
	char line[100];
	char suffix[] = ".dat"; 
	char num[5];
	FILE* fp;
	
	//determine file to read
	snprintf(filename, sizeof filename, "%s%d%s", control->name, *count, suffix);
	
	//printf("file name to read is %s:\n", filename);
	//open file
	fp = fopen(filename, "rt");
	//printf("file open for %s with result %p\n", filename, fp);
	
	//read file into array
	for(i = 0; i < control->num_frames; i++)
		{
		//printf("read value to site %d\n", i);
		fgets(line,100,fp);
		//printf("READ:%s\n", line);
		sscanf(line, "%d", &(frame_order)[i] );
		//printf("read %d as %d\n", i, frame_order[i]);
		}
	
	//close file
	fclose(fp);
	//printf("file closed\n");
}


  ////////////////////////
 //   read_logfile	  ///
////////////////////////

void read_logfile(Controller* control, FILE* lf, double* out_val, int* flag)
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
		eof_exit2(control, lf, flag);
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
			//fgets(line, 100, lf);
			
			//check if we are at EOF
			if( fgets (line, 100, lf) == NULL )
			{
				*flag = 0;
				printf("end of file reached in logfile!\n");
				eof_exit2(control, lf, flag);
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
	
	if(control->log_type == 0) *out_val = val1;
	else if(control->log_type == 1) *out_val = val2;
	
	control->timestep = junk;
	control->volume = junk2;
}

/////////////////////////
///   read_guess	  ///
////////////////////////

void read_guess(Controller* control, FILE* gf, int* flag)
{
	char line[100];
	int i;
	
	if(control->guess_type == 0)
		{
		//check if we are at EOF
		if( fgets (line, 100, gf) == NULL )
			{
			*flag = 0;
			printf("end of file reached in GUESS file!\n");
			eof_exit2(control, gf, flag);
			return;
			}

		//read guess value
		sscanf(line, "%lf", &control->guess);
		}

	else if(control->guess_type == 1)
		{
		//check if we are at EOF
		if( fgets (line, 100, gf) == NULL )
			{
			*flag = 0;
			printf("end of file reached in GUESS file!\n");
			eof_exit2(control, gf, flag);
			return;
			}

		//check if we need to skip header info
		if(control->frame == 1)
			{
			char test[8];
			double junk, junk1, junk2, val1, val2;
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
				fgets(line, 100, gf);
			
				//check if we are at EOF
				if( fgets (line, 100, gf) == NULL )
					{
					*flag = 0;
					printf("end of file reached in GUESS!\n");
					eof_exit2(control, gf, flag);
					return;
					}
				}

			//see if this is a WARNING line or actual content
			//printf("\nwarningtest  logfile line reads %s\n", line);
			memcpy(test, &line[0], 7);
			test[7] = '\0';
			if( strcmp(test, "WARNING") == 0) fgets(line, 100, gf);
	
			//read actual content for energy value
			//printf("acutal logfile line reads %s\n", line);
			sscanf(line, "%lf %lf %lf %lf %lf", &junk, &junk1, &val1, &val2, &junk2);
			//printf("values read are %lf %lf %lf %lf %lf\n", junk, junk1, val1, val2, junk2); 
			if(control->log_type == 0) control->guess = val1;
			else if(control->log_type == 1) control->guess = val2;
			}
		}
}

//////////////////////////////
///   read_force_file	  ///
////////////////////////////


void read_force_file(Controller* control, char* file, double* distance, double* force,  int* number_of_lines)
{
	//declare variables
	int i;
	int num_lines = 0;
	int flag = 1;
	char line[100];
	FILE* fp = fopen(file, "rt");
	
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
	//num_lines--;
		
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
	
	*number_of_lines = num_lines;
	printf("*number of lines %d = %d num lines\n", *number_of_lines, num_lines);

}

//////////////////////////////////
///   read_number_in_line	  ///
////////////////////////////////


void read_number_in_line(int num_vals, char* line, int* vals)
{
	switch(num_vals)
		{
		case 1:
			sscanf(line, "%d", &vals[0]);
			break;
		
		case 2:
			sscanf(line, "%d %d", &vals[0], &vals[1]);
			break;
		
		case 3:
			sscanf(line, "%d %d %d", &vals[0], &vals[1], &vals[2]);
			break;
		
		case 4:
			sscanf(line, "%d %d %d %d", &vals[0], &vals[1], &vals[2], &vals[3]);
			break;
		
		case 5:
			sscanf(line, "%d %d %d %d %d", &vals[0], &vals[1], &vals[2], &vals[3], &vals[4]);
			break;
		
		case 6:
			sscanf(line, "%d %d %d %d %d %d", &vals[0], &vals[1], &vals[2], &vals[3], &vals[4], &vals[5]);
			break;
		
		case 7:
			sscanf(line, "%d %d %d %d %d %d %d", &vals[0], &vals[1], &vals[2], &vals[3], &vals[4], &vals[5], &vals[6]);
			break;
		
		case 8:
			sscanf(line, "%d %d %d %d %d %d %d %d", &vals[0], &vals[1], &vals[2], &vals[3], &vals[4], &vals[5], &vals[6], &vals[7]);
			break;
		
		case 9:
			sscanf(line, "%d %d %d %d %d %d %d %d %d", &vals[0], &vals[1], &vals[2], &vals[3], &vals[4], &vals[5], &vals[6], &vals[7], &vals[8]);
			break;
		
		case 0:	
		default:
			printf("zero or unsupported number of elements to read in\n");
			break;
		}
}

void read_number_in_line_float(int num_vals, char* line, double* vals)
{
	switch(num_vals)
		{
		case 1:
			sscanf(line, "%lf", &vals[0]);
			break;
		
		case 2:
			sscanf(line, "%lf %lf", &vals[0], &vals[1]);
			break;
		
		case 3:
			sscanf(line, "%lf %lf %lf", &vals[0], &vals[1], &vals[2]);
			break;
		
		case 4:
			sscanf(line, "%lf %lf %lf %lf", &vals[0], &vals[1], &vals[2], &vals[3]);
			break;
		
		case 5:
			sscanf(line, "%lf %lf %lf %lf %lf", &vals[0], &vals[1], &vals[2], &vals[3], &vals[4]);
			break;
		
		case 6:
			sscanf(line, "%lf %lf %lf %lf %lf %lf", &vals[0], &vals[1], &vals[2], &vals[3], &vals[4], &vals[5]);
			break;
		
		case 7:
			sscanf(line, "%lf %lf %lf %lf %lf %lf %lf", &vals[0], &vals[1], &vals[2], &vals[3], &vals[4], &vals[5], &vals[6]);
			break;
		
		case 8:
			sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf", &vals[0], &vals[1], &vals[2], &vals[3], &vals[4], &vals[5], &vals[6], &vals[7]);
			break;
		
		case 9:
			sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", &vals[0], &vals[1], &vals[2], &vals[3], &vals[4], &vals[5], &vals[6], &vals[7], &vals[8]);
			break;
		
		case 0:	
		default:
			printf("zero or unsupported number of elements to read in\n");
			break;
		}
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

void eof_exit(Controller* control, Frame* frame)
{
	//see if we need to allocate space so normal free is okay
	if(control->frame <= 1)
		{
		frame->atoms = malloc(1 * sizeof(ATOM));
		frame->atoms[0].observables = malloc(1 * sizeof(double));
		}
}

void eof_exit2(Controller* control, FILE* lf, int* flag)
{
}
