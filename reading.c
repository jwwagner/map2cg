//reading.c
#include "headers.h"

//main functions called from other files
void parse_command_line_arguments(int, char**, char*, char*, char*);
void read_topology_file(Controller*, char*);
void read_frame(Controller*, Frame*, FILE*, int*);
void read_frame_minimal(Controller*, Frame*, FILE*, int*);
void read_frames_and_log(Controller*, Frame*, Frame*, FILE*, FILE*, FILE*, FILE*, int*);
void read_charge_frames(Controller*, Frame*, int*);

//slave functions
void read_logfile(Controller*, FILE*, double*, int*);
void read_guess(Controller*, FILE*, int*);
void read_number_in_line(int, char*, int*);
void read_number_in_line_float(int, char*, double*);
void check_file_extension(char*, char*);
void report_traj_input_suffix_error(char*);
void report_usage_error(char*);
void eof_exit(Controller*, Frame*);

//function_pointed reading function
void read0full(Frame*, FILE*);
void read3full(Frame*, FILE*);
void read0min(Frame*, FILE*);
void read3min(Frame*, FILE*);

//////////////////////////////////////////
///   parse_command_line_arguments	  ///
////////////////////////////////////////

void parse_command_line_arguments(int num_arg, char** arg, char* filename1, char* filename2, char* outfile)
{
    if (num_arg != 1 && num_arg != 3 && num_arg != 5 && num_arg != 7) {
    	report_usage_error(arg[0]);
    }
    else if (num_arg == 3) {
        if (strcmp(arg[1], "-f") == 0) {
        	sscanf(arg[2], "%s", filename1);
        	check_file_extension(arg[2], "dat");
        	strcpy(filename2, "map.top");
        	strcpy(outfile, "out.dat");
        } else if  (strcmp(arg[1], "-f1") == 0) {
        	sscanf(arg[2], "%s", filename2);
        	check_file_extension(arg[2], "top");
        	strcpy(filename1, "dump.dat");
        	strcpy(outfile, "out.dat");
        } else if (strcmp(arg[1], "-o") == 0) {
        	sscanf(arg[2], "%s", outfile);
        	strcpy(filename1, "dump.dat");
        	strcpy(filename2, "map.top");
        } else  report_usage_error(arg[0]);
    } else if (num_arg == 5) {
        if ( ((strcmp(arg[1], "-f") != 0) || (strcmp(arg[3], "-f1") != 0)) && ((strcmp(arg[1], "-f1") != 0) || (strcmp(arg[3], "-o") != 0)) && ((strcmp(arg[1], "-f") != 0) || (strcmp(arg[3], "-o") != 0)) ) {
        	report_usage_error(arg[0]);
        } else if ( strcmp(arg[1], "-f") == 0 && strcmp(arg[3], "-f1") == 0 ) {
        	sscanf(arg[2], "%s", filename1);
        	sscanf(arg[4], "%s", filename2);
        	check_file_extension(arg[2], "dat");
        	check_file_extension(arg[4], "top");
        	strcpy(outfile, "out.dat");
        } else if ( strcmp(arg[1], "-f") == 0 && strcmp(arg[3], "-o") == 0 ) {
        	sscanf(arg[2], "%s", filename1);
        	sscanf(arg[4], "%s", outfile);
        	check_file_extension(arg[2], "dat");
        	strcpy(filename2, "map.top");
        } else if (strcmp(arg[1], "-f1") == 0 && strcmp(arg[3], "-o") == 0 ) {
     		sscanf(arg[2], "%s", filename1);
        	sscanf(arg[4], "%s", outfile);
        	check_file_extension(arg[2], "top");
     		strcpy(filename2, "map.top");
     	}
     } else if (num_arg == 7) {
        if (strcmp(arg[1], "-f") != 0 || strcmp(arg[3], "-f1") != 0 || strcmp(arg[5], "-o") != 0) {
        	report_usage_error(arg[0]);
		} else {
			sscanf(arg[2], "%s", filename1);
        	sscanf(arg[4], "%s", filename2);
        	sscanf(arg[6], "%s", outfile);
        	check_file_extension(arg[2], "dat");
        	check_file_extension(arg[4], "top");
 		}
	} else if (num_arg == 1) {
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
	int i=0;
	int j=0;
	FILE* fr = fopen(topfile, "rt");
    char line[100];
    
    //initialize
    control->num_cg_types = 0;
    control->num_fg_types = 0;
    control->map_style_flag = 0;
    control->sensitivity_flag = 0;
    control->num_charges = 0;
	control->num_files = 0;
	control->num_observables = 0;
    
    //read_input
    fgets(line,100,fr);//1
    sscanf(line, "%d", &control->num_cg_sites);
    
    fgets(line,100,fr);//2
    sscanf(line, "%d", &control->num_cg_types);
    
    fgets(line,100,fr);//3
     
    fgets(line,100,fr);//4
    sscanf(line, "%d", &control->num_observables);
    
    fgets(line,100,fr);//5
    sscanf(line, "%d", &control->output_flag);
    
    fgets(line,100,fr);//6
    sscanf(line, "%d", &control->sensitivity_flag);
    
    fgets(line,100,fr);//14+
    
    //print information read in file
    printf("num_cg_sites %d\n", control->num_cg_sites);
    printf("num_cg_types %d\n", control->num_cg_types);
    printf("\n");
    printf("num_observables %d\n", control->num_observables);
    printf("output_flag %d \n", control->output_flag);
    printf("sensitivity_flag %d\n", control->sensitivity_flag);
    
    //set function pointer for reading function
    control->header_function = &default_func;
    switch(control->num_observables) {
		case 0:
    		control->read_function = &read0full;
    		control->output_function = &out0full;
    		control->header_function = &header0full;
			break;
		case 3:
    		control->read_function = &read3full;
    		control->output_function = &out3full;
    		control->header_function = &header3full;
			break;
		default:
			printf("number of observables requested %d is not supported \n", control->num_observables);
			break;
	}
    	
    if(control->sensitivity_flag == 1) {
    	control->debug_flag = 0;
    	
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
    		
    	//set scaleU value
    	double temp = 300.0 * 0.00198720414; //kcal/(mol K)
    	control->scaleU = temp * (double) control->num_cg_sites;
		
    	printf("control->num_cg_sites = %d\n", control->num_cg_sites);
    	printf("SCALE U = %lf \n", control->scaleU);
    
	} else if(control->sensitivity_flag == 6) { //no parameters need for control->sensitivity_flag == 4  
		//routine specific variables
		char name[64];
		double* temp;
	    
	    printf("\nreading for sensitivity_flag = 6\n");
    	fgets(line,100,fr);//blank line
		sscanf(line, "%s", name);

    	//read in charge values
    	control->num_charges = 1;
    	control->charge = malloc( control->num_charges * sizeof(double) );
    	fgets(line,100,fr);//read charge
    	sscanf(line,"%lf", &control->charge[0]);
    	printf("line for charge is %s = %lf\n", line, control->charge[0]);
    	
    	//read in filenames
    	control->num_files = 1;
    	control->file_point = malloc( control->num_files * sizeof(FILE*) );
    	for(i = 0; i < control->num_files; i++){
    		fgets(line, 100,fr);
    		sscanf(line,"%s", name);
    		control->file_point[i] = fopen(name, "rt");
    		printf("infile: %s = %s\n", name, line);
    	}
    		    	
    	//allocate space for output file pointers and read in (and open files)
    	control->num_outfile = 1;
    	control->outfile = malloc( control->num_outfile * sizeof(FILE*) );
    	for(i = 0; i < control->num_outfile; i++) {
    		fgets(line, 100, fr);
    		sscanf(line, "%s", name);
    		control->outfile[i] = fopen( name, "w+");
    		printf("outfile: %s = %s\n", name, line);
    	}
    } else if(control->sensitivity_flag == 7) {
		//routine specific variables
		char name[64];
		double* temp;
	    
    	fgets(line,100,fr);//blank line
    	//read in charge values
    	control->num_charges = 2;
    	control->charge = malloc( control->num_charges * sizeof(double) );
    	fgets(line,100,fr);//read charge
    	sscanf(line,"%lf %lf", &control->charge[0], &control->charge[1]);
      	printf("line for charge is %s = %lf and %lf\n", line, control->charge[0], control->charge[1]);
 	
    	//read in filenames
    	control->num_files = 3;
    	control->file_point = malloc( control->num_files * sizeof(FILE*) );
    	for(i = 0; i < control->num_files; i++) {
    		fgets(line, 100,fr);
    		sscanf(line,"%s", name);
    		control->file_point[i] = fopen(name, "rt");
    		printf("infile: %s\n", name);
    	}
    		    	
    	//allocate space for output file pointers and read in (and open files)
    	control->num_outfile = 2;
    	control->outfile = malloc( control->num_outfile * sizeof(FILE*) );
    	for(i = 0; i < control->num_outfile; i++) {
    		fgets(line, 100, fr);
    		sscanf(line, "%s", name);
    		control->outfile[i] = fopen( name, "w+");
    		printf("outfile: %s\n", name);
    	}	 	
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
	int i=0;
	int j=0;
	char line[100];
	char test[15];
	
	frame->num_observables = control->num_observables;
	
	//check if EOF
	if( fgets (line, 100, df) == NULL ) { //1
		*flag = 0;
		printf("end of file reached in read_frame\n");
		eof_exit(control, frame);
		printf("finished eof_exit\n");
		return;
	}
	
	//check if content matches expected 1st line of frame
	memcpy( test, &line[0], 14);
	test[14] = '\0';
	if( strcmp(test, "ITEM: TIMESTEP") != 0) {
		printf("error in frame format\n");
		printf("line found was %s\n", test);
		printf("expected line is ITEM: TIMESTEP\n");
		*flag = 0;
		eof_exit(control, frame);
		return;
	}
	
	//so, we assume content is correct now for next frame and begin reading read in timestep
	fgets(line,100,df); //2
    sscanf(line, "%d", &frame->timestep);
    
    //read in number of atoms
    fgets(line,100,df); //3 
	fgets(line,100,df); //4
    sscanf(line, "%d", &i);
    
    printf("timestep is %d\n", frame->timestep);
    
    //check if num_atoms changed
    if(i != frame->num_atoms) {
    	if(frame->num_atoms == 0) {	//do initial allocation
    		if(i != control->num_fg_sites) {
    			if(control->sens_map_flag == 1) {
    				if(i != control->num_cg_sites) {
    						printf("num cg sites = %d does not match num atoms %d\n", control->num_cg_sites, i);
    				}
    			} else {
    				printf("num fg sites = %d does not match num atoms %d\n", control->num_fg_sites, i);
    			}
    		}
    		frame->num_atoms = i;
			frame->atoms = malloc(frame->num_atoms * sizeof(ATOM));
			frame->num_observables = control->num_observables;
			for(j = 0; j < frame->num_atoms; j++) { 
				frame->atoms[j].mol = 0; 
			}
		} else { //number of atoms has changed
			printf("number of atoms has  changed");
			//free all existing data for atoms
			for(j = 0; j < frame->num_atoms; j++) {
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
		for(i=0; i < frame->num_atoms; i++) {
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
    
    //read in atoms and observables for frame
	fgets(line,100,df);
	(*control->read_function)(frame, df);		
}

//////////////////////////////////
///   read_frame_minimal	  ///
/////////////////////////////////

void read_frame_minimal(Controller* control, Frame* frame, FILE* df, int* flag)
{
	int i=0;
	int j=0;
	char line[100];
	char test[15];
	
	frame->num_observables = control->num_observables;
	
	//check if EOF
	if( fgets (line, 100, df) == NULL ) { //1
		*flag = 0;
		printf("end of file reached in read_frame\n");
		eof_exit(control, frame);
		printf("finished eof_exit\n");
		return;
	}
	
	//check if content matches expected 1st line of frame
	memcpy( test, &line[0], 14);
	test[14] = '\0';
	if( strcmp(test, "ITEM: TIMESTEP") != 0) {
			printf("error in frame format\n");
			printf("line found was %s\n", test);
			printf("expected line is ITEM: TIMESTEP\n");
			*flag = 0;
			eof_exit(control, frame);
			return;
	}
	
	//so, we assume content is correct now for next frame and begin reading
	//read in timestep
	fgets(line,100,df); //2
    sscanf(line, "%d", &frame->timestep);
    
    //read in number of atoms
    fgets(line,100,df); //3 
	fgets(line,100,df); //4
    sscanf(line, "%d", &i);
    
    printf("timestep is %d\n", frame->timestep);
    
    //check if num_atoms changed
    if(i != frame->num_atoms) {
    	if(frame->num_atoms == 0) {	//do initial allocation
    		if(i != control->num_fg_sites) {
    			if(control->sens_map_flag == 1) {
					if(i != control->num_cg_sites) {
    					printf("num cg sites = %d does not match num atoms %d\n", control->num_cg_sites, i);
    				}
    			} else {
					printf("num fg sites = %d does not match num atoms %d\n", control->num_fg_sites, i);
				}
    		}
			frame->num_atoms = i;
			frame->atoms = malloc(frame->num_atoms * sizeof(ATOM));
			frame->num_observables = control->num_observables;
			for(j = 0; j < frame->num_atoms; j++) frame->atoms[j].mol = 0;
		} else { //number of atoms has changed
			printf("number of atoms has  changed");
			//free all existing data for atoms
			for(j = 0; j < frame->num_atoms; j++) {
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
		for(i=0; i < frame->num_atoms; i++) {
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

    //read in atoms and observables for frame
	fgets(line,100,df);
	(*control->read_function)(frame, df);		
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
	read_frame(control, inframe1, df1, &test1);
	if(test1 == 0) {
		*flag = 0;
		return;	
	}
	
	read_frame(control, inframe2, df2, &test2);
	if(test2 == 0) {
		*flag = 0;
		return;
	}
		
	printf("1:inframe1->num_atoms = %d and inframe2->num_atoms = %d\n", inframe1->num_atoms, inframe2->num_atoms);
	read_logfile(control, lf, &double_val, &test_log);
	control->log_value = double_val;
	
	if(test_log == 0)  {
		*flag = 0;
		return;
	}
	
	read_guess(control, gf, &guess_read);
	if(guess_read == 0) {
		*flag = 0;
		return;
	}
	
	//do test on flags to determine composite (0 in any is a fail)
	*flag = test1 && test2 && test_log && guess_read;	
	printf("2:inframe1->num_atoms = %d and inframe2->num_atoms = %d\n", inframe1->num_atoms, inframe2->num_atoms);
}

  ////////////////////////////////
 //   read_charge_frames	  ///
////////////////////////////////

void read_charge_frames(Controller* control, Frame* inframes, int* flag)
{	//declare variables
	int i;
	int temp = 1;
	
	//read each frame in turn
	for(i = 0; i < control->num_files; i++) {
		printf("read frame %d\n", i);
		read_frame(control, &(inframes[i]), (control->file_point[i]), &temp);
		printf("result is %d\n", temp);
		if(temp == 0) {
			*flag = 0;
			return;
		}
	}
}

  ////////////////////////
 //   read_logfile	  ///
////////////////////////

void read_logfile(Controller* control, FILE* lf, double* out_val, int* flag)
{	//declare varaibles
	int i;
	double time, pe, vdwl, coul, volume;
	char line[100];
	char test[8];
	
	//check if we are at EOF
	if( fgets(line, 100, lf) == NULL ) {
		*flag = 0;
		printf("end of file reached in logfile!\n");
		return;
	}

	//check if we need to skip header info
	if(control->frame == 1) {
		i = 1;
		char test2[5];
		
		while(i == 1) {
			//check line
			memcpy(test2, &line[0], 4);
			test2[4] = '\0';
			if( strcmp(test2, "Step") == 0) {
				i = 0;			
			}
			
			//check if we are at EOF
			if( fgets (line, 100, lf) == NULL ) {
				*flag = 0;
				printf("end of file reached in logfile!\n");
				return;
			}
		}
	}
	
	//see if this is a WARNING line or actual content
	memcpy(test, &line[0], 7);
	test[7] = '\0';
	if( strcmp(test, "WARNING") == 0) {
		fgets(line, 100, lf);
	}
	
	//read actual content for energy value
	sscanf(line, "%lf %lf %lf %lf %lf", &time, &pe, &vdwl, &coul, &volume);
	
	if(control->log_type == 0) {
		*out_val = vdwl;
	} else if(control->log_type == 1) {
		*out_val = coul;
	}
		
	control->timestep = time;
	control->volume = volume;
}

/////////////////////////
///   read_guess	  ///
////////////////////////

void read_guess(Controller* control, FILE* gf, int* flag)
{
	char line[100];
	int i;
	
	if(control->guess_type == 0) {
		//check if we are at EOF
		if( fgets (line, 100, gf) == NULL ) {
			*flag = 0;
			printf("end of file reached in GUESS file!\n");
			return;
		}

		//read guess value
		sscanf(line, "%lf", &control->guess);
	} else if(control->guess_type == 1) {
		//check if we are at EOF
		if( fgets (line, 100, gf) == NULL ) {
			*flag = 0;
			printf("end of file reached in GUESS file!\n");
			return;
		}

		//check if we need to skip header info
		if(control->frame == 1) {
			char test[8];
			char test2[5];
			double time, pe, vdwl, coul, volume;
			i = 1;
		
			while(i == 1) {
				//check line
				memcpy(test2, &line[0], 4);
				test2[4] = '\0';
				if( strcmp(test2, "Step") == 0) i = 0;			
				//read next line
				fgets(line, 100, gf);
			
				//check if we are at EOF
				if( fgets (line, 100, gf) == NULL ) {
					*flag = 0;
					printf("end of file reached in GUESS!\n");
					return;
				}
			}

			//see if this is a WARNING line or actual content
			memcpy(test, &line[0], 7);
			test[7] = '\0';
			if( strcmp(test, "WARNING") == 0) {
				fgets(line, 100, gf);
			}
				
			//read actual content for energy value
			sscanf(line, "%lf %lf %lf %lf %lf", &time, &pe, &vdwl, &coul, &volume);
			if(control->log_type == 0) {
				control->guess = vdwl;
			} else if(control->log_type == 1) {
				control->guess = coul;
			}
		}
	}
}

//////////////////////////////////
///   read_number_in_line	  ///
////////////////////////////////

void read_number_in_line(int num_vals, char* line, int* vals)
{
	switch(num_vals) {
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
	switch(num_vals) {
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
        if (name[i] == '.') {
			pos = i;
        }
    }
    if (pos < 0) {
    	report_traj_input_suffix_error(suffix);
    	}
    if (strcmp(&name[pos + 1], suffix) != 0) 
    	{
    	report_traj_input_suffix_error(suffix);
		}
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
	if(control->frame <= 1) {
		frame->atoms = malloc(1 * sizeof(ATOM));
		frame->atoms[0].observables = malloc(1 * sizeof(double));
	}
}


//////////////////////////////////////////
///		FUNCTION POINTER READ FUNCTION ///
//////////////////////////////////////////

void read0full( Frame* frame, FILE* df)
{
	int i;
	char line[100];
	for(i=0; i < frame->num_atoms; i++) {
		fgets(line,100,df);
		sscanf(line, "%d %d %d %lf %lf %lf %lf %lf", &frame->atoms[i].id, \
    	&frame->atoms[i].mol, &frame->atoms[i].type, &frame->atoms[i].q, &frame->atoms[i].mass, \
   		&frame->atoms[i].x, &frame->atoms[i].y, &frame->atoms[i].z);		
		if( (frame->atoms[i].mol > frame->num_mol) || (frame->num_mol > frame->num_atoms) )  {
			frame->num_mol = frame->atoms[i].mol;
		}
	}
}
void read3full( Frame* frame, FILE* df)
{
	int i;
	char line[100];
	for(i=0; i < frame->num_atoms; i++) {
		fgets(line,100,df);
		sscanf(line, "%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf", &frame->atoms[i].id, \
		&frame->atoms[i].mol, &frame->atoms[i].type, &frame->atoms[i].q, &frame->atoms[i].mass, \
		&frame->atoms[i].x, &frame->atoms[i].y, &frame->atoms[i].z, &frame->atoms[i].observables[0], \ 
		&frame->atoms[i].observables[1], &frame->atoms[i].observables[2]);		
		if( (frame->atoms[i].mol > frame->num_mol) || (frame->num_mol > frame->num_atoms) )  {
			frame->num_mol = frame->atoms[i].mol;
		}
	}
}
void read0min( Frame* frame, FILE* df)
{
	int i;
	char line[100];
	for(i=0; i < frame->num_atoms; i++) {
		fgets(line,100,df);
		sscanf(line, "%lf %lf %lf", \
		&frame->atoms[i].x, &frame->atoms[i].y, &frame->atoms[i].z);		
		if( (frame->atoms[i].mol > frame->num_mol) || (frame->num_mol > frame->num_atoms) ) {
			frame->num_mol = frame->atoms[i].mol;
		}
	}
}		
void read3min( Frame* frame, FILE* df)
{
	int i;
	char line[100];
	for(i=0; i < frame->num_atoms; i++) {
		fgets(line,100,df);
		sscanf(line, "%lf %lf %lf %lf %lf %lf", \
		&frame->atoms[i].x, &frame->atoms[i].y, &frame->atoms[i].z, &frame->atoms[i].observables[0], \ 
		&frame->atoms[i].observables[1], &frame->atoms[i].observables[2]);		
		if( (frame->atoms[i].mol > frame->num_mol) || (frame->num_mol > frame->num_atoms) ) {
			frame->num_mol = frame->atoms[i].mol;
		}
	}
}