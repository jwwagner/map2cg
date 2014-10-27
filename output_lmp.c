//output_lmp.c

#include "headers.h"

void output_frame(Controller*, Frame*, char*);
void output_frame_all(Controller*, Frame*, char*);
void output_frame_minimal(Controller*, Frame*, char*);
void output_frame_minimal_charge(Controller*, Frame*, FILE*);
void output_frame_4(Controller*, Frame*, char*);
void output_frame_second3(Controller*, Frame*, char*);
void output_topology(Controller*, Frame*);
void output_force_file(Controller*, int, double*, double*, double*, char*);
void output_charge_frames(Controller*, Frame*);
void create_and_output_bootstrapping_file(Controller*, Frame*, int*, int*);
//pointer output functions
void out0full(Frame*, FILE*);
void out1full(Frame*, FILE*);
void out2full(Frame*, FILE*);
void out3full(Frame*, FILE*);
void out4full(Frame*, FILE*);
void out5full(Frame*, FILE*);
void out6full(Frame*, FILE*);
void out7full(Frame*, FILE*);
void out8full(Frame*, FILE*);
void out9full(Frame*, FILE*);
void out0min(Frame*, FILE*);
void out1min(Frame*, FILE*);
void out2min(Frame*, FILE*);
void out3min(Frame*, FILE*);
void out4min(Frame*, FILE*);
void out5min(Frame*, FILE*);
void out6min(Frame*, FILE*);
void out7min(Frame*, FILE*);
void out8min(Frame*, FILE*);
void out9min(Frame*, FILE*);
void out0minid(Frame*, FILE*);
void out1minid(Frame*, FILE*);
void out2minid(Frame*, FILE*);
void out3minid(Frame*, FILE*);
void out4minid(Frame*, FILE*);
void out5minid(Frame*, FILE*);
void out6minid(Frame*, FILE*);
void out7minid(Frame*, FILE*);
void out8minid(Frame*, FILE*);
void out9minid(Frame*, FILE*);
void default_func(FILE* of);
void header0(FILE*);
void header1(FILE*);
void header3(FILE*);
void header4(FILE*);
void header6(FILE*);
void header7(FILE*);
void header2(FILE*);
void header5(FILE*);
void header8(FILE*);
void header9(FILE*);
void header0id(FILE*);
void header1id(FILE*);
void header3id(FILE*);
void header4id(FILE*);
void header6id(FILE*);
void header7id(FILE*);
void header2id(FILE*);
void header5id(FILE*);
void header8id(FILE*);
void header9id(FILE*);
void header0full(FILE*);
void header1full(FILE*);
void header3full(FILE*);
void header4full(FILE*);
void header6full(FILE*);
void header7full(FILE*);
void generic_frame_header(Frame*, FILE*);
//////////////////////////
///   output_frame	  ///
/////////////////////////
void output_frame(Controller* control, Frame* outframe, char* outfile)
{
	switch(control->output_flag)
		{
		case 1:
			output_frame_minimal(control, outframe, outfile);
			break;
		case 2:
			output_frame_4(control, outframe, outfile);
			break;
		case 3:
			output_frame_second3(control, outframe, outfile);
		default:
		case 0:
			output_frame_all(control, outframe, outfile);
			break;
		}
}

//////////////////////////////
///   output_frame_all	  ///
////////////////////////////

void output_frame_all(Controller* control, Frame* outframe, char* outfile)
{
	int i, j;
	if(outframe->sites[0].type == -1)
		{
		for( i = 0; i < outframe->num_atoms; i++)
			{	
			printf("%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf\n", outframe->sites[i].id, outframe->sites[i].mol, \
			outframe->sites[i].type, outframe->sites[i].q, outframe->sites[i].mass, outframe->sites[i].x, \
			outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0], \
			outframe->sites[i].observables[1], outframe->sites[i].observables[2]);
			}
		}
	FILE* of = fopen(outfile, "a+");
	generic_frame_header(outframe, of); 	//output frame header
	(*control->header_function)(of);
	(*control->output_function)(outframe, of);
	fclose(of);
}

/////////////////////////////////
///   output_frame_minimal	  ///
/////////////////////////////////

void output_frame_minimal(Controller* control, Frame* outframe, char* outfile)
{
	int i, j;
	FILE* of = fopen(outfile, "a+");
	
	generic_frame_header(outframe, of); 	//output frame header
	(*control->header_function)(of);
		
	for( i = 0; i < outframe->num_atoms; i++)
		{
		switch(outframe->num_observables)
			{
			case 0:
				fprintf(of, "%lf %lf %lf\n", outframe->sites[i].x, \
				outframe->sites[i].y, outframe->sites[i].z);
				break;
			
			case 1:
				fprintf(of, "%lf %lf %lf %lf\n", outframe->sites[i].x, \
				outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0]);
				break;
				
			case 2:
				fprintf(of, "%lf %lf %lf %lf %lf\n", outframe->sites[i].x, \
				outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0], \
				outframe->sites[i].observables[1]);
				break;

			case 3:
				fprintf(of, "%lf %lf %lf %lf %lf %lf\n", outframe->sites[i].x, \
				outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0], \
				outframe->sites[i].observables[1], outframe->sites[i].observables[2]);
				break;
				
			case 4:
				fprintf(of, "%lf %lf %lf %lf %lf %lf %lf\n", outframe->sites[i].x, \
				outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0], \
				outframe->sites[i].observables[1], outframe->sites[i].observables[2], \
				outframe->sites[i].observables[3]);
				break;
			
			case 5:
				fprintf(of, "%lf %lf %lf %lf %lf %lf %lf %lf\n", outframe->sites[i].x, \
				outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0], \
				outframe->sites[i].observables[1], outframe->sites[i].observables[2], \
				outframe->sites[i].observables[3], outframe->sites[i].observables[4]);
				break;
				
			case 6:
				fprintf(of, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", outframe->sites[i].x, \
				outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0], \
				outframe->sites[i].observables[1], outframe->sites[i].observables[2], \
				outframe->sites[i].observables[3], outframe->sites[i].observables[4], \
				outframe->sites[i].observables[5]);
				break;
				
			case 7:
				fprintf(of, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", outframe->sites[i].x, \
				outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0], \
				outframe->sites[i].observables[1], outframe->sites[i].observables[2], \
				outframe->sites[i].observables[3], outframe->sites[i].observables[4], \
				outframe->sites[i].observables[5], outframe->sites[i].observables[6]);
				break;
				
			case 8:
				fprintf(of, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", outframe->sites[i].x, \
				outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0], \
				outframe->sites[i].observables[1], outframe->sites[i].observables[2], \
				outframe->sites[i].observables[3], outframe->sites[i].observables[4], \
				outframe->sites[i].observables[5], outframe->sites[i].observables[6], \
				outframe->sites[i].observables[7]);
				break;
				
			case 9:
				fprintf(of, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", outframe->sites[i].x, \
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

/////////////////////////////////////
///   output_frame_minimal_charge ///
/////////////////////////////////////

void output_frame_minimal_charge(Controller* control, Frame* outframe, FILE* of)
{
	int i, j;
	generic_frame_header(outframe, of);	//output frame header
	(*control->header_function)(of);		
	(*control->output_function)(outframe, of);
}

/////////////////////////////
///   output_frame_4	  ///
/////////////////////////////

void output_frame_4(Controller* control, Frame* outframe, char* outfile)
{
	int i, j;
	FILE* of = fopen(outfile, "a+");
	generic_frame_header(outframe, of); //output frame header
	header4full(of);
	out4full(outframe, of);	
	fclose(of);
}

/////////////////////////////////
///   output_frame_second3	  ///
/////////////////////////////////

void output_frame_second3(Controller* control, Frame* outframe, char* outfile)
{
	int i, j;
	FILE* of = fopen(outfile, "a+");
	generic_frame_header(outframe, of); //output frame header
	header6full(of);
	out6full(outframe, of);
	fclose(of);
}

//////////////////////////
///   output_topology  ///
/////////////////////////
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

/////////////////////////////
///   output_force_file  ///
////////////////////////////

void output_force_file(Controller* control, int num, double* distance, double* potential, double* force, char* file)
{
	//declare variables
	int i;
	FILE* fp = fopen(file, "w+");
	
	//write header information
	printf("write header information\n");
	fprintf(fp, "# Header information on force file\n");
	fprintf(fp, "\n");
	fprintf(fp, "%s\n", control->files.guess);
	
	printf("checkpoint #1\n");
	printf("num is %d\n", num);
	printf("distance[0] is %lf\n", distance[0]);
	printf("distance[n-1] is %lf\n", distance[num-1] );
	
	fprintf(fp, "N %d R %lf %lf\n", num, distance[0], distance[num-1]);
	fprintf(fp, "\n");
	
	printf("write body information\n");
	
	//write body information
	for(i = 0; i < num; i++)
		{
		fprintf(fp, "%d %lf %lf %lf\n", (i+1), distance[i], potential[i], force[i]);
		}
	
	printf("close output file\n");
	//close file
	fclose(fp);
}

///////////////////////////////
///   output_charge_frames  ///
///////////////////////////////

void output_charge_frames(Controller* control, Frame* outframe)
{
	int i=0;	
	for(i = 0; i < control->num_outfile; i++)
		{
		printf("output %d\n", i);
		output_frame_minimal_charge(control, &outframe[i], control->outfile[i]);
		printf("output %d finished\n", i);
		}
}

///////////////////////////////////////////////
///   create_and_output_bootstrapping_file  ///
///////////////////////////////////////////////

void create_and_output_bootstrapping_file(Controller* control, Frame* frames, int* order, int* count)
{
	//declare variables
	int i;
	char filename[35] = "";
	char suffix[] = "map.dat"; 
	char num[5];
	FILE* fp;

	//determine name of output file 
	snprintf(filename, sizeof filename, "%s%d%s", control->name, *count, suffix);
	printf("open output file named %s\n", filename);
	
	//toggle output file to restart
	fp = fopen(filename, "w+");
	fclose(fp);
	
	//write each frame in order to file
	for(i=0; i < control->num_frames; i++)
		{
		output_frame(control, &(frames[ order[i] ]), filename);
		}
}

////////////////////////////////////////////
/////	output_pointer_functions		////
///////////////////////////////////////////
void out0full(Frame* outframe, FILE* of)
{
	int i;
	for( i = 0; i < outframe->num_atoms; i++)
		{	
		fprintf(of, "%d %d %d %lf %lf %lf %lf %lf\n", outframe->sites[i].id, outframe->sites[i].mol, \
		outframe->sites[i].type, outframe->sites[i].q, outframe->sites[i].mass, outframe->sites[i].x, \
		outframe->sites[i].y, outframe->sites[i].z);
		}
}
void out1full(Frame* outframe, FILE* of)
{
	int i;
	for( i = 0; i < outframe->num_atoms; i++)
		{	
		fprintf(of, "%d %d %d %lf %lf %lf %lf %lf %lf\n", outframe->sites[i].id, outframe->sites[i].mol, \
		outframe->sites[i].type, outframe->sites[i].q, outframe->sites[i].mass, outframe->sites[i].x, \
		outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0]);
		}
}
void out2full(Frame* outframe, FILE* of)
{
	int i;
	for( i = 0; i < outframe->num_atoms; i++)
		{	
		fprintf(of, "%d %d %d %lf %lf %lf %lf %lf %lf %lf\n", outframe->sites[i].id, outframe->sites[i].mol, \
		outframe->sites[i].type, outframe->sites[i].q, outframe->sites[i].mass, outframe->sites[i].x, \
		outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0], \
		outframe->sites[i].observables[1]);
		}
}
void out3full(Frame* outframe, FILE* of)
{
	int i;
	for( i = 0; i < outframe->num_atoms; i++)
		{	
		fprintf(of, "%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf\n", outframe->sites[i].id, outframe->sites[i].mol, \
		outframe->sites[i].type, outframe->sites[i].q, outframe->sites[i].mass, outframe->sites[i].x, \
		outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0], \
		outframe->sites[i].observables[1], outframe->sites[i].observables[2]);
		}
}
void out4full(Frame* outframe, FILE* of)
{
	int i;
	for( i = 0; i < outframe->num_atoms; i++)
		{	
		fprintf(of, "%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", outframe->sites[i].id, outframe->sites[i].mol, \
		outframe->sites[i].type, outframe->sites[i].q, outframe->sites[i].mass, outframe->sites[i].x, \
		outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0], \
		outframe->sites[i].observables[1], outframe->sites[i].observables[2], \
		outframe->sites[i].observables[3]);
		}
}	
void out5full(Frame* outframe, FILE* of)
{
	int i;
	for( i = 0; i < outframe->num_atoms; i++)
		{	
		fprintf(of, "%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", outframe->sites[i].id, outframe->sites[i].mol, \
		outframe->sites[i].type, outframe->sites[i].q, outframe->sites[i].mass, outframe->sites[i].x, \
		outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0], \
		outframe->sites[i].observables[1], outframe->sites[i].observables[2], \
		outframe->sites[i].observables[3], outframe->sites[i].observables[4]);
		}
}
void out6full(Frame* outframe, FILE* of)
{
	int i;
	for( i = 0; i < outframe->num_atoms; i++)
		{	
		fprintf(of, "%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", outframe->sites[i].id, outframe->sites[i].mol, \
		outframe->sites[i].type, outframe->sites[i].q, outframe->sites[i].mass, outframe->sites[i].x, \
		outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0], \
		outframe->sites[i].observables[1], outframe->sites[i].observables[2], \
		outframe->sites[i].observables[3], outframe->sites[i].observables[4], \
		outframe->sites[i].observables[5]);
		}
}
void out7full(Frame* outframe, FILE* of)
{
	int i;
	for( i = 0; i < outframe->num_atoms; i++)
		{	
		fprintf(of, "%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", outframe->sites[i].id, outframe->sites[i].mol, \
		outframe->sites[i].type, outframe->sites[i].q, outframe->sites[i].mass, outframe->sites[i].x, \
		outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0], \
		outframe->sites[i].observables[1], outframe->sites[i].observables[2], \
		outframe->sites[i].observables[3], outframe->sites[i].observables[4], \
		outframe->sites[i].observables[5], outframe->sites[i].observables[6]);
		}
}
void out8full(Frame* outframe, FILE* of)
{
	int i;
	for( i = 0; i < outframe->num_atoms; i++)
		{	
		fprintf(of, "%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", outframe->sites[i].id, outframe->sites[i].mol, \
		outframe->sites[i].type, outframe->sites[i].q, outframe->sites[i].mass, outframe->sites[i].x, \
		outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0], \
		outframe->sites[i].observables[1], outframe->sites[i].observables[2], \
		outframe->sites[i].observables[3], outframe->sites[i].observables[4], \
		outframe->sites[i].observables[5], outframe->sites[i].observables[6], \
		outframe->sites[i].observables[7]);
		}
}
void out9full(Frame* outframe, FILE* of)
{
	int i;
	for( i = 0; i < outframe->num_atoms; i++)
		{	
		fprintf(of, "%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", outframe->sites[i].id, outframe->sites[i].mol, \
		outframe->sites[i].type, outframe->sites[i].q, outframe->sites[i].mass, outframe->sites[i].x, \
		outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0], \
		outframe->sites[i].observables[1], outframe->sites[i].observables[2], \
		outframe->sites[i].observables[3], outframe->sites[i].observables[4], \
		outframe->sites[i].observables[5], outframe->sites[i].observables[6], \
		outframe->sites[i].observables[7], outframe->sites[i].observables[8]);
		}
}
void out0min(Frame* outframe, FILE* of)
{
	int i;
	for( i = 0; i < outframe->num_atoms; i++)
		{	
		fprintf(of, "%d %lf %lf %lf %lf %lf\n", outframe->sites[i].id, \
		outframe->sites[i].q, outframe->sites[i].mass, outframe->sites[i].x, \
		outframe->sites[i].y, outframe->sites[i].z);
		}
}
void out1min(Frame* outframe, FILE* of)
{
	int i;
	for( i = 0; i < outframe->num_atoms; i++)
		{	
		fprintf(of, "%d %lf %lf %lf %lf %lf %lf\n", outframe->sites[i].id, \
		outframe->sites[i].q, outframe->sites[i].mass, outframe->sites[i].x, \
		outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0]);
		}
}
void out2min(Frame* outframe, FILE* of)
{
	int i;
	for( i = 0; i < outframe->num_atoms; i++)
		{	
		fprintf(of, "%d %lf %lf %lf %lf %lf %lf %lf\n", outframe->sites[i].id, \
		outframe->sites[i].q, outframe->sites[i].mass, outframe->sites[i].x, \
		outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0], \
		outframe->sites[i].observables[1]);
		}
}
void out3min(Frame* outframe, FILE* of)
{
	int i;
	for( i = 0; i < outframe->num_atoms; i++)
		{	
		fprintf(of, "%d %lf %lf %lf %lf %lf %lf %lf %lf\n", outframe->sites[i].id, \
		outframe->sites[i].q, outframe->sites[i].mass, outframe->sites[i].x, \
		outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0], \
		outframe->sites[i].observables[1], outframe->sites[i].observables[2]);
		}
}
void out4min(Frame* outframe, FILE* of)
{
	int i;
	for( i = 0; i < outframe->num_atoms; i++)
		{	
		fprintf(of, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", outframe->sites[i].id, \
		outframe->sites[i].q, outframe->sites[i].mass, outframe->sites[i].x, \
		outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0], \
		outframe->sites[i].observables[1], outframe->sites[i].observables[2], \
		outframe->sites[i].observables[3]);
		}
}	
void out5min(Frame* outframe, FILE* of)
{
	int i;
	for( i = 0; i < outframe->num_atoms; i++)
		{	
		fprintf(of, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", outframe->sites[i].id, \
		outframe->sites[i].q, outframe->sites[i].mass, outframe->sites[i].x, \
		outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0], \
		outframe->sites[i].observables[1], outframe->sites[i].observables[2], \
		outframe->sites[i].observables[3], outframe->sites[i].observables[4]);
		}
}
void out6min(Frame* outframe, FILE* of)
{
	int i;
	for( i = 0; i < outframe->num_atoms; i++)
		{	
		fprintf(of, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", outframe->sites[i].id, \
		outframe->sites[i].q, outframe->sites[i].mass, outframe->sites[i].x, \
		outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0], \
		outframe->sites[i].observables[1], outframe->sites[i].observables[2], \
		outframe->sites[i].observables[3], outframe->sites[i].observables[4], \
		outframe->sites[i].observables[5]);
		}
}
void out7min(Frame* outframe, FILE* of)
{
	int i;
	for( i = 0; i < outframe->num_atoms; i++)
		{	
		fprintf(of, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", outframe->sites[i].id, \
		outframe->sites[i].q, outframe->sites[i].mass, outframe->sites[i].x, \
		outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0], \
		outframe->sites[i].observables[1], outframe->sites[i].observables[2], \
		outframe->sites[i].observables[3], outframe->sites[i].observables[4], \
		outframe->sites[i].observables[5], outframe->sites[i].observables[6]);
		}
}
void out8min(Frame* outframe, FILE* of)
{
	int i;
	for( i = 0; i < outframe->num_atoms; i++)
		{	
		fprintf(of, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", outframe->sites[i].id, \
		outframe->sites[i].q, outframe->sites[i].mass, outframe->sites[i].x, \
		outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0], \
		outframe->sites[i].observables[1], outframe->sites[i].observables[2], \
		outframe->sites[i].observables[3], outframe->sites[i].observables[4], \
		outframe->sites[i].observables[5], outframe->sites[i].observables[6], \
		outframe->sites[i].observables[7]);
		}
}
void out9min(Frame* outframe, FILE* of)
{
	int i;
	for( i = 0; i < outframe->num_atoms; i++)
		{	
		fprintf(of, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", outframe->sites[i].id, \
		outframe->sites[i].q, outframe->sites[i].mass, outframe->sites[i].x, \
		outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0], \
		outframe->sites[i].observables[1], outframe->sites[i].observables[2], \
		outframe->sites[i].observables[3], outframe->sites[i].observables[4], \
		outframe->sites[i].observables[5], outframe->sites[i].observables[6], \
		outframe->sites[i].observables[7], outframe->sites[i].observables[8]);
		}
}
void out0minid(Frame* outframe, FILE* of)
{
	int i;
	for( i = 0; i < outframe->num_atoms; i++)
		{	
		fprintf(of, "%d %lf %lf %lf\n", outframe->sites[i].id, \
		outframe->sites[i].x, \
		outframe->sites[i].y, outframe->sites[i].z);
		}
}
void out1minid(Frame* outframe, FILE* of)
{
	int i;
	for( i = 0; i < outframe->num_atoms; i++)
		{	
		fprintf(of, "%d %lf %lf %lf %lf\n", outframe->sites[i].id, \
		outframe->sites[i].x, \
		outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0]);
		}
}
void out2minid(Frame* outframe, FILE* of)
{
	int i;
	for( i = 0; i < outframe->num_atoms; i++)
		{	
		fprintf(of, "%d %lf %lf %lf %lf %lf\n", outframe->sites[i].id, \
		outframe->sites[i].x, \
		outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0], \
		outframe->sites[i].observables[1]);
		}
}
void out3minid(Frame* outframe, FILE* of)
{
	//printf("out3minid\n");
	int i;
	for( i = 0; i < outframe->num_atoms; i++)
		{	
		fprintf(of, "%d %lf %lf %lf %lf %lf %lf\n", outframe->sites[i].id, \
		outframe->sites[i].x, \
		outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0], \
		outframe->sites[i].observables[1], outframe->sites[i].observables[2]);
		}
}
void out4minid(Frame* outframe, FILE* of)
{
	int i;
	for( i = 0; i < outframe->num_atoms; i++)
		{	
		fprintf(of, "%d %lf %lf %lf %lf %lf %lf %lf\n", outframe->sites[i].id, \
		outframe->sites[i].x, \
		outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0], \
		outframe->sites[i].observables[1], outframe->sites[i].observables[2], \
		outframe->sites[i].observables[3]);
		}
}	
void out5minid(Frame* outframe, FILE* of)
{
	int i;
	for( i = 0; i < outframe->num_atoms; i++)
		{	
		fprintf(of, "%d %lf %lf %lf %lf %lf %lf %lf %lf\n", outframe->sites[i].id, \
		outframe->sites[i].x, \
		outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0], \
		outframe->sites[i].observables[1], outframe->sites[i].observables[2], \
		outframe->sites[i].observables[3], outframe->sites[i].observables[4]);
		}
}
void out6minid(Frame* outframe, FILE* of)
{
	int i;
	for( i = 0; i < outframe->num_atoms; i++)
		{	
		fprintf(of, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", outframe->sites[i].id, \
		outframe->sites[i].x, \
		outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0], \
		outframe->sites[i].observables[1], outframe->sites[i].observables[2], \
		outframe->sites[i].observables[3], outframe->sites[i].observables[4], \
		outframe->sites[i].observables[5]);
		}
}
void out7minid(Frame* outframe, FILE* of)
{
	int i;
	for( i = 0; i < outframe->num_atoms; i++)
		{	
		fprintf(of, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", outframe->sites[i].id, \
		outframe->sites[i].x, \
		outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0], \
		outframe->sites[i].observables[1], outframe->sites[i].observables[2], \
		outframe->sites[i].observables[3], outframe->sites[i].observables[4], \
		outframe->sites[i].observables[5], outframe->sites[i].observables[6]);
		}
}
void out8minid(Frame* outframe, FILE* of)
{
	int i;
	for( i = 0; i < outframe->num_atoms; i++)
		{	
		fprintf(of, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", outframe->sites[i].id, \
		outframe->sites[i].x, \
		outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0], \
		outframe->sites[i].observables[1], outframe->sites[i].observables[2], \
		outframe->sites[i].observables[3], outframe->sites[i].observables[4], \
		outframe->sites[i].observables[5], outframe->sites[i].observables[6], \
		outframe->sites[i].observables[7]);
		}
}
void out9minid(Frame* outframe, FILE* of)
{
	int i;
	for( i = 0; i < outframe->num_atoms; i++)
		{	
		fprintf(of, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", outframe->sites[i].id, \
		outframe->sites[i].x, \
		outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0], \
		outframe->sites[i].observables[1], outframe->sites[i].observables[2], \
		outframe->sites[i].observables[3], outframe->sites[i].observables[4], \
		outframe->sites[i].observables[5], outframe->sites[i].observables[6], \
		outframe->sites[i].observables[7], outframe->sites[i].observables[8]);
		}
}
void default_func(FILE* of)
{
	printf("output number of observables not supported so giving generic 3 observable header\n");
	fprintf(of, "ITEM: ATOMS id mol type q mass x y z fx fy fz\n");
}
void header0(FILE* of)
{
	fprintf(of, "ITEM: ATOMS x y z \n");
}
void header1(FILE* of)
{
	fprintf(of, "ITEM: ATOMS x y z U \n");
}
void header3(FILE* of)
{
	fprintf(of, "ITEM: ATOMS x y z fx fy fz \n");
}
void header4(FILE* of)
{
	fprintf(of, "ITEM: ATOMS x y z dfx dfy dfz dU \n");
}
void header6(FILE* of)
{
	fprintf(of, "ITEM: ATOMS x y z fx fy fz dfx dfy dfz \n");
}
void header7(FILE* of)
{
	fprintf(of, "ITEM: ATOMS x y z fx fy fz dfx dfy dfz dU \n");
}
void header2(FILE* of)
{
	fprintf(of, "ITEM: ATOMS x y z fx fy \n");
}
void header5(FILE* of)
{
	fprintf(of, "ITEM: ATOMS x y z U O \n");
}
void header8(FILE* of)
{
	fprintf(of, "ITEM: ATOMS x y z fx fy fz dfx dfy dfz U O \n");
}
void header9(FILE* of)
{
	fprintf(of, "ITEM: ATOMS x y z fx fy fz dfx dfy dfz qfx qfy qfz \n");
}
void header0id(FILE* of)
{
	fprintf(of, "ITEM: ATOMS x y z \n");
}
void header1id(FILE* of)
{
	fprintf(of, "ITEM: ATOMS x y z U \n");
}
void header3id(FILE* of)
{
	fprintf(of, "ITEM: ATOMS x y z fx fy fz \n");
}
void header4id(FILE* of)
{
	fprintf(of, "ITEM: ATOMS x y z dfx dfy dfz dU \n");
}
void header6id(FILE* of)
{
	fprintf(of, "ITEM: ATOMS x y z fx fy fz dfx dfy dfz \n");
}
void header7id(FILE* of)
{
	fprintf(of, "ITEM: ATOMS x y z fx fy fz dfx dfy dfz dU \n");
}
void header2id(FILE* of)
{
	fprintf(of, "ITEM: ATOMS x y z fx fy \n");
}
void header5id(FILE* of)
{
	fprintf(of, "ITEM: ATOMS x y z U O \n");
}
void header8id(FILE* of)
{
	fprintf(of, "ITEM: ATOMS x y z fx fy fz dfx dfy dfz U O \n");
}
void header9id(FILE* of)
{
	fprintf(of, "ITEM: ATOMS x y z fx fy fz dfx dfy dfz qfx qfy qfz \n");
}
void header0full(FILE* of)
{
	fprintf(of, "ITEM: ATOMS id mol type q mass x y z \n");
}
void header1full(FILE* of)
{
	fprintf(of, "ITEM: ATOMS id mol type q mass x y z U \n");
}
void header3full(FILE* of)
{
	fprintf(of, "ITEM: ATOMS id mol type q mass x y z fx fy fz \n");
}
void header4full(FILE* of)
{
	fprintf(of, "ITEM: ATOMS id mol type q mass x y z dfx dfy dfz dU \n");
}
void header6full(FILE* of)
{
	fprintf(of, "ITEM: ATOMS id mol type q mass x y z fx fy fz dfx dfy dfz \n");
}
void header7full(FILE* of)
{
	fprintf(of, "ITEM: ATOMS id mol type q mass x y z fx fy fz dfx dfy dfz dU \n");
}
void generic_frame_header(Frame* outframe, FILE* of)
{
	fprintf(of, "ITEM: TIMESTEP\n");
	fprintf(of, "%d\n", outframe->timestep);
	fprintf(of, "ITEM: NUMBER OF ATOMS\n");
	fprintf(of, "%d\n", outframe->num_atoms); //could also consider using outframe->num_mol
	fprintf(of, "ITEM: BOX BOUNDS pp pp pp\n");
	fprintf(of, "%lf %lf\n%lf %lf\n%lf %lf\n", outframe->xmin, outframe->xmax, outframe->ymin, outframe->ymax, outframe->zmin, outframe->zmax);
}