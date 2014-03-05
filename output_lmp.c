//output_lmp.c

#include "headers.h"

void output_frame(Controller*, Frame*, char*);
void output_frame_all(Controller*, Frame*, char*);
void output_frame_minimal(Controller*, Frame*, char*);
void output_frame_minimal_charge(Controller*, Frame*, FILE*);
void output_topology(Controller*, Frame*);
void output_charge_frames(Controller*, Frame*);
//pointer output functions
void out0full(Frame*, FILE*);
void out3full(Frame*, FILE*);
void default_func(FILE* of);
void header0(FILE*);
void header3(FILE*);
void header0full(FILE*);
void header3full(FILE*);
void generic_frame_header(Frame*, FILE*);
//////////////////////////
///   output_frame	  ///
/////////////////////////
void output_frame(Controller* control, Frame* outframe, char* outfile)
{
	switch(control->output_flag) {
		case 1:
			output_frame_minimal(control, outframe, outfile);
			break;
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
	if(outframe->sites[0].type == -1) {
		for( i = 0; i < outframe->num_atoms; i++) {	
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
		
	for( i = 0; i < outframe->num_atoms; i++) {
		switch(outframe->num_observables) {
			case 0:
				fprintf(of, "%lf %lf %lf\n", outframe->sites[i].x, \
				outframe->sites[i].y, outframe->sites[i].z);
				break;
			case 3:
				fprintf(of, "%lf %lf %lf %lf %lf %lf\n", outframe->sites[i].x, \
				outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0], \
				outframe->sites[i].observables[1], outframe->sites[i].observables[2]);
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

//////////////////////////
///   output_topology  ///
/////////////////////////
void output_topology(Controller* control, Frame* outframe)
{
	FILE* fp = fopen("top.in", "w+");
	int i, j;
	
	fprintf(fp, "cgsites %d\n", outframe->num_atoms);
	fprintf(fp, "cgtypes %d\n", outframe->type_count);
	for(i = 0; i < outframe->type_count; i++) {
		fprintf(fp, "%d\n", i+1);
	}
	
	fprintf(fp, "moltypes %d\n", outframe->type_count);
	for(i = 0; i < outframe->type_count; i++) {
		fprintf(fp, "mol %d %d\n", 1, 3); 
	}
	
	fprintf(fp, "sitetypes\n");
	for(i = 0; i < 1; i++) {
		fprintf(fp, "%d\n", i+1);
	}
	
	fprintf(fp, "bonds %d\n", 0);
	
	fprintf(fp, "system %d\n", outframe->type_count); 
	for(i = 0; i < outframe->type_count; i++) {
		fprintf(fp, "%d %d\n", i+1, outframe->type_num[i]);
	}
	fprintf(fp, "\n");	
	fclose(fp);
}

///////////////////////////////
///   output_charge_frames  ///
///////////////////////////////

void output_charge_frames(Controller* control, Frame* outframe)
{
	int i=0;	
	for(i = 0; i < control->num_outfile; i++) {
		printf("output %d\n", i);
		output_frame_minimal_charge(control, &outframe[i], control->outfile[i]);
		printf("output %d finished\n", i);
	}
}

////////////////////////////////////////////
/////	output_pointer_functions		////
///////////////////////////////////////////
void out0full(Frame* outframe, FILE* of)
{
	int i;
	for( i = 0; i < outframe->num_atoms; i++) {	
		fprintf(of, "%d %d %d %lf %lf %lf %lf %lf\n", outframe->sites[i].id, outframe->sites[i].mol, \
		outframe->sites[i].type, outframe->sites[i].q, outframe->sites[i].mass, outframe->sites[i].x, \
		outframe->sites[i].y, outframe->sites[i].z);
	}
}
void out3full(Frame* outframe, FILE* of)
{
	int i;
	for( i = 0; i < outframe->num_atoms; i++) {	
		fprintf(of, "%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf\n", outframe->sites[i].id, outframe->sites[i].mol, \
		outframe->sites[i].type, outframe->sites[i].q, outframe->sites[i].mass, outframe->sites[i].x, \
		outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[0], \
		outframe->sites[i].observables[1], outframe->sites[i].observables[2]);
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
void header3(FILE* of)
{
	fprintf(of, "ITEM: ATOMS x y z fx fy fz \n");
}
void header0full(FILE* of)
{
	fprintf(of, "ITEM: ATOMS id mol type q mass x y z \n");
}
void header3full(FILE* of)
{
	fprintf(of, "ITEM: ATOMS id mol type q mass x y z fx fy fz \n");
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