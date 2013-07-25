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
void output_charge_log(Controller*);
void create_and_output_bootstrapping_file(Controller*, Frame*, int*, int*);

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
	FILE* of = fopen(outfile, "a+");
	
	//output frame header
	fprintf(of, "ITEM: TIMESTEP\n");
	fprintf(of, "%d\n", outframe->timestep);
	fprintf(of, "ITEM: NUMBER OF ATOMS\n");
	fprintf(of, "%d\n", outframe->num_atoms); //could also consider using outframe->num_mol
	fprintf(of, "ITEM: BOX BOUNDS pp pp pp\n");
	fprintf(of, "%lf %lf\n%lf %lf\n%lf %lf\n", outframe->xmin, outframe->xmax, outframe->ymin, outframe->ymax, outframe->zmin, outframe->zmax);

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

/////////////////////////////////
///   output_frame_minimal	  ///
/////////////////////////////////

void output_frame_minimal(Controller* control, Frame* outframe, char* outfile)
{
	int i, j;
	FILE* of = fopen(outfile, "a+");
	
	//output frame header
	fprintf(of, "ITEM: TIMESTEP\n");
	fprintf(of, "%d\n", outframe->timestep);
	fprintf(of, "ITEM: NUMBER OF ATOMS\n");
	fprintf(of, "%d\n", outframe->num_atoms); //could also consider using outframe->num_mol
	fprintf(of, "ITEM: BOX BOUNDS pp pp pp\n");
	fprintf(of, "%lf %lf\n%lf %lf\n%lf %lf\n", outframe->xmin, outframe->xmax, outframe->ymin, outframe->ymax, outframe->zmin, outframe->zmax);

	switch(outframe->num_observables)
		{	
		case 0:
			fprintf(of, "ITEM: ATOMS x y z \n");
			break;
		
		case 1:
			fprintf(of, "ITEM: ATOMS x y z U \n");
			break;
			
		case 3:
			fprintf(of, "ITEM: ATOMS x y z fx fy fz \n");
			break;
		
		case 4:
			fprintf(of, "ITEM: ATOMS x y z dfx dfy dfz dU \n");
			break;
			
		case 6:
			fprintf(of, "ITEM: ATOMS x y z fx fy fz dfx dfy dfz \n");
			break;
			
		case 7:
			fprintf(of, "ITEM: ATOMS x y z fx fy fz dfx dfy dfz dU \n");
			break;
			
		default:
			printf("output number of observables not supported so giving generic 3 observable header\n");
			fprintf(of, "ITEM: ATOMS x y z fx fy fz\n");
			//for(i = 0; i < outframe->num_atoms; i++)
		}
		
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

	//output frame header
	fprintf(of, "ITEM: TIMESTEP\n");
	fprintf(of, "%d\n", outframe->timestep);
	fprintf(of, "ITEM: NUMBER OF ATOMS\n");
	fprintf(of, "%d\n", outframe->num_atoms); //could also consider using outframe->num_mol
	fprintf(of, "ITEM: BOX BOUNDS pp pp pp\n");
	fprintf(of, "%lf %lf\n%lf %lf\n%lf %lf\n", outframe->xmin, outframe->xmax, outframe->ymin, outframe->ymax, outframe->zmin, outframe->zmax);

	//printf("finished header\n");
	
	switch(outframe->num_observables)
		{	
		case 0:
			fprintf(of, "ITEM: ATOMS x y z \n");
			break;
		
		case 1:
			fprintf(of, "ITEM: ATOMS x y z U \n");
			break;
			
		case 3:
			fprintf(of, "ITEM: ATOMS x y z fx fy fz \n");
			break;
		
		case 4:
			fprintf(of, "ITEM: ATOMS x y z dfx dfy dfz dU \n");
			break;
			
		case 6:
			fprintf(of, "ITEM: ATOMS x y z fx fy fz dfx dfy dfz \n");
			break;
			
		case 7:
			fprintf(of, "ITEM: ATOMS x y z fx fy fz dfx dfy dfz dU \n");
			break;
			
		default:
			printf("output number of observables not supported so giving generic 3 observable header\n");
			fprintf(of, "ITEM: ATOMS x y z fx fy fz\n");
			//for(i = 0; i < outframe->num_atoms; i++)
		}
		
	//printf("writing data\n");
	//printf("num_observables is %d\n", outframe->num_observables);
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
}

/////////////////////////////
///   output_frame_4	  ///
/////////////////////////////

void output_frame_4(Controller* control, Frame* outframe, char* outfile)
{
	int i, j;
	FILE* of = fopen(outfile, "a+");
	
	//output frame header
	fprintf(of, "ITEM: TIMESTEP\n");
	fprintf(of, "%d\n", outframe->timestep);
	fprintf(of, "ITEM: NUMBER OF ATOMS\n");
	fprintf(of, "%d\n", outframe->num_atoms); //could also consider using outframe->num_mol
	fprintf(of, "ITEM: BOX BOUNDS pp pp pp\n");
	fprintf(of, "%lf %lf\n%lf %lf\n%lf %lf\n", outframe->xmin, outframe->xmax, outframe->ymin, outframe->ymax, outframe->zmin, outframe->zmax);

	fprintf(of, "ITEM: ATOMS id mol type q mass x y z U \n");
		
	for( i = 0; i < outframe->num_atoms; i++)
		{
		fprintf(of, "%d %d %d %lf %lf %lf %lf %lf %lf\n", outframe->sites[i].id, outframe->sites[i].mol, \
			outframe->sites[i].type, outframe->sites[i].q, outframe->sites[i].mass, outframe->sites[i].x, \
			outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[3]);
		}
	fclose(of);
}

/////////////////////////////////
///   output_frame_second3	  ///
/////////////////////////////////

void output_frame_second3(Controller* control, Frame* outframe, char* outfile)
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

	fprintf(of, "ITEM: ATOMS id mol type q mass x y z dfx dfy dfz \n");
	
	for( i = 0; i < outframe->num_atoms; i++)
		{
		
		fprintf(of, "%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf\n", outframe->sites[i].id, outframe->sites[i].mol, \
			outframe->sites[i].type, outframe->sites[i].q, outframe->sites[i].mass, outframe->sites[i].x, \
			outframe->sites[i].y, outframe->sites[i].z, outframe->sites[i].observables[3], \
			outframe->sites[i].observables[4], outframe->sites[i].observables[5]);
		}
		
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
	//declare variables
	int i;
	
	//printf("output frames\n");
	for(i = 0; i < control->num_charges; i++)
		{
		//printf("output %d\n", i);
		output_frame_minimal_charge(control, &outframe[i], control->outfile[i]);
		
		//printf("output %d finished\n", i);
		}
}

////////////////////////////
///   output_charge_log  ///
////////////////////////////

void output_charge_log(Controller* control)
{
	//declare variables
	int i;
	
	if(control->frame <= 1)
		{
		for(i = 0; i < control->num_charges; i++) 
			{
			fprintf(control->outfile[i], "Step PotEng E_vdwl E_coul Volume\n");
			}
		}
	
	for(i = 0; i < control->num_charges; i++)
		{
		fprintf(control->outfile[i], "%lf\t%lf\t%lf\t%lf\t%lf\n", control->timestep, control->guesses[i], 0.0, control->guesses[i], control->volume);
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
		//printf("output_frame #%d as value frames[ %d ]\n", i, order[i]);
		output_frame(control, &(frames[ order[i] ]), filename);
		}
	//printf("finished create_and_output_bootstrapping_file\n");
}