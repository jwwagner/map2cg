//process_frame.c

#include "headers.h"

//sensitivity functions
void process_minimal_frame(Controller*, Frame*, Frame*);
void process_no_map_frames_and_log(Controller*, Frame*, Frame*, Frame*);
void combine_sensitivity_data(Controller*, Frame*, Frame*, Frame*);
void combine_cg_sensitivity_data(Controller*, Frame*, Frame*, Frame*);
void process_ij_charge_frames(Controller*, Frame*, Frame*);
void process_ii_charge_frames(Controller*, Frame*, Frame*);
//free functions
void free_sensitivity_intermediates(Controller*, Frame*, Frame*);
void free_charge_intermediates(Controller*, Frame*);

//////////////////////////////////
///   process_minimal_frame	  ///
////////////////////////////////

void process_minimal_frame(Controller* control, Frame* inframe, Frame* outframe)
{
	//declare variables
	int i, j;
	int l = 0;
	int type, type_spot;
	
	//PREPARE INFRAME FORMAT: determine type and type_num information missing by not doing in-house mapping
	if(control->frame == 1) {
		outframe->type = malloc(control->num_cg_types * sizeof(int));
		outframe->type_num = malloc(control->num_cg_types * sizeof(int));
	}

	//populate general information in outframes
	outframe->type_count = 0;
	outframe->xmin = inframe->xmin;
	outframe->xmax = inframe->xmax;
	outframe->ymin = inframe->ymin;
	outframe->ymax = inframe->ymax;
	outframe->zmin = inframe->zmin;
	outframe->zmax = inframe->zmax;
	outframe->timestep = inframe->timestep;
	outframe->num_observables = inframe->num_observables;	

	//copy types
	outframe->type[0] = 1;
	outframe->type_num[0] = inframe->num_atoms;
	
	//allocate site space
	if(outframe->num_atoms != control->num_cg_sites) {
		outframe->num_atoms = control->num_cg_sites;
		outframe->sites = malloc(outframe->num_atoms * sizeof(SITE));		
		for(i=0; i < outframe->num_atoms; i++) {
			outframe->sites[i].observables = malloc(outframe->num_observables * sizeof(double));
		}
	}
	for(j = 0; j < outframe->num_atoms; j++) {
		outframe->sites[j].num_in_site = 1;
		outframe->sites[j].id = j;
		outframe->sites[j].mol = j;
		outframe->sites[j].type = 1;
		outframe->sites[j].mass = control->mass_full;
		outframe->sites[j].q = control->charge_full;
		outframe->sites[j].x = inframe->atoms[j].x;
		outframe->sites[j].y = inframe->atoms[j].y;
		outframe->sites[j].z = inframe->atoms[j].z;
		for(l = 0; l < outframe->num_observables; l++) {
			outframe->sites[j].observables[l] = inframe->atoms[j].observables[l];
		}
	}
}

////////////////////////////////////////
///   process_no_map_frames_and_log  ///
////////////////////////////////////////

void process_no_map_frames_and_log(Controller* control, Frame* inframe1, Frame* inframe2, Frame* outframe)
{
	//declare variables
	int i, j;
	int type, type_spot;
	
	//SKIP normal mapping for the 2 frames
	
	//determine type and type_num information missing by not doing in-house mapping
	if(control->frame == 1) {
		inframe1->type = malloc(control->num_cg_types * sizeof(int));
		inframe1->type_num = malloc(control->num_cg_types * sizeof(int));
		
		inframe2->type = malloc(control->num_cg_types * sizeof(int));
		inframe2->type_num = malloc(control->num_cg_types * sizeof(int));
	}

	//initailize OR reset values
	for(i = 0; i< control->num_cg_types; i++) {
		inframe1->type[i] = -1;
		inframe1->type_num[i] = 0;
		
		inframe2->type[i] = -1;
		inframe2->type_num[i] = 0;
	}
	inframe1->type_count = 0;
	inframe2->type_count = 0;
	
	//tabulate values
	for(i = 0; i< inframe1->num_atoms; i++) {
		//FRAME 1
		type = inframe1->atoms[i].type;
		type_spot = -1;
		//check if this is new type (if not, increment counter)
		for(j = 0; j <= inframe1->type_count; i++) {
			if( inframe1->type[j] == type) {
				inframe1->type_num[j]++;
				type_spot = j;
			}
		}
		//handle new type
		if(type_spot == -1) {
			inframe1->type_count++;
			inframe1->type[inframe1->type_count] = type;
			inframe1->type_num[inframe1->type_count]++;
		}
		
		//FRAME 2
		type = inframe2->atoms[i].type;
		type_spot = -1;
		//check if this is new type (if not, increment counter)
		for(j = 0; j <= inframe2->type_count; i++) {
			if( inframe2->type[j] == type) {
				inframe2->type_num[j]++;
				type_spot = j;
			}
		}
		//handle new type
		if(type_spot == -1) {
			inframe2->type_count++;
			inframe2->type[inframe2->type_count] = type;
			inframe2->type_num[inframe2->type_count]++;
		}	
	}
	//do combination of files to get final output
	combine_cg_sensitivity_data(control, inframe1, inframe2, outframe);
	//SKIP: free map data since "map frames were not needed"
}

//////////////////////////////////////
///   combine_sensitivity_data	  ///
/////////////////////////////////////

void combine_sensitivity_data(Controller* control, Frame* inframe1, Frame* inframe2, Frame* outframe)
{
	//declare variables
	int i, j, k;
	double scalar = (control->log_value - control->guess) / control->scaleU;
	int paired[inframe2->num_atoms];
	for(i = 0; i < inframe2->num_atoms; i++) {
		paired[i] = 0;
	}
		
	//popluate outframe with necessary general information
	outframe->type_count = 0;
	outframe->xmin = inframe1->xmin;
	outframe->xmax = inframe1->xmax;
	outframe->ymin = inframe1->ymin;
	outframe->ymax = inframe1->ymax;
	outframe->zmin = inframe1->zmin;
	outframe->zmax = inframe1->zmax;
	outframe->timestep = inframe1->timestep;
	outframe->num_observables = inframe1->num_observables;
	
	//allocate type space and initalize
	if(control->frame == 1) {
		outframe->type = malloc(control->num_cg_types * sizeof(int));
		outframe->type_num = malloc(control->num_cg_types * sizeof(int));
	}
		
	outframe->type_count = inframe1->type_count;
	for(i = 0; i< control->num_cg_types; i++) {
		outframe->type[i] = inframe1->type[i];
		outframe->type_num[i] = inframe1->type_num[i];
	}

	//allocate site space
	if(outframe->num_atoms != control->num_cg_sites) {
		outframe->num_atoms = control->num_cg_sites;
		outframe->sites = malloc(outframe->num_atoms * sizeof(SITE));		
	}
	
	if(control->frame == 1) {
		for(i=0; i < outframe->num_atoms; i++) {
			outframe->sites[i].observables = malloc(outframe->num_observables * sizeof(double));
		}
	}
	for(i=0; i < outframe->num_atoms; i++) {
		for(j = 0; j < outframe->num_observables; j++) outframe->sites[i].observables[j] = 0.0;
	}

	//copy frame info on basic info
	for(i = 0; i < outframe->num_atoms; i++) {
		outframe->sites[i].num_in_site = 1;
		outframe->sites[i].id = inframe1->sites[i].id;
		outframe->sites[i].mol = inframe1->sites[i].mol;
		outframe->sites[i].type = inframe1->sites[i].type;
		outframe->sites[i].mass = inframe1->sites[i].mass;
		outframe->sites[i].q = inframe1->sites[i].q;
		outframe->sites[i].x = inframe1->sites[i].x;
		outframe->sites[i].y = inframe1->sites[i].y;
		outframe->sites[i].z = inframe1->sites[i].z;
	
		//look for correct match based on position
		for(j = 0; j < inframe2->num_atoms; j++) {
			//only consider sites that are not yet paired
			if(paired[j] == 1) continue;
			//check each position and skip rest if failed (allow for some rounding error in mapping)
			if( abs(outframe->sites[i].x - inframe2->sites[j].x) > 0.001 ) continue;
			if( abs(outframe->sites[i].y - inframe2->sites[j].y) > 0.001 ) continue;
			if( abs(outframe->sites[i].z - inframe2->sites[j].z) > 0.001 ) continue;
			
			//for matched sites
			paired[j] = 1;
			//combine "observable information"
			for(k = 0; k < outframe->num_observables; k++) {
				outframe->sites[i].observables[k] = inframe2->sites[j].observables[k] - ( scalar * inframe1->sites[i].observables[k];
			}
			break;
		}
	}
}

///////////////////////////////////////
///   combine_cg_sensitivity_data 	///
///////////////////////////////////////

void combine_cg_sensitivity_data(Controller* control, Frame* inframe1, Frame* inframe2, Frame* outframe)
{
	int i, j, k;
	double scalar = control->log_value - control->guess;
	int paired[inframe2->num_atoms];
	
	for(i = 0; i < inframe2->num_atoms; i++) {
		paired[i] = 0;
	}
		
	//popluate outframe with necessary general information
	outframe->type_count = 0;
	outframe->xmin = inframe1->xmin;
	outframe->xmax = inframe1->xmax;
	outframe->ymin = inframe1->ymin;
	outframe->ymax = inframe1->ymax;
	outframe->zmin = inframe1->zmin;
	outframe->zmax = inframe1->zmax;
	outframe->timestep = inframe1->timestep;
	outframe->num_observables = inframe1->num_observables;
	
	//allocate type space and initalize
	if(control->frame == 1) {
		outframe->type = malloc(control->num_cg_types * sizeof(int));
		outframe->type_num = malloc(control->num_cg_types * sizeof(int));
	}
		
	outframe->type_count = inframe1->type_count;
	for(i = 0; i< control->num_cg_types; i++) {
		outframe->type[i] = inframe1->type[i];
		outframe->type_num[i] = inframe1->type_num[i];
	}

	//allocate site space
	if(outframe->num_atoms != control->num_cg_sites) {
		outframe->num_atoms = control->num_cg_sites;
		outframe->sites = malloc(outframe->num_atoms * sizeof(SITE));		
	}
	
	if(control->frame == 1) {
		for(i=0; i < outframe->num_atoms; i++) {
			outframe->sites[i].observables = malloc(outframe->num_observables * sizeof(double));
		}
	}
	for(i=0; i < outframe->num_atoms; i++) {
		for(j = 0; j < outframe->num_observables; j++) outframe->sites[i].observables[j] = 0.0;
	}

	//copy frame info on basic info
	for(i = 0; i < outframe->num_atoms; i++) {
		outframe->sites[i].num_in_site = 1;
		outframe->sites[i].id = inframe1->atoms[i].id;
		outframe->sites[i].mol = inframe1->atoms[i].mol;
		outframe->sites[i].type = inframe1->atoms[i].type;
		outframe->sites[i].mass = inframe1->atoms[i].mass;
		outframe->sites[i].q = inframe1->atoms[i].q;
		outframe->sites[i].x = inframe1->atoms[i].x;
		outframe->sites[i].y = inframe1->atoms[i].y;
		outframe->sites[i].z = inframe1->atoms[i].z;
	
		//look for correct match based on position
		for(j = 0; j < inframe2->num_atoms; j++) {			
			//only consider sites that are not yet paired
			if(paired[j] == 1) continue;
			//check each position and skip rest if failed (allow for some rounding error in mapping)
			if( abs(outframe->sites[i].x - inframe2->atoms[j].x) > 0.001 ) continue;
			if( abs(outframe->sites[i].y - inframe2->atoms[j].y) > 0.001 ) continue;
			if( abs(outframe->sites[i].z - inframe2->atoms[j].z) > 0.001 ) continue;
			
			//for matched sites
			paired[j] = 1;
			//combine "observable information"
			for(k = 0; k < outframe->num_observables; k++) {
				outframe->sites[i].observables[k] = inframe2->atoms[j].observables[k] - scalar * inframe1->atoms[i].observables[k];
			}
			break;
		}
	}
}

/////////////////////////////////////
///   process_ij_charge_frames	 ///
///////////////////////////////////

void process_ij_charge_frames(Controller* control, Frame* inframes, Frame* outframes)
{
	int i, j, k, l;
	int type, type_spot;
	double charge = 0.0;
	int* paired1 = malloc(inframes[1].num_atoms * sizeof(int));
	int* paired2 = malloc(inframes[2].num_atoms * sizeof(int));
	int match1 = -1;
	int match2 = -1;
	for(i = 0; i < inframes[1].num_atoms; i++) {
		paired1[i] = 0;
		paired2[i] = 0;
	}

	//PREPARE INFRAME FORMAT: determine type and type_num information missing by not doing in-house mapping
	if(control->frame == 1) {
		inframes[0].type = malloc(control->num_cg_types * sizeof(int));
		inframes[0].type_num = malloc(control->num_cg_types * sizeof(int));
		inframes[1].type = malloc(control->num_cg_types * sizeof(int));
		inframes[1].type_num = malloc(control->num_cg_types * sizeof(int));
		inframes[2].type = malloc(control->num_cg_types * sizeof(int));
		inframes[2].type_num = malloc(control->num_cg_types * sizeof(int));
		
		outframes[0].type = malloc(control->num_cg_types * sizeof(int));
		outframes[0].type_num = malloc(control->num_cg_types * sizeof(int));
		outframes[1].type = malloc(control->num_cg_types * sizeof(int));
		outframes[1].type_num = malloc(control->num_cg_types * sizeof(int));
	}
	for(k = 0; k< control->num_files; k++) {
		for(i = 0; i< control->num_cg_types; i++) {
			inframes[k].type[i] = -1;
			inframes[k].type_num[i] = 0;
		}
		inframes[k].type_count = 0;
		
		for(i = 0; i< inframes[k].num_atoms; i++) {
			type = inframes[k].atoms[i].type;
			type_spot = -1;
			//check if this is new type (if not, increment counter)
			for(j = 0; j <= inframes[k].type_count; j++) {
				if( inframes[k].type[j] == type) {
					inframes[k].type_num[j]++;
					type_spot = j;
				}
			}
			//handle new type
			if(type_spot == -1) {
				inframes[k].type_count++;
				inframes[k].type[inframes[k].type_count] = type;
				inframes[k].type_num[inframes[k].type_count]++;
			}
		}
	}
			
	//populate general information in outframes
	for(i = 0; i < control->num_outfile; i++) {
		outframes[i].type_count = 0;
		outframes[i].xmin = inframes[0].xmin;
		outframes[i].xmax = inframes[0].xmax;
		outframes[i].ymin = inframes[0].ymin;
		outframes[i].ymax = inframes[0].ymax;
		outframes[i].zmin = inframes[0].zmin;
		outframes[i].zmax = inframes[0].zmax;
		outframes[i].timestep = inframes[0].timestep;
		outframes[i].num_observables = inframes[0].num_observables;	

		//copy types
		for(j = 0; j< control->num_cg_types; j++) {
			outframes[i].type[j] = inframes[0].type[j];
			outframes[i].type_num[j] = inframes[0].type_num[j];
		}
		//allocate site space
		if(outframes[i].num_atoms != inframes[0].num_atoms) {
			outframes[i].num_atoms = control->num_cg_sites;
			outframes[i].sites = malloc(outframes[i].num_atoms * sizeof(SITE));		
		
			for(j=0; j < outframes[i].num_atoms; j++) {
				outframes[i].sites[j].observables = malloc(outframes[i].num_observables * sizeof(double));
			}
		}
		for(j = 0; j < outframes[i].num_atoms; j++) {
			outframes[i].sites[j].num_in_site = 1;
			outframes[i].sites[j].id = inframes[0].atoms[j].id;
			outframes[i].sites[j].mol = inframes[0].atoms[j].mol;
			outframes[i].sites[j].type = inframes[0].atoms[j].type;
			outframes[i].sites[j].mass = inframes[0].atoms[j].mass;
			outframes[i].sites[j].q = inframes[0].atoms[j].q;
			outframes[i].sites[j].x = inframes[0].atoms[j].x;
			outframes[i].sites[j].y = inframes[0].atoms[j].y;
			outframes[i].sites[j].z = inframes[0].atoms[j].z;
		}
	}
	
	//look for correct match based on position
	for(j = 0; j < inframes[0].num_atoms; j++) { //loop over destination sites
		//reset match each time
		match1 = -1;
		match2 = -1;
		
		for(k = 0; k < inframes[1].num_atoms; k++) { //loop over potential inframe#1 site to find match
			//only consider sites that are not yet paired
			if(paired1[k] == 1) continue;
		
			//check each position and skip rest if failed (allow for some rounding error in mapping)
			if( abs(outframes[0].sites[j].x - inframes[1].atoms[k].x) > 0.001 ) continue;
			if( abs(outframes[0].sites[j].y - inframes[1].atoms[k].y) > 0.001 ) continue;
			if( abs(outframes[0].sites[j].z - inframes[1].atoms[k].z) > 0.001 ) continue;
			//for matched sites
			paired1[k] = 1;
			match1 = k;
			break;
		}
		for(k = 0; k < inframes[1].num_atoms; k++) { //loop over potential inframe#1 site to find match	
			//only consider sites that are not yet paired
			if(paired2[k] == 1) continue;
		
			//check each position and skip rest if failed (allow for some rounding error in mapping)
			if( abs(outframes[0].sites[j].x - inframes[2].atoms[k].x) > 0.001 ) continue;
			if( abs(outframes[0].sites[j].y - inframes[2].atoms[k].y) > 0.001 ) continue;
			if( abs(outframes[0].sites[j].z - inframes[2].atoms[k].z) > 0.001 ) continue;
			//for matched sites
			paired2[k] = 1;
			match2 = k;
			break;
		}
		
		//see if a match was found
		if( (match1 == -1) || (match1 == -1) ) {
			printf("no match found for site %d with values 1=%d and 2=%d\n", j, match1, match2);
		} else {
			for(l = 0; l < outframes[0].num_observables; l++) {
				for(k = 0; k < control->num_outfile; k++) {
					outframes[k].sites[j].observables[l] = control->charge[k] * (inframes[0].atoms[j].observables[l] - inframes[1].atoms[match1].observables[l] - inframes[2].atoms[match2].observables[l]);
				}
			}
		}
	}
	free(paired1);
	free(paired2);	
}

/////////////////////////////////////
///   process_ii_charge_frames	 ///
///////////////////////////////////

void process_ii_charge_frames(Controller* control, Frame* inframe, Frame* outframe)
{
	//declare variables
	int i, j;
	int l = 0;
	int type, type_spot;
	double charge = 0.0;
	int* paired = malloc(inframe->num_atoms * sizeof(int));
	int match1 = -1;
	for(i = 0; i < inframe->num_atoms; i++) paired[i] = 0;

	//PREPARE INFRAME FORMAT: determine type and type_num information missing by not doing in-house mapping
	if(control->frame == 1) {
		inframe->type = malloc(control->num_cg_types * sizeof(int));
		inframe->type_num = malloc(control->num_cg_types * sizeof(int));
		outframe->type = malloc(control->num_cg_types * sizeof(int));
		outframe->type_num = malloc(control->num_cg_types * sizeof(int));
	} for(i = 0; i< control->num_cg_types; i++) {
		inframe->type[i] = -1;
		inframe->type_num[i] = 0;
	}
	inframe->type_count = 0;
	
	for(i = 0; i< inframe->num_atoms; i++) {
		//FRAME 1
		type = inframe->atoms[i].type;
		type_spot = -1;
		//check if this is new type (if not, increment counter)
		for(j = 0; j <= inframe->type_count; j++) {
			if( inframe->type[j] == type) {
				inframe->type_num[j]++;
				type_spot = j;
			}
		}
		//handle new type
		if(type_spot == -1) {
			inframe->type_count++;
			inframe->type[inframe->type_count] = type;
			inframe->type_num[inframe->type_count]++;
		}
	}	

	//populate general information in outframes
	outframe->type_count = 0;
	outframe->xmin = inframe->xmin;
	outframe->xmax = inframe->xmax;
	outframe->ymin = inframe->ymin;
	outframe->ymax = inframe->ymax;
	outframe->zmin = inframe->zmin;
	outframe->zmax = inframe->zmax;
	outframe->timestep = inframe->timestep;
	outframe->num_observables = inframe->num_observables;	

	//copy types
	for(j = 0; j< control->num_cg_types; j++) {
		outframe->type[j] = inframe->type[j];
		outframe->type_num[j] = inframe->type_num[j];
	}
	//allocate site space
	if(outframe->num_atoms != control->num_cg_sites) {
		outframe->num_atoms = control->num_cg_sites;
		outframe->sites = malloc(outframe->num_atoms * sizeof(SITE));		
		for(i=0; i < outframe->num_atoms; i++) {
			outframe->sites[i].observables = malloc(outframe->num_observables * sizeof(double));
		}
	}
	for(j = 0; j < outframe->num_atoms; j++) {
		outframe->sites[j].num_in_site = 1;
		outframe->sites[j].id = inframe->atoms[j].id;
		outframe->sites[j].mol = inframe->atoms[j].mol;
		outframe->sites[j].type = inframe->atoms[j].type;
		outframe->sites[j].mass = inframe->atoms[j].mass;
		outframe->sites[j].q = inframe->atoms[j].q;
		outframe->sites[j].x = inframe->atoms[j].x;
		outframe->sites[j].y = inframe->atoms[j].y;
		outframe->sites[j].z = inframe->atoms[j].z;
		for(l = 0; l < outframe->num_observables; l++) {
			outframe->sites[j].observables[l] = 2.0 * control->charge[0] * inframe->atoms[j].observables[l];
		}
		}
	free(paired);
}

//////////////////////////////////////////////
///   free_sensitivity_intermediates	  ///
/////////////////////////////////////////////

void free_sensitivity_intermediates(Controller* control, Frame* map1, Frame* map2)
{
	int i;
	
	//free all data for atoms
	for(i = 0; i < map1->num_atoms; i++) {
		free(map1->sites[i].observables);
		free(map2->sites[i].observables);
		
		free(map1->sites[i].coord);
		free(map2->sites[i].coord);
		
		free(map1->sites[i].matches);
		free(map2->sites[i].matches);
	}
		
	free(map1->sites);
	free(map2->sites);
	
	free(map1->type);
	free(map2->type);
	
	free(map1->type_num);
	free(map2->type_num);
}

void free_charge_intermediates(Controller* control, Frame* mapframes)
{
	int i, j;
	
	for(i = 0; i < control->num_files; i++) {
		//free all data for atoms
		for(j = 0; j < mapframes[0].num_atoms; j++) {
			free(mapframes[i].sites[j].observables);	
			free(mapframes[i].sites[j].coord);
			free(mapframes[i].sites[j].matches);
		}
		
		free(mapframes[i].sites);
		
		free(mapframes[i].type);
		free(mapframes[i].type_num);
	}
}