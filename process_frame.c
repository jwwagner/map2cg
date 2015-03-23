//process_frame.c

#include "headers.h"

//original functions
void process_frame(Controller*, Frame*, Frame*);
void process_minimal_frames(Controller*, Frame*, Frame*);
void process_frames_and_log(Controller*, Frame*, Frame*, Frame*);
void process_no_map_frames_and_log(Controller*, Frame*, Frame*, Frame*);
void map_all_atoms(Controller*, Frame*, Frame*); 
void map_some_atoms(Controller*, Frame*, Frame*);
void map_beginning_atoms(Controller*, Frame*, Frame*);
//sensitivity functions
void combine_sensitivity_data(Controller*, Frame*, Frame*, Frame*);
void combine_cg_sensitivity_data(Controller*, Frame*, Frame*, Frame*);
void process_frame_order(Controller*, Frame*, Frame*);
void process_ij_charge_frames(Controller*, Frame*, Frame*);
void process_ii_charge_frames(Controller*, Frame*, Frame*);
void process_file_merge(Controller*, Frame*, double*);
//slave functions (newer modularization)
void determine_type(Controller*, Frame*, Frame*, int*, int*);
int compatible_type_test(Controller*, Frame*, Frame*, int, int);
void geometry_mapping(Controller*, Frame*, Frame*);
//free functions
void free_sensitivity_intermediates(Controller*, Frame*, Frame*);
void free_charge_intermediates(Controller*, Frame*);

//////////////////////////
///   process_frame	  ///
/////////////////////////

void process_frame(Controller* control, Frame* inframe, Frame* outframe)
{
	int i, j;
	outframe->type_count = 0;
	
	if( (control->frame == 1) || (control->sensitivity_flag == 1) )
		{
		outframe->type = malloc(control->num_cg_types * sizeof(int));
		outframe->type_num = malloc(control->num_cg_types * sizeof(int));
		}
	
	//intitalize type array
	for(i = 0; i< control->num_cg_types; i++)
		{
		outframe->type[i] = -1;
		outframe->type_num[i] = 0;
		}

	//copy basic information from inframe to outframe
	outframe->xmin = inframe->xmin;
	outframe->xmax = inframe->xmax;
	outframe->ymin = inframe->ymin;
	outframe->ymax = inframe->ymax;
	outframe->zmin = inframe->zmin;
	outframe->zmax = inframe->zmax;
	
	outframe->timestep = inframe->timestep;
	outframe->num_observables = inframe->num_observables;
	
	//determine if we need to allocate sites
	if(control->num_cg_sites != outframe->num_atoms) //assumes unitialized (0) although flexibility could be added to change number of sites (e.g. Grand Canonical Ensemble)
		{
		printf("allocate sites for outframe to size %d\n", control->num_cg_sites);
		outframe->num_atoms = control->num_cg_sites;
		outframe->sites = malloc(outframe->num_atoms * sizeof(SITE));
		
		//also need to allocate observables
		for(i=0; i < outframe->num_atoms; i++)
			{
			outframe->sites[i].observables = malloc(outframe->num_observables * sizeof(double));
			outframe->sites[i].coord = malloc(control->max_to_map *sizeof(COORD));
			outframe->sites[i].matches = malloc(control->num_cg_types * sizeof(int));
			}
		printf("observables and coord space allocated\n");
		}

	//reset frame info
	for(i = 0; i < outframe->num_atoms; i++)
		{
		outframe->sites[i].q = 0.0;
		outframe->sites[i].num_in_site = 0;
		outframe->sites[i].type = -1;
		
		for(j = 0; j < outframe->num_observables; j++)
			{
			outframe->sites[i].observables[j] = 0.0;
			}
		
		for(j = 0; j < control->max_to_map; j++)
			{
			outframe->sites[i].coord[j].x = 0.0;
			outframe->sites[i].coord[j].y = 0.0;
			outframe->sites[i].coord[j].z = 0.0;
			}
			
		for(j = 0; j < control->num_cg_types; j++)
			{
			outframe->sites[i].matches[j] = 1;
			}
		}
	
	//accumulate information from fine_grained sites as we go (assuming molecule to molecule mapping)
	if(control->map_style_flag == 0) 
		{
		map_all_atoms(control, inframe, outframe);
		}
	else if(control->map_style_flag == 1)
		{
		map_some_atoms(control, inframe, outframe);
		}
	else if(control->map_style_flag == 2)
		{
		map_beginning_atoms(control, inframe, outframe);
		}
}

//////////////////////////////////
///   process_minimal_frame	  ///
////////////////////////////////

void process_minimal_frame(Controller* control, Frame* inframe, Frame* outframe)
{
	//declare variables
	int i, j;
	int l = 0;
	int type, type_spot;
	
	//PREPARE INFRAME FORMAT
	//determine type and type_num information missing by not doing in-house mapping
	//allocate space
	if(control->frame == 1)
		{
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
	if(outframe->num_atoms != control->num_cg_sites)
		{
		outframe->num_atoms = control->num_cg_sites;
		outframe->sites = malloc(outframe->num_atoms * sizeof(SITE));		
		for(i=0; i < outframe->num_atoms; i++)
			{
			outframe->sites[i].observables = malloc(outframe->num_observables * sizeof(double));
			}
		}
	for(j = 0; j < outframe->num_atoms; j++)
		{
		outframe->sites[j].num_in_site = 1;
		outframe->sites[j].id = j;
		outframe->sites[j].mol = j;
		outframe->sites[j].type = 1;
		outframe->sites[j].mass = control->scaleU2;
		outframe->sites[j].q = control->scaleU1;
		outframe->sites[j].x = inframe->atoms[j].x;
		outframe->sites[j].y = inframe->atoms[j].y;
		outframe->sites[j].z = inframe->atoms[j].z;
		for(l = 0; l < outframe->num_observables; l++)
			{
			outframe->sites[j].observables[l] = inframe->atoms[j].observables[l];
			}
		}
}

/////////////////////////////////
///   process_frames_and_log  ///
/////////////////////////////////

void process_frames_and_log(Controller* control, Frame* inframe1, Frame* inframe2, Frame* outframe)
{
	Frame map1;
	Frame map2;
	map1.num_atoms = 0;
	map2.num_atoms = 0;
	
	//do normal mapping for the 2 frames
	process_frame(control, inframe1, &map1);
	process_frame(control, inframe2, &map2);
	
	//do combination of files to get final output
	combine_sensitivity_data(control, &map1, &map2, outframe);

	//free map data
	free_sensitivity_intermediates(control, &map1, &map2);
}

void process_no_map_frames_and_log(Controller* control, Frame* inframe1, Frame* inframe2, Frame* outframe)
{
	//declare variables
	int i, j;
	int type, type_spot;
	
	//SKIP normal mapping for the 2 frames
	
	//determine type and type_num information missing by not doing in-house mapping
	if(control->frame == 1)
		{
		inframe1->type = malloc(control->num_cg_types * sizeof(int));
		inframe1->type_num = malloc(control->num_cg_types * sizeof(int));
		
		inframe2->type = malloc(control->num_cg_types * sizeof(int));
		inframe2->type_num = malloc(control->num_cg_types * sizeof(int));
		}

	//initailize OR reset values
	for(i = 0; i< control->num_cg_types; i++)
		{
		inframe1->type[i] = -1;
		inframe1->type_num[i] = 0;
		
		inframe2->type[i] = -1;
		inframe2->type_num[i] = 0;
		}
	inframe1->type_count = 0;
	inframe2->type_count = 0;
	
	//tabulate values
	for(i = 0; i< inframe1->num_atoms; i++)
		{
		//FRAME 1
		type = inframe1->atoms[i].type;
		
		type_spot = -1;
		//check if this is new type (if not, increment counter)
		for(j = 0; j <= inframe1->type_count; j++)
			{
			if( inframe1->type[j] == type)
				{
				inframe1->type_num[j]++;
				type_spot = j;
				}
			}
		
		//handle new type
		if(type_spot == -1)
			{
			inframe1->type_count++;
			inframe1->type[inframe1->type_count] = type;
			inframe1->type_num[inframe1->type_count]++;
			}
		
		//FRAME 2
		type = inframe2->atoms[i].type;
		type_spot = -1;
		//check if this is new type (if not, increment counter)
		for(j = 0; j <= inframe2->type_count; j++)
			{
			if( inframe2->type[j] == type)
				{
				inframe2->type_num[j]++;
				type_spot = j;
				}
			}
		//handle new type
		if(type_spot == -1)
			{
			inframe2->type_count++;
			inframe2->type[inframe2->type_count] = type;
			inframe2->type_num[inframe2->type_count]++;
			}	
		}
	//do combination of files to get final output
	combine_cg_sensitivity_data(control, inframe1, inframe2, outframe);
	//SKIP: free map data since "map frames were not needed"
}

//////////////////////////
///   map_all_atoms	  ///
/////////////////////////

void map_all_atoms(Controller* control, Frame* inframe, Frame* outframe)
{
	int i, j, k, l;
	int mol_count = 0;
	int mol_val = 0;
	int real_mol_val = 0;
	int site_count;
	int num_key = inframe->num_mol * control->num_cg_types;
	int key[num_key];
	int assigned;
	int testid = 0;
	
	//intialize key
	//printf("process frame with inframe->num_mol %d \n", inframe->num_mol);
	for(i = 0; i < num_key; i++) {
		key[i] = -1;
		
		for(j=0; j < control->num_cg_types; j++) {
			outframe->sites[i].matches[j] = 1;
		}	
	}

	if(control->observable_map_flag == 0) { //sum observables
		//process and sort all FG atoms
		for(i = 0; i < inframe->num_atoms; i++) {
			//check to see if mol key is set
			mol_val = (inframe->atoms[i].mol - 1) * control->num_cg_types;
			//printf("atom number is %d in mol %d with mol_val %d and key %d\n", i, inframe->atoms[i].mol, mol_val, key[mol_val]);
			if(key[mol_val] == -1) {
				//printf("assign new type first\n");
				key[mol_val] = mol_count;
				outframe->sites[key[mol_val]].id = key[mol_val] + 1;
				outframe->sites[key[mol_val]].mol = real_mol_val + 1;
				outframe->sites[key[mol_val]].num_in_site = 1;
				mol_count++;
				real_mol_val++;
			} else if( compatible_type_test(control, inframe, outframe, inframe->atoms[i].type, key[mol_val]) == -1 ) {
				//this new atom is NOT compatible with original CG type
				testid = mol_val;
				
				//printf("comp: check assignment\n");
				//test if we have already assigned another type (that this would be compatible with)
				assigned = -1;
				for(j = 1; j < control->num_cg_types; j++) {
					//check if next type is allocated
					if( key[mol_val + j] == -1 ) {
						testid = mol_val + j - 1;
						break;
					} else {
						//key is assigned
						testid = mol_val + j; 
						//test if key is compatible
						//printf("test next key\n");
						if(compatible_type_test(control, inframe, outframe, inframe->atoms[i].type, key[testid]) == 1) {		
							assigned = 1;
							mol_val = testid;
							break;
						}
						//otherwise, keep looking with next key
					}
				}
				if(assigned == -1) {
					//So, we still need to create a new type
					//printf("create new mol_val from %d and testid %d\n", mol_val, testid);
					
					mol_val = testid + 1;
					
					key[mol_val] = mol_count;
					outframe->sites[key[mol_val]].id = key[mol_val] + 1;
					outframe->sites[key[mol_val]].mol = outframe->sites[key[(inframe->atoms[i].mol - 1) * control->num_cg_types] ].mol;
					outframe->sites[key[mol_val]].num_in_site = 1;
					mol_count++;
				}
			} else {
				//printf("type is compatible\n");
			}

			//transfer information to CG site
			site_count = outframe->sites[ key[mol_val] ].num_in_site;		 
			outframe->sites[key[mol_val]].coord[site_count].x = inframe->atoms[i].x;
			outframe->sites[key[mol_val]].coord[site_count].y = inframe->atoms[i].y;
			outframe->sites[key[mol_val]].coord[site_count].z = inframe->atoms[i].z;
			outframe->sites[key[mol_val]].coord[site_count].mass = inframe->atoms[i].mass;
			outframe->sites[key[mol_val]].coord[site_count].type = inframe->atoms[i].type;
			outframe->sites[key[mol_val]].q += inframe->atoms[i].q;
			outframe->sites[ key[mol_val] ].num_in_site++;

			determine_type(control, inframe, outframe, &i, &key[mol_val]);
	
			for(j = 0; j < outframe->num_observables; j++) {
				outframe->sites[key[mol_val]].observables[j] += inframe->atoms[i].observables[j];
			}
		}			
		geometry_mapping(control, inframe, outframe);
	}
	outframe->num_mol = mol_count;	
}

/////////////////////////////
///   map_some_atoms	  ///
/////////////////////////////

void map_some_atoms(Controller* control, Frame* inframe, Frame* outframe)
{
	int i, j, k;
	int mol_count = 0;
	int mol_val;
	int site_count;
	int go_flag;
	int key[inframe->num_mol];
	
	//intialize key
	for(i = 0; i < inframe->num_mol; i++)
		{
		key[i] = -1;
		}
		
	if(control->observable_map_flag == 0) //sum observables
		{
		//process and sort all FG atoms
		for(i = 0; i < inframe->num_atoms; i++)
			{			
			//see if molecule is in include list
			go_flag = 0;
			for(j = 0; j < control->num_map; j++)
				{
				if(control->map[j] == inframe->atoms[i].type)
					{
					go_flag = 1;
					break;
					}
				}
			
			//everything from this point forward is the same as in map_all_atoms
			if(go_flag == 0) continue;
			
			//check to see if mol key is set
			mol_val = inframe->atoms[i].mol - 1;
			if(key[mol_val] == -1)
				{
				key[mol_val] = mol_count;
				outframe->sites[key[mol_val]].id = key[mol_val] + 1;
				outframe->sites[key[mol_val]].mol = key[mol_val] + 1;
				mol_count++;
				}	
			determine_type(control, inframe, outframe, &i, &key[mol_val]);
				
			//transfer information to CG site
			site_count = outframe->sites[ key[mol_val] ].num_in_site;
			
			outframe->sites[key[mol_val]].coord[site_count].x = inframe->atoms[i].x;
			outframe->sites[key[mol_val]].coord[site_count].y = inframe->atoms[i].y;
			outframe->sites[key[mol_val]].coord[site_count].z = inframe->atoms[i].z;
			outframe->sites[key[mol_val]].coord[site_count].mass = inframe->atoms[i].mass;
			outframe->sites[key[mol_val]].coord[site_count].type = inframe->atoms[i].type;
			outframe->sites[key[mol_val]].q += inframe->atoms[i].q;
			
			for(j = 0; j < outframe->num_observables; j++)
				{
				outframe->sites[key[mol_val]].observables[j] += inframe->atoms[i].observables[j];
				}
			outframe->sites[ key[mol_val] ].num_in_site++;
			}
				
		geometry_mapping(control, inframe, outframe);
		}
	outframe->num_mol = mol_count;		
}

/////////////////////////////
///   map_beginning_atoms ///
/////////////////////////////

void map_beginning_atoms(Controller* control, Frame* inframe, Frame* outframe) 
{
	int i, j, k, l;
	int mol_count = 0;
	int mol_val;
	int site_count;
	int key[inframe->num_mol];
	
	//intialize key
	for(i = 0; i < inframe->num_mol; i++)
		{
		key[i] = -1;
		}

	if(control->observable_map_flag == 0) //sum observables
		{
		//process and sort all FG atoms
		for(i = 0; i < inframe->num_atoms; i++)
			{
			//check if atom is within truncation range
			if(inframe->atoms[i].id > control->num_truncate) continue;
			
			//check to see if mol key is set
			mol_val = inframe->atoms[i].mol - 1;
			if(key[mol_val] == -1)
				{
				key[mol_val] = mol_count;
				outframe->sites[key[mol_val]].id = key[mol_val] + 1;
				outframe->sites[key[mol_val]].mol = key[mol_val] + 1;
				mol_count++;
				}												
		
			//transfer information to CG site
			site_count = outframe->sites[ key[mol_val] ].num_in_site;		 
			outframe->sites[key[mol_val]].coord[site_count].x = inframe->atoms[i].x;
			outframe->sites[key[mol_val]].coord[site_count].y = inframe->atoms[i].y;
			outframe->sites[key[mol_val]].coord[site_count].z = inframe->atoms[i].z;
			outframe->sites[key[mol_val]].coord[site_count].mass = inframe->atoms[i].mass;
			outframe->sites[key[mol_val]].coord[site_count].type = inframe->atoms[i].type;
			outframe->sites[key[mol_val]].q += inframe->atoms[i].q;
			outframe->sites[ key[mol_val] ].num_in_site++;
			
			determine_type(control, inframe, outframe, &i, &key[mol_val]);
	
			for(j = 0; j < outframe->num_observables; j++)
				{
				outframe->sites[key[mol_val]].observables[j] += inframe->atoms[i].observables[j];
				}
			}
				
		geometry_mapping(control, inframe, outframe);
		}
	outframe->num_mol = mol_count;	
}

//////////////////////////////////////
///   combine_sensitivity_data	  ///
/////////////////////////////////////

void combine_sensitivity_data(Controller* control, Frame* inframe1, Frame* inframe2, Frame* outframe)
{
	//declare variables
	int i, j, k;
	double scalar = control->log_value / control->scaleU1 - control->guess / control->scaleU2;
	int paired[inframe2->num_atoms];
	for(i = 0; i < inframe2->num_atoms; i++) 
		{
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
	if(control->frame == 1)
		{
		outframe->type = malloc(control->num_cg_types * sizeof(int));
		outframe->type_num = malloc(control->num_cg_types * sizeof(int));
		}
		
	outframe->type_count = inframe1->type_count;
	for(i = 0; i < control->num_cg_types; i++)
		{
		outframe->type[i] = inframe1->type[i];
		outframe->type_num[i] = inframe1->type_num[i];
		}

	//allocate site space
	if(outframe->num_atoms != control->num_cg_sites)
		{
		outframe->num_atoms = control->num_cg_sites;
		outframe->sites = malloc(outframe->num_atoms * sizeof(SITE));		
		}
	
	if(control->frame == 1)
		{
		for(i=0; i < outframe->num_atoms; i++)
			{
			outframe->sites[i].observables = malloc(outframe->num_observables * sizeof(double));
			}
		}
	for(i=0; i < outframe->num_atoms; i++)
		{
		for(j = 0; j < outframe->num_observables; j++) outframe->sites[i].observables[j] = 0.0;
		}

	//copy frame info on basic info
	for(i = 0; i < outframe->num_atoms; i++)
		{
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
		for(j = 0; j < inframe2->num_atoms; j++)
			{
			//only consider sites that are not yet paired
			if(paired[j] == 1) continue;
			//check each position and skip rest if failed (allow for some rounding error in mapping)
			if( abs(outframe->sites[i].x - inframe2->sites[j].x) > 0.001 ) continue;
			if( abs(outframe->sites[i].y - inframe2->sites[j].y) > 0.001 ) continue;
			if( abs(outframe->sites[i].z - inframe2->sites[j].z) > 0.001 ) continue;
			
			//for matched sites
			paired[j] = 1;
			//combine "observable information"
			for(k = 0; k < outframe->num_observables; k++)
				{
				outframe->sites[i].observables[k] = inframe2->sites[j].observables[k] + ((double)control->sign_flag) * scalar * inframe1->sites[i].observables[k];
				outframe->sites[i].observables[k] /= control->scaleF;
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
	//declare variables
	int i, j, k;
	double scalar = control->log_value - control->guess;
	int paired[inframe2->num_atoms];
	
	for(i = 0; i < inframe2->num_atoms; i++) 
		{
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
	if(control->frame == 1)
		{
		outframe->type = malloc(control->num_cg_types * sizeof(int));
		outframe->type_num = malloc(control->num_cg_types * sizeof(int));
		}
		
	outframe->type_count = inframe1->type_count;
	for(i = 0; i< control->num_cg_types; i++)
		{
		outframe->type[i] = inframe1->type[i];
		outframe->type_num[i] = inframe1->type_num[i];
		}

	//allocate site space
	if(outframe->num_atoms != control->num_cg_sites)
		{
		outframe->num_atoms = control->num_cg_sites;
		outframe->sites = malloc(outframe->num_atoms * sizeof(SITE));		
		}
	
	if(control->frame == 1)
		{
		for(i=0; i < outframe->num_atoms; i++)
			{
			outframe->sites[i].observables = malloc(outframe->num_observables * sizeof(double));
			}
		}
	for(i=0; i < outframe->num_atoms; i++)
		{
		for(j = 0; j < outframe->num_observables; j++) outframe->sites[i].observables[j] = 0.0;
		}

	//copy frame info on basic info
	for(i = 0; i < outframe->num_atoms; i++)
		{
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
		for(j = 0; j < inframe2->num_atoms; j++)
			{			
			//only consider sites that are not yet paired
			if(paired[j] == 1) continue;
			//check each position and skip rest if failed (allow for some rounding error in mapping)
			if( abs(outframe->sites[i].x - inframe2->atoms[j].x) > 0.001 ) continue;
			if( abs(outframe->sites[i].y - inframe2->atoms[j].y) > 0.001 ) continue;
			if( abs(outframe->sites[i].z - inframe2->atoms[j].z) > 0.001 ) continue;
			
			//for matched sites
			paired[j] = 1;
			//combine "observable information"
			for(k = 0; k < outframe->num_observables; k++)
				{
				outframe->sites[i].observables[k] = inframe2->atoms[j].observables[k] - scalar * inframe1->atoms[i].observables[k];
				}
			break;
			}
		}
}

/////////////////////////////////
///   process_frame_order	 ///
////////////////////////////////

void process_frame_order(Controller* control, Frame* inframe, Frame* outframe)
{
	int i, j, l, id, spot, prev;
	int type, type_spot;
		
	if(control->frame == 1)
		{
		outframe->type = malloc(control->num_cg_types * sizeof(int));
		outframe->type_num = malloc(control->num_cg_types * sizeof(int));
		inframe->type = malloc(control->num_cg_types * sizeof(int));
		inframe->type_num = malloc(control->num_cg_types * sizeof(int));
		control->order = malloc(control->num_cg_types * sizeof(int));
		
		inframe->type_list = malloc(control->num_cg_types * sizeof(int*));
		for(i = 0; i < control->num_cg_types; i++)
			{
			inframe->type_list[i] = malloc( inframe->num_atoms * sizeof(int) );
			}
		}
	
	//intitalize type array
	outframe->type_count = 0;
	inframe->type_count = 0;
	for(i = 0; i < control->num_cg_types; i++)
		{
		outframe->type[i] = -1;
		outframe->type_num[i] = 0;
		inframe->type[i] = -1;
		inframe->type_num[i] = 0;
		}
		
	//determine if we need to allocate sites
	if(control->num_cg_sites != outframe->num_atoms) //assumes unitialized (0) although flexibility could be added to change number of sites (e.g. Grand Canonical Ensemble)
		{
		printf("allocate sites for outframe to size %d\n", control->num_cg_sites);
		outframe->num_atoms = control->num_cg_sites;
		outframe->sites = malloc(outframe->num_atoms * sizeof(SITE));
		
		//also need to allocate observables
		for(i=0; i < outframe->num_atoms; i++)
			{
			outframe->sites[i].observables = malloc(outframe->num_observables * sizeof(double));
			outframe->sites[i].coord = malloc(control->max_to_map *sizeof(COORD));
			outframe->sites[i].matches = malloc(control->num_cg_types * sizeof(int));
			}
		printf("observables and coord space allocated\n");
		}

	//reset frame info
	for(i = 0; i < outframe->num_atoms; i++)
		{
		outframe->sites[i].q = 0.0;
		outframe->sites[i].num_in_site = 0;
		outframe->sites[i].type = -1;	
		
		for(j = 0; j < outframe->num_observables; j++)
			{
			outframe->sites[i].observables[j] = 0.0;
			}
		
		for(j = 0; j < control->max_to_map; j++)
			{
			outframe->sites[i].coord[j].x = 0.0;
			outframe->sites[i].coord[j].y = 0.0;
			outframe->sites[i].coord[j].z = 0.0;
			}
			
		for(j = 0; j < control->num_cg_types; j++)
			{
			outframe->sites[i].matches[j] = 1;
			}
		}

	//CREATE OUTPUT
	//copy basic information from inframe to outframe
	outframe->xmin = inframe->xmin;
	outframe->xmax = inframe->xmax;
	outframe->ymin = inframe->ymin;
	outframe->ymax = inframe->ymax;
	outframe->zmin = inframe->zmin;
	outframe->zmax = inframe->zmax;
	outframe->timestep = inframe->timestep;
	outframe->num_observables = inframe->num_observables;
	
	//SORT ATOM_ID BY TYPE INTO TYPE_LIST BINS
	for(i = 0; i < inframe->num_atoms; i++)
		{
		//look-up type
		type = inframe->atoms[i].type;
		type_spot = -1;
		for(j = 0; j < inframe->type_count; j++)
			{
			if( inframe->type[j] == type)
				{
				inframe->type_num[j]++;
				type_spot = j;
				}
			}
		//handle new type
		if(type_spot == -1)
			{
			type_spot = inframe->type_count;
			inframe->type[inframe->type_count] = type;
			inframe->type_num[inframe->type_count] = 1;
			inframe->type_count++;
			}
		//catalog type into type_list array
		inframe->type_list[ type_spot ][ inframe->type_num[type_spot] - 1] = i;
		}
	
	//GET ORDER FOR OUTPUT
	prev = -1;
	//int type_val = control->num_cg_types;
	for(i = 0; i < control->num_cg_types; i++)
		{
		int type_val = control->num_cg_types;
		for(j = 0; j < control->num_cg_types; j++)
			{
			if (inframe->type[j] <= type_val) //get the lowest type value as "low water mark"
				{
				if(inframe->type[j] > prev) //exclude the same type from matching repeatedly
					{
					type = j;
					type_val = inframe->type[j];
					}
				}
			}
		control->order[i] = type_val;
		prev = inframe->type[type];
		}
	
	spot = 0;
	for(i = 0; i < control->num_cg_types; i++)
		{
		//copy type info in order it is placed in frame
		prev = -1;
		for(j = 0; j < control->num_cg_types; j++)
			{
			//get type index
			if( inframe->type[j] == control->order[i] )
				{
				prev = j;
				}
			}
		if(prev == -1) 
			{
			printf("ERROR: NO TYPE MATCH FOUND FOR ORDER #%d TYPE %d (in process frame order)\n", i, inframe->type[i]);
			}
		outframe->type[i] = inframe->type[ prev ];
		outframe->type_num[i] = inframe->type_num[ prev ];
		
		//copy frame info for all entries of this type
		for(j = 0; j < inframe->type_num[ prev ]; j++)
			{
			id = inframe->type_list[ prev ][j];
			outframe->sites[spot].num_in_site = 1;
			outframe->sites[spot].id = spot+1; //inframe->atoms[id].id;
			outframe->sites[spot].mol = spot+1; //inframe->atoms[id].mol;
			outframe->sites[spot].type = inframe->atoms[id].type;
			outframe->sites[spot].q = inframe->atoms[id].q;
			outframe->sites[spot].mass = inframe->atoms[id].mass;
			outframe->sites[spot].x = inframe->atoms[id].x;
			outframe->sites[spot].y = inframe->atoms[id].y;
			outframe->sites[spot].z = inframe->atoms[id].z;
			for(l = 0; l < outframe->num_observables; l++)
				{
				outframe->sites[spot].observables[l] = inframe->atoms[id].observables[l];
				}
			spot++;
			}
		}
}

/////////////////////////////////////
///   process_ij_charge_frames	 ///
///////////////////////////////////

void process_ij_charge_frames(Controller* control, Frame* inframes, Frame* outframes)
{
	//declare variables
	int i, j, k, l;
	int type, type_spot;
	double charge = 0.0;
	int* paired1 = malloc(inframes[1].num_atoms * sizeof(int));
	int* paired2 = malloc(inframes[2].num_atoms * sizeof(int));
	int match1 = -1;
	int match2 = -1;
	for(i = 0; i < inframes[1].num_atoms; i++)
		{
		paired1[i] = 0;
		paired2[i] = 0;
		}

	//PREPARE INFRAME FORMAT
	//determine type and type_num information missing by not doing in-house mapping
	//allocate space
	if(control->frame == 1)
		{
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
	for(k = 0; k< control->num_files; k++)
		{
		for(i = 0; i< control->num_cg_types; i++)
			{
			inframes[k].type[i] = -1;
			inframes[k].type_num[i] = 0;
			}
		inframes[k].type_count = 0;
		
		for(i = 0; i< inframes[k].num_atoms; i++)
			{
			type = inframes[k].atoms[i].type;
			type_spot = -1;
			//check if this is new type (if not, increment counter)
			for(j = 0; j <= inframes[k].type_count; j++)
				{
				if( inframes[k].type[j] == type)
					{
					inframes[k].type_num[j]++;
					type_spot = j;
					}
				}
			//handle new type
			if(type_spot == -1)
				{
				inframes[k].type_count++;
				inframes[k].type[inframes[k].type_count] = type;
				inframes[k].type_num[inframes[k].type_count]++;
				}
			}
		}
			
	//populate general information in outframes
	for(i = 0; i < control->num_outfile; i++)
		{
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
		for(j = 0; j< control->num_cg_types; j++)
			{
			outframes[i].type[j] = inframes[0].type[j];
			outframes[i].type_num[j] = inframes[0].type_num[j];
			}
		//allocate site space
		if(outframes[i].num_atoms != inframes[0].num_atoms)
			{
			outframes[i].num_atoms = control->num_cg_sites;
			outframes[i].sites = malloc(outframes[i].num_atoms * sizeof(SITE));		
		
			for(j=0; j < outframes[i].num_atoms; j++)
				{
				outframes[i].sites[j].observables = malloc(outframes[i].num_observables * sizeof(double));
				}
			}
		for(j = 0; j < outframes[i].num_atoms; j++)
			{
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
	for(j = 0; j < inframes[0].num_atoms; j++) //loop over destination sites
		{
		//reset match each time
		match1 = -1;
		match2 = -1;
		
		for(k = 0; k < inframes[1].num_atoms; k++) //loop over potential inframe#1 site to find match
			{	
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
		for(k = 0; k < inframes[1].num_atoms; k++) //loop over potential inframe#1 site to find match
			{	
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
		if( (match1 == -1) || (match1 == -1) )
			{
			printf("no match found for site %d with values 1=%d and 2=%d\n", j, match1, match2);
			}
		else
			{
			for(l = 0; l < outframes[0].num_observables; l++)
				{
				for(k = 0; k < control->num_outfile; k++)
					{
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

	//PREPARE INFRAME FORMAT
	//determine type and type_num information missing by not doing in-house mapping
	//allocate space
	if(control->frame == 1)
		{
		inframe->type = malloc(control->num_cg_types * sizeof(int));
		inframe->type_num = malloc(control->num_cg_types * sizeof(int));
		outframe->type = malloc(control->num_cg_types * sizeof(int));
		outframe->type_num = malloc(control->num_cg_types * sizeof(int));
		}
	for(i = 0; i< control->num_cg_types; i++)
		{
		inframe->type[i] = -1;
		inframe->type_num[i] = 0;
		}
	inframe->type_count = 0;
	
	for(i = 0; i< inframe->num_atoms; i++)
		{
		//FRAME 1
		type = inframe->atoms[i].type;
		type_spot = -1;
		//check if this is new type (if not, increment counter)
		for(j = 0; j <= inframe->type_count; j++)
			{
			if( inframe->type[j] == type)
				{
				inframe->type_num[j]++;
				type_spot = j;
				}
			}
		//handle new type
		if(type_spot == -1)
			{
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
	for(j = 0; j< control->num_cg_types; j++)
		{
		outframe->type[j] = inframe->type[j];
		outframe->type_num[j] = inframe->type_num[j];
		}
	//allocate site space
	if(outframe->num_atoms != control->num_cg_sites)
		{
		outframe->num_atoms = control->num_cg_sites;
		outframe->sites = malloc(outframe->num_atoms * sizeof(SITE));		
		for(i=0; i < outframe->num_atoms; i++)
			{
			outframe->sites[i].observables = malloc(outframe->num_observables * sizeof(double));
			}
		}
	for(j = 0; j < outframe->num_atoms; j++)
		{
		outframe->sites[j].num_in_site = 1;
		outframe->sites[j].id = inframe->atoms[j].id;
		outframe->sites[j].mol = inframe->atoms[j].mol;
		outframe->sites[j].type = inframe->atoms[j].type;
		outframe->sites[j].mass = inframe->atoms[j].mass;
		outframe->sites[j].q = inframe->atoms[j].q;
		outframe->sites[j].x = inframe->atoms[j].x;
		outframe->sites[j].y = inframe->atoms[j].y;
		outframe->sites[j].z = inframe->atoms[j].z;
		for(l = 0; l < outframe->num_observables; l++)
			{
			outframe->sites[j].observables[l] = 2.0 * control->charge[0] * inframe->atoms[j].observables[l];
			}
		}
	free(paired);
}

/////////////////////////////
///   process_file_merge  ///
/////////////////////////////

void process_file_merge(Controller* control, Frame* frame, double* observable)
{
	int i;
	
	for(i=0; i < control->num_cg_sites; i++) 
		{
		frame->sites[i].observables[0] = observable[i];
		}
}

/////////////////////////////
///   determine_type	 ///
///////////////////////////

void determine_type(Controller* control, Frame* inframe, Frame* outframe, int* current, int* map)
{
	int j, k, l;
	int found = 0;	
	int num_matches = 0;
	int min_loc = -1;
	int tie = 0;
	int min_val = 99;
	
	//see if type is determined 
	//printf("in determine_type() type is %d for *map %d\n", outframe->sites[*map].type, *map);
	if(outframe->sites[*map].type == -1)
		{
		//see how many matches remain
		//printf("num_cg_types is %d\n", control->num_cg_types);
		for(j = 0; j < control->num_cg_types; j++) //loop through all prototypes
			{			
			//printf("test type %d with match value %d\n", j, outframe->sites[*map].matches[j]);
			if(outframe->sites[*map].matches[j] == 1) //only consider if they are still in the running
				{
				//printf("possible match for type %d\n", j);
				
				//verify if newest type is still a match
				found = 0;
				
				if( control->prototype[j].num == min_val) 
					{
					tie++;
					}
				if( control->prototype[j].num < min_val)
					{
					min_loc = j;
					min_val = control->prototype[j].num;
					tie = 0;
					}
				
				//printf("prototype[%d].num is %d\n", j, control->prototype[j].num);
				for(k = 0; k < control->prototype[j].num; k++) //check all prototypes types
					{	
					//printf("k is %d: compare inframe->atoms[*current].type  %d ==  %d control->prototype[j].num_list[k]	\n", k, inframe->atoms[*current].type, control->prototype[j].num_list[k]);	
					if(inframe->atoms[*current].type == control->prototype[j].num_list[k])
						{
						found=1;
						break;
						}
					}
				}	
				
			if(found == 0)  
				{
				outframe->sites[*map].matches[j] = 0;
				}
			else 
				{
				num_matches++;
				}
			}
			
		//printf("%d matches possible\n", num_matches);
				
		//see if we have a unique match (or need evaluation or no match)
		if(num_matches > 1) //we need to force evaluation in tie-break
			{
			//printf("force evaluation for %d matches\n", num_matches);
			//check smallest
			if( outframe->sites[*map].num_in_site > (min_val - 1) ) //we have to fail because the molecule is too large
				{
				outframe->sites[*map].matches[min_loc] = 0;
				
				//see if we need to do this for others
				if(tie > 0)
					{
					for(j = min_loc; j < control->num_cg_types; j++)
						{
						if(outframe->sites[*map].matches[j] == min_val)
							{
							outframe->sites[*map].matches[j] = 0;
							num_matches--;
							}
						}
					}
				}
			
			if(outframe->sites[*map].num_in_site == (min_val - 1)) //actual situations may arise where we need to use both this and the previous statements
				{
				//evaluate to determine if total match
				//declare and initialize tesiting prototype
				int evaluate[ control->prototype[min_loc].num ];
				for(j = 0; j < control->prototype[min_loc].num; j++) evaluate[j] = -1;
				l = 1;
				
				//loop through site-type and record matches
				for(j = 0; j < outframe->sites[*map].num_in_site; j++) //j is for site types
					{
					//try to match with each prototype type
					for(k = 0; k < control->prototype[min_loc].num; k++)	//k is for prototype type
						{
						if( evaluate[k] == -1) //only consider if it is not yet claimed
							{
							if( outframe->sites[*map].coord[j].type == control->prototype[min_loc].num_list[k])
								{
								evaluate[k] = j;
								break;	//this leaves the for loop in k and moves on to next type
								}
							}
						//fail if we are still here at the last point and have not already last
						if( k == (control->prototype[min_loc].num - 1) ) l = 0;
						}
					}
				//also check new value
				if(l == 1)
					{
					for(k = 0; k < control->prototype[min_loc].num; k++)
						{
						if (evaluate[k] == -1)
							{
							if( inframe->atoms[*current].type == control->prototype[min_loc].num_list[k] )
								{
								evaluate[k] = outframe->sites[*map].num_in_site + 1;
								break;
								}
							}
						}
					}
				
				//if we have failed the evaluate
				if(l == 0)
					{
					outframe->sites[*map].matches[min_loc] = 0;
					num_matches--;
					}
					
				//what if it passes and there a tie?	
				//if there are no unknown sites (temporary pass)
				else if(l == 1)
					{
					//check that all prototypes entries are used
					for(k = 0; k < control->prototype[min_loc].num; k++)
						{
						if(evaluate[k] == -1) l = 0;
						}
									
					if(l == 0)
						{
						outframe->sites[*map].matches[min_loc] = 0;
						num_matches--;
						}
					else if(l == 1)
						{
						num_matches = 1; 
						//set all other matches values to zero
						for(k = 0; k < control->num_cg_types; k++)
							{
							outframe->sites[*map].matches[min_loc] = 0;
							}
						outframe->sites[*map].matches[min_loc] = 1;
						}
					}
				}
			}

		if(num_matches == 1) //we have a unique winner
			{
			//printf("unique winner\n");
			//find and set type
			for(j = 0; j < control->num_cg_types; j++)
				{
				if(outframe->sites[*map].matches[j] == 1)
					{
					outframe->sites[*map].type = j + 1;
					//if( (j + 1) > outframe->type_count) outframe->type_count = j + 1;
					outframe->type_num[j] += 1;
					
					//check all known types to see if this is a match
					l = 0;
					for(k = 0; k < outframe->type_count; k++) 
						{
						if(outframe->type[k] == outframe->sites[*map].type)
							{
							l = 1;
							}
						}
					
					if(l == 0)
						{
						outframe->type[outframe->type_count] = outframe->sites[*map].type;
						outframe->type_count++;
						}
					}
				}
			}
		else if(num_matches == 0)
			{
			//there is either a new molecule or an error in the inputs
			 printf("ERROR: NO MOLECULE TYPE MATCH !!! (in determine_type())\n");
			}				
		//} else {
		//	printf("DETERMINE_TYPE SAYS THAT TYPE IS DETERMINED!!!\n");
		}
}	

///////////////////////////////
///   compatible_type_test  ///
///////////////////////////////

int compatible_type_test(Controller* control, Frame* inframe, Frame* outframe, int current_type, int map)
{
	//printf("compatible_type_test for type %d\n", current_type);
	int j, k, l, m;
	int check = 0;
	
	//see how many matches remain
	for(j = 0; j < control->num_cg_types; j++) { //loop through all prototypes			
		//printf("type %d with matches of map %d for %d\n", j, map, outframe->sites[map].matches[j]);
		//printf("type %d with matches %d of map %d\n", j, outframe->sites[map].matches[j], map);
		//if(outframe->sites[map].matches[j] == 1) { //only consider if they are still in the running
			
			//see if the combined molecule is of reasonable size
			//printf("num_in_site is %d\n", outframe->sites[map].num_in_site);
			//printf("prototype %d has num %d\n", j, control->prototype[j].num);
			
			if( outframe->sites[map].num_in_site > control->prototype[j].num ) {
				//skip this type since it is impossible
				//printf("skip since this type is impossible\n");
				continue;
			}
			
			//see if this type is part of the potential prototype
			int evaluate[ control->prototype[j].num ];		
			check = 0;
			for(k = 0; k < control->prototype[j].num; k++) evaluate[k] = -1;
			
			//specifically test new type
			//printf("check evaluate with j %d  and num %d\n", j, control->prototype[j].num);
			
			for(k = 0; k < control->prototype[j].num; k++) {
				//printf(" k %d with current_type %d and num_list type %d\n", k, current_type, control->prototype[j].num_list[k]);
				if( current_type == control->prototype[j].num_list[k] )
				{
					evaluate[k] = outframe->sites[map].num_in_site + 1;
					check = 1;
					//printf("temp pass to check 1\n");
					return 1;
					break;
				}
			}
/*			
			printf("at if of evaluate with check %d\n", check);
			
			if (check == 0) {
				continue;
			} else {
			
				printf("else for check test\n");
				for(l = 0; l < outframe->sites[map].num_in_site; l++) { //l is for site types
					//try to match with each prototype type
					for(m = 0; m < control->prototype[j].num; m++) { //m is for prototype type
						if( evaluate[m] == -1) { //only consider if it is not yet claimed
						
							printf("  test COORD %d type %d vs num %d type %d\n", l, outframe->sites[map].coord[l].type, m, control->prototype[j].num_list[m]);
							if( outframe->sites[map].coord[l].type == control->prototype[j].num_list[m]) {
								evaluate[m] = l;
								break;	//this leaves the for loop in k and moves on to next type
							}
						}
						//fail if we are still here at the last point and have not already last
						if( m == (control->prototype[j].num - 1) ) {
							printf("inner inner fail\n");
							check = 0;
							break;
						}
					}
					if (check == 0) {
						printf("inner fail\n");
						break;
					}
				}
				if(check == 0) {
					continue;
				}
			
			//we have found a match that is still viable
			printf("pass\n");
			return 1;
			}
*/		//}
	}
	//otherwise, we have failed to find a match
	//printf("fail\n");
	return -1;
}	

/////////////////////////////
///   geometry_mapping	  ///
/////////////////////////////

void geometry_mapping(Controller* control, Frame* inframe, Frame* outframe)
{
	int i, j;
	double outx, outy, outz;
	double tot_mass;
	double dist, box;
	
	//apply averaging and processing
	if(control->geometry_map_flag == 0) //map to CoM
		{
		for(i = 0; i < outframe->num_atoms; i++)
			{
			outx = 0.0;
			outy = 0.0;
			outz = 0.0;
			tot_mass = 0.0;
				
			for(j = 0; j < outframe->sites[i].num_in_site; j++)
				{
				//check to see if coordinates are wrapped (and reset out* values)
				dist = outframe->sites[i].coord[j].x - outframe->sites[i].coord[0].x;
				box = outframe->xmax - outframe->xmin;
				if( abs(dist) >= (box/2.0) )
					{
					if(dist > 0) 
						{
						outframe->sites[i].coord[j].x -= box;
						}
					else
						{
						outframe->sites[i].coord[j].x += box;
						}
					}					
					
				dist = outframe->sites[i].coord[j].y - outframe->sites[i].coord[0].y;
				box = outframe->ymax - outframe->ymin;
				if( abs(dist) >= (box/2.0) )
					{
					if(dist > 0) 
						{
						outframe->sites[i].coord[j].y -= box;
						}
					else
						{
						outframe->sites[i].coord[j].y += box;
						}
					}	
					
				dist = outframe->sites[i].coord[j].z - outframe->sites[i].coord[0].z;
				box = outframe->zmax - outframe->zmin;
				if( abs(dist) >= (box/2.0) )
					{
					if(dist > 0) 
						{
						outframe->sites[i].coord[j].z -= box;
						}
					else
						{
						outframe->sites[i].coord[j].z += box;
						}
					}
						
				//add value to weighted sum
				outx += outframe->sites[i].coord[j].x * outframe->sites[i].coord[j].mass;
				outy += outframe->sites[i].coord[j].y * outframe->sites[i].coord[j].mass;
				outz += outframe->sites[i].coord[j].z * outframe->sites[i].coord[j].mass;
				tot_mass += outframe->sites[i].coord[j].mass;
				}
				
			//calculate average positions
			//if(i == 1) printf("outx = %lf, tot_mass = %lf from outframe->sites[%d].num_in_site = %d\n", outx, tot_mass, i, outframe->sites[i].num_in_site);
			outframe->sites[i].x = outx / tot_mass;
			outframe->sites[i].y = outy / tot_mass;
			outframe->sites[i].z = outz / tot_mass;

			//check if final coordinate is in box (wrap)
			if( outframe->sites[i].x > outframe->xmax ) 
				{
				outframe->sites[i].x -= outframe->xmax - outframe->xmin;
				}
			if( outframe->sites[i].x < outframe->xmin ) 
				{
				outframe->sites[i].x += outframe->xmax - outframe->xmin;
				}
			if( outframe->sites[i].y > outframe->ymax ) 
				{
				outframe->sites[i].y -= outframe->ymax - outframe->ymin;
				}
			if( outframe->sites[i].y < outframe->ymin ) 
				{
				outframe->sites[i].y += outframe->ymax - outframe->ymin;
				}
			if( outframe->sites[i].z > outframe->zmax ) 
				{
				outframe->sites[i].z -= outframe->zmax - outframe->zmin;
				}
			if( outframe->sites[i].z < outframe->zmin ) 
				{
				outframe->sites[i].z += outframe->zmax - outframe->zmin;
				}
					
			//set total mass
			outframe->sites[i].mass = tot_mass;
			}							
		}
	else if(control->geometry_map_flag == 1) //map to CoG
		{
		for(i = 0; i < outframe->num_atoms; i++)
			{
			outx = 0.0;
			outy = 0.0;
			outz = 0.0;
			tot_mass = 0.0;
				
			for(j = 0; j < outframe->sites[i].num_in_site; j++)
				{
				//check to see if coordinates are wrapped (and reset out* values)
				dist = outframe->sites[i].coord[j].x - outframe->sites[i].coord[0].x;
				box = outframe->xmax - outframe->xmin;
				if( abs(dist) >= (box/2.0) )
					{
					if(dist > 0) 
						{
						outframe->sites[i].coord[j].x -= box;
						}
					else
						{
						outframe->sites[i].coord[j].x += box;
						}
					}
				
				dist = outframe->sites[i].coord[j].y - outframe->sites[i].coord[0].y;
				box = outframe->ymax - outframe->ymin;
				if( abs(dist) >= (box/2.0) )
					{
					if(dist > 0) 
						{
						outframe->sites[i].coord[j].y -= box;
						}
					else
						{
						outframe->sites[i].coord[j].y += box;
						}
					}
				
				dist = outframe->sites[i].coord[j].z - outframe->sites[i].coord[0].z;
				box = outframe->zmax - outframe->zmin;
				if( abs(dist) >= (box/2.0) )
					{
					if(dist > 0) 
						{
						outframe->sites[i].coord[j].z -= box;
						}
					else
						{
						outframe->sites[i].coord[j].z += box;
						}
					}
				
				//add value to weighted sum
				outx += outframe->sites[i].coord[j].x;
				outy += outframe->sites[i].coord[j].y;
				outz += outframe->sites[i].coord[j].z;
				tot_mass += outframe->sites[i].coord[j].mass;
				}
			
			//calculate average positions
			outframe->sites[i].x = outx / outframe->sites[i].num_in_site;
			outframe->sites[i].y = outy / outframe->sites[i].num_in_site;
			outframe->sites[i].z = outz / outframe->sites[i].num_in_site;
			
			//check if final coordinate is in box (wrap)
			if( outframe->sites[i].x > outframe->xmax ) 
				{
				outframe->sites[i].x -= outframe->xmax - outframe->xmin;
				}
			if( outframe->sites[i].x < outframe->xmin ) 
				{
				outframe->sites[i].x += outframe->xmax - outframe->xmin;
				}
			if( outframe->sites[i].y > outframe->ymax ) 
				{
				outframe->sites[i].y -= outframe->ymax - outframe->ymin;
				}
			if( outframe->sites[i].y < outframe->ymin ) 
				{
				outframe->sites[i].y += outframe->ymax - outframe->ymin;
				}
			if( outframe->sites[i].z > outframe->zmax ) 
				{
				outframe->sites[i].z -= outframe->zmax - outframe->zmin;
				}
			if( outframe->sites[i].z < outframe->zmin ) 
				{
				outframe->sites[i].z += outframe->zmax - outframe->zmin;
				}
					
			//set total mass
			outframe->sites[i].mass = tot_mass;
			}
		}
}

//////////////////////////////////////////////
///   free_sensitivity_intermediates	  ///
/////////////////////////////////////////////

void free_sensitivity_intermediates(Controller* control, Frame* map1, Frame* map2)
{
	int i;
	
	//free all data for atoms
	for(i = 0; i < map1->num_atoms; i++)
		{
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
	
	
	for(i = 0; i < control->num_files; i++)
		{
	
		//free all data for atoms
		for(j = 0; j < mapframes[0].num_atoms; j++)
			{
			free(mapframes[i].sites[j].observables);	
			free(mapframes[i].sites[j].coord);
			free(mapframes[i].sites[j].matches);
			}
		
		free(mapframes[i].sites);
		
		free(mapframes[i].type);
		free(mapframes[i].type_num);
		}
}