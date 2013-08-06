//process_frame.c

#include "headers.h"

//original functions
void process_frame(Controller*, Frame*, Frame*);
void process_frames_and_log(Controller*, Frame*, Frame*, Frame*);
void process_no_map_frames_and_log(Controller*, Frame*, Frame*, Frame*);
void map_all_atoms(Controller*, Frame*, Frame*); 
void map_some_atoms(Controller*, Frame*, Frame*);
//sensitivity functions
void combine_sensitivity_data(Controller*, Frame*, Frame*, Frame*);
void combine_cg_sensitivity_data(Controller*, Frame*, Frame*, Frame*);
void process_charge_frames(Controller*, Frame*, Frame*);
void process_charge_log(Controller*);
//slave functions (newer modularization)
void determine_type(Controller*, Frame*, Frame*, int*, int*);
void geometry_mapping(Controller*, Frame*, Frame*);
void composite_charge_frame(Controller*, Frame*, Frame*, int, int*, int*);
void composite_charge_log(Controller*, int, int*, int*);
//key/hash functions
void key_lookup(int, int, int, int*);
void reverse_key_lookup(int, int, int*, int*);
//free functions
void free_sensitivity_intermediates(Controller*, Frame*, Frame*);
void free_charge_intermediates(Controller*, Frame*);

//////////////////////////
///   process_frame	  ///
/////////////////////////

void process_frame(Controller* control, Frame* inframe, Frame* outframe)
{
	//printf("enter process_frame\n");
	int i, j;
	outframe->type_count = 0;
	
	//printf("process_frame: inframe->atoms[1].x = %lf\n", inframe->atoms[1].x);
	//printf("processing frame at step %d\n", control->frame);
	
	if( (control->frame == 1) || (control->sensitivity_flag == 1) )
		{
		//printf("allocate types\n");
		outframe->type = malloc(control->num_cg_types * sizeof(int));
		outframe->type_num = malloc(control->num_cg_types * sizeof(int));
		}
	
	//intitalize type array
	//printf("intializing types\n");
	for(i = 0; i< control->num_cg_types; i++)
		{
		//printf("initalize type %d\n", i);
		outframe->type[i] = -1;
		outframe->type_num[i] = 0;
		}

	//printf("start process frame\n");
	
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
	
	//printf("accumulate information from fine_grained sites as we go for frame at timestep %d\n", inframe->timestep);
	//accumulate information from fine_grained sites as we go (assuming molecule to molecule mapping)

	if(control->map_style_flag == 0) 
		{
		map_all_atoms(control, inframe, outframe);
		}
	else if(control->map_style_flag == 1)
		{
		map_some_atoms(control, inframe, outframe);
		}
		
	//printf("proces_frames: outframe->sites[1].x = %lf\n", outframe->sites[1].x);
}

/////////////////////////////////
///   process_frames_and_log  ///
/////////////////////////////////

void process_frames_and_log(Controller* control, Frame* inframe1, Frame* inframe2, Frame* outframe)
{
	//decclare Frames for intermediates mapped
	Frame map1;
	Frame map2;
	map1.num_atoms = 0;
	map2.num_atoms = 0;
	
	//do normal mapping for the 2 frames
	//printf("PROCESS_FRAME #1 called\n");
	process_frame(control, inframe1, &map1);
	//printf("PROCESS_FRAME #2 called\n");
	process_frame(control, inframe2, &map2);
	
	//printf("COMBINATION called\n");
	//do combination of files to get final output
	combine_sensitivity_data(control, &map1, &map2, outframe);

	//printf("FREE INTERMEDIATES\n");
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
	//allocate space
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
		//printf("initalize type %d\n", i);
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
		for(j = 0; j <= inframe1->type_count; i++)
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
		for(j = 0; j <= inframe2->type_count; i++)
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
		
	//printf("COMBINATION called\n");
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
	int mol_val;
	int site_count;
	int key[inframe->num_mol];
	
	//intialize key
	for(i = 0; i < inframe->num_mol; i++)
		{
		key[i] = -1;
		}

	//printf("map_all_atoms: inframe->atoms[1].x = %lf\n", inframe->atoms[1].x);

	
	if(control->observable_map_flag == 0) //sum obesrvables
		{
		//process and sort all FG atoms
		for(i = 0; i < inframe->num_atoms; i++)
			{
			//check to see if mol key is set
			mol_val = inframe->atoms[i].mol - 1;
			//printf("key[%d] = %d\n", mol_val, key[mol_val]);
			if(key[mol_val] == -1)
				{
				//printf("set information on key[%d] to %d\n", mol_val, mol_count);
				key[mol_val] = mol_count;
				outframe->sites[key[mol_val]].id = key[mol_val] + 1;
				outframe->sites[key[mol_val]].mol = key[mol_val] + 1;
				mol_count++;
				}					
							
			determine_type(control, inframe, outframe, &i, &key[mol_val]);
			
			//if(i == 1)
			//printf("transfer information to CG site\n");
			//transfer information to CG site
			site_count = outframe->sites[ key[mol_val] ].num_in_site;
			//printf("site_count is %d\n", site_count);
		 
			outframe->sites[key[mol_val]].coord[site_count].x = inframe->atoms[i].x;
			outframe->sites[key[mol_val]].coord[site_count].y = inframe->atoms[i].y;
			outframe->sites[key[mol_val]].coord[site_count].z = inframe->atoms[i].z;
			outframe->sites[key[mol_val]].coord[site_count].mass = inframe->atoms[i].mass;
			outframe->sites[key[mol_val]].coord[site_count].type = inframe->atoms[i].type;
			outframe->sites[key[mol_val]].q += inframe->atoms[i].q;
			outframe->sites[ key[mol_val] ].num_in_site++;
			//if(i == 1) printf("num_in_site for 1 is %d\n", outframe->sites[key[mol_val]].num_in_site);

			//printf("outframe->sites[key[%d]].q = %lf\n", mol_val, outframe->sites[key[mol_val]].q);
			
			for(j = 0; j < outframe->num_observables; j++)
				{
				outframe->sites[key[mol_val]].observables[j] += inframe->atoms[i].observables[j];
				//printf("observable[%d] is %lf\n", j, outframe->sites[key[mol_val]].observables[j]);
				}
			//printf("finished observable transfer at num_in_site %d\n", outframe->sites[ key[mol_val] ].num_in_site);
			
			}
				
		geometry_mapping(control, inframe, outframe);
		
		}
	outframe->num_mol = mol_count;	
		
	//printf("map_all_atoms: outframe->sites[1].x = %lf\n", outframe->sites[1].x);

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
		
	if(control->observable_map_flag == 0) //sum obesrvables
		{
		//process and sort all FG atoms
		for(i = 0; i < inframe->num_atoms; i++)
			{
			
			//see if molecule is in include list
			go_flag = 0;
			for(j = 0; j < control->num_map; j++)
				{
				if(control->map[j] == inframe->atoms[i].id)
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
				//printf("set information on key[%d] to %d\n", mol_val, mol_count);
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
			outframe->sites[key[mol_val]].q += inframe->atoms[i].q;
			
			for(j = 0; j < outframe->num_observables; j++)
				{
				outframe->sites[key[mol_val]].observables[j] += inframe->atoms[i].observables[j];
				}
			//printf("finished observable transfer at num_in_site %d\n", outframe->sites[ key[mol_val] ].num_in_site);
			outframe->sites[ key[mol_val] ].num_in_site++;
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

	//printf("in combine\n");
	int paired[inframe2->num_atoms];
	
	for(i = 0; i < inframe2->num_atoms; i++) paired[i] = 0;

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
	
	//printf("allocate type with size %d and %d\n", control->num_cg_types, control->num_cg_types);
	//allocate type space and initalize
	if(control->frame == 1)
		{
		outframe->type = malloc(control->num_cg_types * sizeof(int));
		outframe->type_num = malloc(control->num_cg_types * sizeof(int));
		}
		
	//printf("copy types\n");
	outframe->type_count = inframe1->type_count;
	for(i = 0; i< control->num_cg_types; i++)
		{
		outframe->type[i] = inframe1->type[i];
		outframe->type_num[i] = inframe1->type_num[i];
		}

	//printf("allocate combine space\n");
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

	//printf("start conversion to outframe\n");
	//copy frame info on basic info
	for(i = 0; i < outframe->num_atoms; i++)
		{
		//found = 0;
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
			//printf("testing pairs as %d %d at %lf = %lf vs %lf\n", i, j, outframe->sites[i].x, inframe1->sites[j].x, inframe2->sites[j].x);
			//check each position and skip rest if failed (allow for some rounding error in mapping)
			if( abs(outframe->sites[i].x - inframe2->sites[j].x) > 0.001 ) continue;
			if( abs(outframe->sites[i].y - inframe2->sites[j].y) > 0.001 ) continue;
			if( abs(outframe->sites[i].z - inframe2->sites[j].z) > 0.001 ) continue;
			
			//for matched sites
			paired[j] = 1;
			//found = 1;
			//combine "observable information"
			for(k = 0; k < outframe->num_observables; k++)
				{
				outframe->sites[i].observables[k] = inframe2->sites[j].observables[k] + control->sign_flag * scalar * inframe1->sites[i].observables[k];
				outframe->sites[i].observables[k] /= control->scaleF;
				//printf("pair matched for i = %d j = %d k = %d\n", i, j, k);
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

	//printf("in combine\n");
	int paired[inframe2->num_atoms];
	
	for(i = 0; i < inframe2->num_atoms; i++) paired[i] = 0;

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
	
	//printf("allocate type with size %d and %d\n", control->num_cg_types, control->num_cg_types);
	//allocate type space and initalize
	if(control->frame == 1)
		{
		outframe->type = malloc(control->num_cg_types * sizeof(int));
		outframe->type_num = malloc(control->num_cg_types * sizeof(int));
		}
		
	//printf("copy types\n");
	outframe->type_count = inframe1->type_count;
	for(i = 0; i< control->num_cg_types; i++)
		{
		outframe->type[i] = inframe1->type[i];
		outframe->type_num[i] = inframe1->type_num[i];
		}

	//printf("allocate combine space\n");
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

	//printf("start conversion to outframe\n");
	//copy frame info on basic info
	for(i = 0; i < outframe->num_atoms; i++)
		{
		//found = 0;
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
			//printf("testing pairs as %d %d at %lf = %lf vs %lf\n", i, j, outframe->sites[i].x, inframe1->sites[j].x, inframe2->sites[j].x);
			//check each position and skip rest if failed (allow for some rounding error in mapping)
			if( abs(outframe->sites[i].x - inframe2->atoms[j].x) > 0.001 ) continue;
			if( abs(outframe->sites[i].y - inframe2->atoms[j].y) > 0.001 ) continue;
			if( abs(outframe->sites[i].z - inframe2->atoms[j].z) > 0.001 ) continue;
			
			//for matched sites
			paired[j] = 1;
			//found = 1;
			//combine "observable information"
			for(k = 0; k < outframe->num_observables; k++)
				{
				outframe->sites[i].observables[k] = inframe2->atoms[j].observables[k] - scalar * inframe1->atoms[i].observables[k];
				//printf("pair matched for i = %d j = %d k = %d\n", i, j, k);
				}
			break;
			}
		}
}

/////////////////////////////////
///   process_charge_frames	 ///
////////////////////////////////

void process_charge_frames(Controller* control, Frame* inframes, Frame* outframes)
{
	//declare variables
	int i, j;
	int self;
	int temp, count;
	int* own = malloc( (control->num_charges) * sizeof(int) );
	int* mixed;
	if(control->num_charges > 1) mixed = malloc( (control->num_charges - 1) * sizeof(int) );
	else mixed = malloc( sizeof(int) );
	
	//printf("process_charge_frames: inframes[1].atoms[1].x = %lf\n", inframes[1].atoms[1].x);
	//create temporary frames
	Frame* tempframes = malloc(control->num_files * sizeof(Frame));
	
	//map each frame in turn
	for(i = 0; i < control->num_files; i++)
		{
		//printf("process frame %d\n", i);
		//set-up each outframe as necessary
		tempframes[i].num_atoms = 0;
		
		//process into CG
		process_frame(control, &inframes[i], &tempframes[i] );
		}
	
	//printf("set-up inteactions\n");
	//set-up own-array for like-like interactions
	for(i = 0; i < control->num_charges; i++) own[i] = i;
	
	//create composite outputs for each charge
	for(i = 0; i < control->num_charges; i++)
		{
		//set-up key-list for each run 
		self = i;
		count = 0;
		for(j = 0; j < control->num_charges; j++)
			{
			if(i != j)
				{
				key_lookup(self, j, control->num_charges, &temp);
				//printf("self = %d, j = %d, size=%d, mixed id is %d\n", self, j, control->num_charges, temp);
				mixed[count] = temp;
				//mixed_charge[count] = control->charge[i];
				count++;
				}
			}
		//own array does not change
		
		//call functions with pre-set templates
		//printf("composite!\n");
		composite_charge_frame(control, tempframes, &outframes[i], self, own, mixed);
		}
		
	//free template information
	free(own);
	free(mixed);
	//free intermediate/temporary frames
	free_charge_intermediates(control, tempframes);
	free(tempframes);
}

/////////////////////////////////
///   process_charge_log	 ///
////////////////////////////////

void process_charge_log(Controller* control)
{
		//declare variables
	int i, j;
	int self;
	int temp, count;
	int* own = malloc( (control->num_charges) * sizeof(int) );
	int* mixed;
	if(control->num_charges > 1) mixed = malloc( (control->num_charges - 1) * sizeof(int) );
	else mixed = malloc( sizeof(int) );
	
	//printf("set-up inteactions\n");
	//set-up own-array for like-like interactions
	for(i = 0; i < control->num_charges; i++) own[i] = i;
	
	//create composite outputs for each charge
	for(i = 0; i < control->num_charges; i++)
		{
		//set-up key-list for each run 
		self = i;
		count = 0;
		for(j = 0; j < control->num_charges; j++)
			{
			if(i != j)
				{
				key_lookup(self, j, control->num_charges, &temp);
				//printf("self = %d, j = %d, size=%d, mixed id is %d\n", self, j, control->num_charges, temp);
				mixed[count] = temp;
				//mixed_charge[count] = control->charge[i];
				count++;
				}
			}
		//own array does not change
		
		//call functions with pre-set templates
		//printf("composite!\n");
		composite_charge_log(control, self, own, mixed);
		}
		
	//free template information
	free(own);
	free(mixed);	
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
	if(outframe->sites[*map].type == -1)
		{
		//see how many matches remain
		for(j = 0; j < control->num_cg_types; j++) //loop through all prototypes
			{
			
			if(outframe->sites[*map].matches[j] == 1) //only consider if they are still in the running
				{
				//verify if newest type is still a match
				found = 0;
					
				if( control->prototype[j].num == min_val) tie++;
				if( control->prototype[j].num < min_val)
					{
					min_loc = j;
					min_val = control->prototype[j].num;
					tie = 0;
					}
				
				for(k = 0; k < control->prototype[j].num; k++) //check all prototypes types
					{
					
					if(inframe->atoms[*current].type == control->prototype[j].num_list[k])
						{
						found=1;
						break;
						}
					}
				}	
				
			if(found == 0)  outframe->sites[*map].matches[j] = 0;
			else num_matches++;
			}
					
		//see if we have a unique match (or need evaluation or no match)
		if(num_matches > 1) //we need to force evaluation in tie-break
			{
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
				
				//if there are no unkown sites (temporary pass)
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
						//printf("outframe->type count is now %d after %d atoms read\n", outframe->type_count, *current);
						}
					}
				}
			}
		else if(num_matches == 0)
			{
			//there is either a new molecule or an error in the inputs
			 printf("ERROR: NO MOLECULE TYPE MATCH !!!\n");
			}
				
		}
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
	
	//printf("\ngeometry_mapping: inframe->atoms[1].x = %lf\n", inframe->atoms[1].x);
	//printf("geometry_mapping: outframe->sites[1].coord[1].x = %lf\n", outframe->sites[1].coord[1].x);
	
	//apply averaging and processing
	if(control->geometry_map_flag == 0) //map to CoM
		{
		//printf("COM\n");
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
					if(dist > 0) outframe->sites[i].coord[j].x -= box;
					else		 outframe->sites[i].coord[j].x += box;
					}					
					
				dist = outframe->sites[i].coord[j].y - outframe->sites[i].coord[0].y;
				box = outframe->ymax - outframe->ymin;
				if( abs(dist) >= (box/2.0) )
					{
					if(dist > 0) outframe->sites[i].coord[j].y -= box;
					else		 outframe->sites[i].coord[j].y += box;
					}	
					
				dist = outframe->sites[i].coord[j].z - outframe->sites[i].coord[0].z;
				box = outframe->zmax - outframe->zmin;
				if( abs(dist) >= (box/2.0) )
					{
					if(dist > 0) outframe->sites[i].coord[j].z -= box;
					else		 outframe->sites[i].coord[j].z += box;
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
			if( outframe->sites[i].x > outframe->xmax ) outframe->sites[i].x -= outframe->xmax - outframe->xmin;
			if( outframe->sites[i].x < outframe->xmin ) outframe->sites[i].x += outframe->xmax - outframe->xmin;
			if( outframe->sites[i].y > outframe->ymax ) outframe->sites[i].y -= outframe->ymax - outframe->ymin;
			if( outframe->sites[i].y < outframe->ymin ) outframe->sites[i].y += outframe->ymax - outframe->ymin;
			if( outframe->sites[i].z > outframe->zmax ) outframe->sites[i].z -= outframe->zmax - outframe->zmin;
			if( outframe->sites[i].z < outframe->zmin ) outframe->sites[i].z += outframe->zmax - outframe->zmin;
					
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
					if(dist > 0) outframe->sites[i].coord[j].x -= box;
					else		 outframe->sites[i].coord[j].x += box;
					}
				
				dist = outframe->sites[i].coord[j].y - outframe->sites[i].coord[0].y;
				box = outframe->ymax - outframe->ymin;
				if( abs(dist) >= (box/2.0) )
					{
					if(dist > 0) outframe->sites[i].coord[j].y -= box;
					else		 outframe->sites[i].coord[j].y += box;
					}
				
				dist = outframe->sites[i].coord[j].z - outframe->sites[i].coord[0].z;
				box = outframe->zmax - outframe->zmin;
				if( abs(dist) >= (box/2.0) )
					{
					if(dist > 0) outframe->sites[i].coord[j].z -= box;
					else		 outframe->sites[i].coord[j].z += box;
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
			if( outframe->sites[i].x > outframe->xmax ) outframe->sites[i].x -= outframe->xmax - outframe->xmin;
			if( outframe->sites[i].x < outframe->xmin ) outframe->sites[i].x += outframe->xmax - outframe->xmin;
			if( outframe->sites[i].y > outframe->ymax ) outframe->sites[i].y -= outframe->ymax - outframe->ymin;
			if( outframe->sites[i].y < outframe->ymin ) outframe->sites[i].y += outframe->ymax - outframe->ymin;
			if( outframe->sites[i].z > outframe->zmax ) outframe->sites[i].z -= outframe->zmax - outframe->zmin;
			if( outframe->sites[i].z < outframe->zmin ) outframe->sites[i].z += outframe->zmax - outframe->zmin;
					
			//set total mass
			outframe->sites[i].mass = tot_mass;
			}
		}
	//printf("geometry_mapping: outframe->sites[1].x = %lf\n", outframe->sites[1].x);

}

/////////////////////////////////////
///   composite_charge_frame	  ///
/////////////////////////////////////

void composite_charge_frame(Controller* control, Frame* mapframes, Frame* outframe, int self, int* own, int* mixed)//own is not actually needed assuming ordering
{
	//declare and initalize variables
	int i, j, k, l;
	int paired[mapframes[self].num_atoms];
	double scalar, key;
	for(i = 0; i < mapframes[self].num_atoms; i++) paired[i] = 0;

	//popluate outframe with necessary general information
	outframe->type_count = 0;
	outframe->xmin = mapframes[self].xmin;
	outframe->xmax = mapframes[self].xmax;
	outframe->ymin = mapframes[self].ymin;
	outframe->ymax = mapframes[self].ymax;
	outframe->zmin = mapframes[self].zmin;
	outframe->zmax = mapframes[self].zmax;
	outframe->timestep = mapframes[self].timestep;
	outframe->num_observables = mapframes[self].num_observables;
	
	//allocate type space and initalize
	outframe->type = malloc(control->num_cg_types * sizeof(int));
	outframe->type_num = malloc(control->num_cg_types * sizeof(int));

	//copy types
	for(i = 0; i< control->num_cg_types; i++)
		{
		outframe->type[i] = mapframes[self].type[i];
		outframe->type_num[i] = mapframes[self].type_num[i];
		}

	//allocate site space
	if(outframe->num_atoms != mapframes[self].num_atoms)
		{
		//printf("outframe->num_atoms set as %d compare to mapframes[0] of %d\n", outframe->num_atoms, mapframes[self].num_atoms);
		outframe->num_atoms = control->num_cg_sites;
		outframe->sites = malloc(outframe->num_atoms * sizeof(SITE));		
		
		for(i=0; i < outframe->num_atoms; i++)
			{
			outframe->sites[i].observables = malloc(outframe->num_observables * sizeof(double));
			}
		}
		
	scalar = 2.0 * control->charge[self];
	for(i = 0; i < control->num_charges; i++)
		{
		if(i != self) scalar -= control->charge[i];
		}
	
	//printf("copy basic info\n");
	//copy frame info on basic info
	for(i = 0; i < outframe->num_atoms; i++)
		{
		outframe->sites[i].num_in_site = 1;
		outframe->sites[i].id = mapframes[self].sites[i].id;
		outframe->sites[i].mol = mapframes[self].sites[i].mol;
		outframe->sites[i].type = mapframes[self].sites[i].type;
		outframe->sites[i].mass = mapframes[self].sites[i].mass;
		outframe->sites[i].q = mapframes[self].sites[i].q;
		outframe->sites[i].x = mapframes[self].sites[i].x;
		outframe->sites[i].y = mapframes[self].sites[i].y;
		outframe->sites[i].z = mapframes[self].sites[i].z;
	
		//ADD IN SELF TERM
		for(j = 0; j < outframe->num_observables; j++) outframe->sites[i].observables[j] = scalar * mapframes[self].sites[i].observables[j];
		
		}
		
	//printf("add OWN terms\n");
	//ADD IN OWN TERMS (as they are matched -- one frame source at a time)
	for(i = 0; i < control->num_charges; i++)
		{
		//only looking for not-self "own" terms
		if(i == self) continue;
		
		scalar = control->charge[i];
		
		//reinitalize paired array
		for(j = 0; j < outframe->num_atoms; j++) paired[j] = 0;
		
		//printf("position matching\n");
		//look for correct match based on position
		for(j = 0; j < outframe->num_atoms; j++) //loop over destination sites
			{
			//if (j < 10 && control->frame == 1) printf("at atom %d where own[%d] is %d\n", j, i, own[i]);
			for(k = 0; k < outframe->num_atoms; k++) //loop over potential mapframes site to find match
				{	
				//only consider sites that are not yet paired
				if(paired[k] == 1) continue;
				
				//if( k < 5 && control->frame == 1)
				//	{
				//	printf("\ntest\n");
				//	printf("outframe->sites[%d].x = %lf\n vs mapframes[%d]->sites[%d].x = %lf\n", j, outframe->sites[j].x, self, j, mapframes[self].sites[j].x);
				//	printf("mapframes[ %d ]->num_atoms = %d\n", own[i], mapframes[ own[i] ].num_atoms);
				//	printf("mapframes[ %d ]->sites[%d].x = %lf\n", own[i], k, mapframes[ own[i] ].sites[k].x);
				//	}
				//check each position and skip rest if failed (allow for some rounding error in mapping)
				if( abs(outframe->sites[j].x - mapframes[ own[i] ].sites[k].x) > 0.001 ) continue;
				if( abs(outframe->sites[j].y - mapframes[ own[i] ].sites[k].y) > 0.001 ) continue;
				if( abs(outframe->sites[j].z - mapframes[ own[i] ].sites[k].z) > 0.001 ) continue;
				//if(k < 10 && control->frame == 1) printf("passed position tests\n");
				//for matched sites
				paired[k] = 1;
				//found = 1;
				//combine "observable information"
				for(l = 0; l < outframe->num_observables; l++)
					{
					outframe->sites[j].observables[l] -= scalar * mapframes[ own[i] ].sites[k].observables[l];
					}
				break;
				}
			//if(k < 10 && control->frame == 1) if(paired[k] == 1) printf("new match found for %d at %d\n", j, k);
			}
		}
	
	//printf("add MIXED terms\n");
	//ADD IN MIXED TERMS (as they are matched -- one frame source at a time)
	for(i = 0; i < (control->num_charges - 1); i++)
		{
		//only looking for not-self "own" terms
		if(i == self) continue;
		
		//lookup pair from key
		reverse_key_lookup(mixed[i], control->num_charges, &k, &l);
		if(self == k) scalar = control->charge[l];
		else scalar = control->charge[k];
		
		//reinitalize paired array
		for(j = 0; j < outframe->num_atoms; j++) paired[j] = 0;
		//printf("look for position\n");
		//look for correct match based on position
		for(j = 0; j < outframe->num_atoms; j++) //loop over destination sites
			{
			for(k = 0; k < outframe->num_atoms; k++) //loop over potential mapframes site to find match
				{	
				//only consider sites that are not yet paired
				if(paired[k] == 1) continue;
				
				if( control->frame == 1)
					{
					//printf("\ntest: index values i=%d, j=%d, k=%d\n", i, j, k); 
					//printf("In MIX for charge index %d and base %d \n", mixed[i], self);
					//printf("outframe->sites[%d].x = %lf vs mapframes[%d]->sites[%d].x = %lf\n", j, outframe->sites[j].x, self, j, mapframes[self].sites[j].x);
					//printf("abs(distances) = %lf\n", abs(outframe->sites[j].x - mapframes[ mixed[i] ].sites[k].x));
					//printf("mapframes[ %d ]->sites[%d].x = %lf\n", own[i], k, mapframes[ own[i] ].sites[k].x);
					}
				//check each position and skip rest if failed (allow for some rounding error in mapping)
				if( abs(outframe->sites[j].x - mapframes[ mixed[i] ].sites[k].x) > 0.001 ) continue;	
				if( abs(outframe->sites[j].y - mapframes[ mixed[i] ].sites[k].y) > 0.001 ) continue;
				if( abs(outframe->sites[j].z - mapframes[ mixed[i] ].sites[k].z) > 0.001 ) continue;
			
			//for matched sites
				paired[k] = 1;
				//found = 1;
				//combine "observable information"
				for(l = 0; l < outframe->num_observables; l++)
					{
					outframe->sites[j].observables[l] += scalar * mapframes[ mixed[i] ].sites[k].observables[l];
					}
				break;
				}
			}
		}
	//printf("done MIX\n");
}

//////////////////////////////////
///   composite_charge_log	  ///
/////////////////////////////////

void composite_charge_log(Controller* control, int self, int* own, int* mixed)
{
	//declare and initalize variables
	int i, j, k, l;
	double scalar, key;
	
	//reset output value 
	control->guesses[self] = 0.0;
		
	//printf("copy basic info\n");
	//do SELF term
	scalar = 2.0 * control->charge[self];
	for(i = 0; i < control->num_charges; i++)
		{
		if(i != self) scalar -= control->charge[i];
		}

	control->guesses[self] += scalar * control->log_values[self];
		
	//printf("add OWN terms\n");
	//ADD IN OWN TERMS (as they are matched -- one frame source at a time)
	for(i = 0; i < control->num_charges; i++)
		{
		//only looking for not-self "own" terms
		if(i == self) continue;
		
		scalar = control->charge[i];
		control->guesses[self] -= scalar * control->log_values[ own[i] ];
		}
	
	//printf("add MIXED terms\n");
	//ADD IN MIXED TERMS (as they are matched -- one frame source at a time)
	for(i = 0; i < (control->num_charges - 1); i++)
		{
		//only looking for not-self "own" terms
		if(i == self) continue;
		
		//lookup pair from key
		reverse_key_lookup(mixed[i], control->num_charges, &k, &l);
		if(self == k) scalar = control->charge[l];
		else scalar = control->charge[k];
		
		control->guesses[self] += scalar * control->log_values[ mixed[i] ];
		}	
	//printf("done MIX\n");	
}

//////////////////////////////////////////////
///   key_lookup and reverse_key_lookup	  ///
/////////////////////////////////////////////

void key_lookup(int i, int j, int charge, int* key)
{
	int k;
	*key = 0;
	
	if(i == j) *key = i;
	else if(i < j)
		{
		for(k = 0; k <= i; k++)
			{
			*key += charge - k -1;
			}
		*key += j;
		}
	else if(i > j)
		{
		for(k = 0; k <= j; k++)
			{
			*key += charge - k -1;
			}
		*key += i;
		}
}

void reverse_key_lookup(int key, int charge, int* i, int* j)
{
	int k, temp, count;
	
	if( key < charge)	//must be self pairing
		{
		*i = key;
		*j = key;
		return;
		}
		
	//find smallest value for offset
	temp = key - charge;
	count = 0;
	while(temp > 0)
		{
		temp -= (charge - count);
		count++;
		}
	
	*i = count; //we have identified the smallest of the pair
	
	//find other value
	temp = key - 1;
	for(k = 0; k <= count; k++) temp -= (charge + k);
	*j = temp;
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