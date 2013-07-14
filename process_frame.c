//process_frame.c

#include "headers.h"

void process_frame(Controller*, Frame*, Frame*);
void process_frames_and_log(Controller*, Frame*, Frame*, Frame*);
void map_all_atoms(Controller*, Frame*, Frame*); 
void map_some_atoms(Controller*, Frame*, Frame*);
void combine_sensitivity_data(Controller*, Frame*, Frame*, Frame*);
void free_sensitivity_intermediates(Controller*, Frame*, Frame*);

//////////////////////////
///   process_frame	  ///
/////////////////////////

void process_frame(Controller* control, Frame* inframe, Frame* outframe)
{
	int i, j;
	outframe->type_count = 0;
	
	printf("processing frame at step %d\n", control->frame);
	
	if( (control->frame == 1) || (control->sensitivity_flag != 0) )
		{
		//printf("allocate types\n");
		outframe->type = malloc(control->num_cg_types * sizeof(int));
		outframe->type_num = malloc(control->num_cg_types * sizeof(int));
		}
		
	for(i = 0; i< control->num_cg_types; i++)
		{
		//printf("initalize types\n");
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
		//printf("allocate sites for outframe to size %d\n", control->num_cg_sites);
		outframe->num_atoms = control->num_cg_sites;
		outframe->sites = malloc(outframe->num_atoms * sizeof(SITE));
		
		//also need to allocate observables
		for(i=0; i < outframe->num_atoms; i++)
			{
			outframe->sites[i].observables = malloc(outframe->num_observables * sizeof(double));
			outframe->sites[i].coord = malloc(control->max_to_map *sizeof(COORD));
			}
		//printf("observables and coord space allocated\n");
		}
	
	//intialize type array
	for(i = 0; i < control->num_cg_types; i++)
		{
		outframe->type[i] = -1;
		}

	//reset frame info
	for(i = 0; i < outframe->num_atoms; i++)
		{
		outframe->sites[i].q = 0.0;
		outframe->sites[i].num_in_site = 0;
		
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

//////////////////////////
///   map_all_atoms	  ///
/////////////////////////

void map_all_atoms(Controller* control, Frame* inframe, Frame* outframe)
{
	int i, j;
	int mol_count = 0;
	int mol_val;
	int site_count;
	double outx, outy, outz, tot_mass;
	double dist, box;
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
			//check to see if mol key is set
			mol_val = inframe->atoms[i].mol - 1;
			//printf("key[%d] = %d\n", mol_val, key[mol_val]);
			if(key[mol_val] == -1)
				{
				//printf("set information on key[%d] to %d\n", mol_val, mol_count);
				key[mol_val] = mol_count;
				outframe->sites[key[mol_val]].id = key[mol_val] + 1;
				outframe->sites[key[mol_val]].mol = key[mol_val] + 1;
				
				//also set molecule type by looking to match first id
				outframe->sites[key[mol_val]].type = -1;
				for(j = 0; j < outframe->type_count; j++)
					{
					if( outframe->type[j] == inframe->atoms[i].type) outframe->sites[key[mol_val]].type = j + 1;
					}
				//assign type if not set before
				if(outframe->sites[key[mol_val]].type == -1)
					{
					outframe->sites[key[mol_val]].type = outframe->type_count + 1;
					outframe->type[outframe->type_count] = inframe->atoms[i].type;
					outframe->type_count++;
					//printf("TYPE_COUNT IS %d ON TYPE %d\n", outframe->type_count, inframe->atoms[i].type);
					}
				mol_count++;
				outframe->type_num[ outframe->sites[key[mol_val]].type - 1]++;
				}					
		
			//printf("transfer information to CG site\n");
			//transfer information to CG site
			site_count = outframe->sites[ key[mol_val] ].num_in_site;
			//printf("site_count is %d\n", site_count);
		 
			outframe->sites[key[mol_val]].coord[site_count].x = inframe->atoms[i].x;
			outframe->sites[key[mol_val]].coord[site_count].y = inframe->atoms[i].y;
			outframe->sites[key[mol_val]].coord[site_count].z = inframe->atoms[i].z;
			outframe->sites[key[mol_val]].coord[site_count].mass = inframe->atoms[i].mass;
			outframe->sites[key[mol_val]].q += inframe->atoms[i].q;
			//printf("outframe->sites[key[%d]].q = %lf\n", mol_val, outframe->sites[key[mol_val]].q);
			
			for(j = 0; j < outframe->num_observables; j++)
				{
				outframe->sites[key[mol_val]].observables[j] += inframe->atoms[i].observables[j];
				//printf("observable[%d] is %lf\n", j, outframe->sites[key[mol_val]].observables[j]);
				}
			//printf("finished observable transfer at num_in_site %d\n", outframe->sites[ key[mol_val] ].num_in_site);
			outframe->sites[ key[mol_val] ].num_in_site++;
			}
				
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
				outframe->sites[i].x = outx / tot_mass;
				outframe->sites[i].y = outy / tot_mass;
				outframe->sites[i].z = outz / tot_mass;
				//printf("done finding center of mass for molecules %d \n", i);
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
		}
		outframe->num_mol = mol_count;		
}

/////////////////////////////
///   map_some_atoms	  ///
/////////////////////////////

void map_some_atoms(Controller* control, Frame* inframe, Frame* outframe)
{
	int i, j;
	int mol_count = 0;
	int mol_val;
	int site_count;
	double outx, outy, outz, tot_mass;
	double dist, box;
	int key[inframe->num_mol];
	
	//intialize key
	for(i = 0; i < inframe->num_mol; i++)
		{
		key[i] = -1;
		}

	int go_flag;
		
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

				//also set molecule type by looking to match first id
				outframe->sites[key[mol_val]].type = -1;
				for(j = 0; j < outframe->type_count; j++)
					{
					if( outframe->type[j] == inframe->atoms[i].type) outframe->sites[key[mol_val]].type = j + 1;
					}
				//assign type if not set before
				if(outframe->sites[key[mol_val]].type == -1)
					{
					outframe->sites[key[mol_val]].type = outframe->type_count + 1;
					outframe->type[outframe->type_count] = inframe->atoms[i].type;
					outframe->type_count++;
					//printf("TYPE_COUNT IS %d ON TYPE %d\n", outframe->type_count, inframe->atoms[i].type);
					}
				mol_count++;
				outframe->type_num[ outframe->sites[key[mol_val]].type  - 1]++;
				}					
				
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
		}
		
	outframe->num_mol = mol_count;		
}

void combine_sensitivity_data(Controller* control, Frame* inframe1, Frame* inframe2, Frame* outframe)
{
	//declare variables
	int i, j;
	double scalar = control->log_value - control->guess;

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
	outframe->type = malloc(control->num_cg_types * sizeof(int));
	outframe->type_num = malloc(control->num_cg_types * sizeof(int));

	for(i = 0; i< control->num_cg_types; i++)
		{
		outframe->type[i] = inframe1->type[i];
		outframe->type_num[i] = inframe1->type_num[i];
		}

	//allocate site space
	outframe->num_atoms = control->num_cg_sites;
	outframe->sites = malloc(outframe->num_atoms * sizeof(SITE));		
	for(i=0; i < outframe->num_atoms; i++)
		{
		outframe->sites[i].observables = malloc(outframe->num_observables * sizeof(double));
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
		
		//combine "observable information"
		for(j = 0; j < outframe->num_observables; j++)
			{
			outframe->sites[i].observables[j] = inframe2->sites[i].observables[j] - scalar * inframe1->sites[i].observables[j];
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
		}
		
	free(map1->sites);
	free(map2->sites);
	
	free(map1->type);
	free(map2->type);
	
	free(map1->type_num);
	free(map2->type_num);
}