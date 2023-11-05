/***** includes *********************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <limits.h>

#define MAX_N_RESIDUES 9999999
/***** local structures *********************/
struct System { 
               int n_atoms;
               struct Atom *atoms;
               int n_residues;
               struct Residue *residues;
               struct Valid_res *valid_res;
               int n_chains;
               struct Chain *chains;
               int n_segments;
               struct Segment *segments;
               int n_models;
               struct Model *models;
              } System;

struct Atom {
/************** fields in PDB ********************/ 
		char at_name[5]; 
		char alt_loc;
		char res_name[4];
		char chain; 
                char element[3]; 
		int model;
                int at_n;
		int res_n;
		char res_ins; 
                double coor[3]; /* queste sostituiranno x, y, z */
                double ref[3];
		double occ,temp;
		char segid[5];
     		char pdb_chrg[3]; /* stringa con la carica formale in PDB */
/************* additional features **********************/
    		double mass;
    		char chem_group[5];  /* chemical group e.g. peptide methyl.... */
    		char group[5];  /* if sidechain or backbone or other */
    		double radius;
    		double charge;
    		int atom_type;
		int i;
	    } Atom;

struct Valid_res {
                 int val_aa;
                 int val_ca;
                 int val_atom_bb; /* = 0 if BB atom missing */
                 int val_atom_sc; /* = 0 if SC atom missing */
/* I segnuenti sono da migliorare */
                 int val_h; /* = 0 if a hydrogen is missing */
                 int val_bonds; /* = 0 if bonds > 10% wrt db */
                 int val_angles; /* = 0 if angles > 10% wrt db */
               } Valid_res;


struct Trj {
           int nf, naxf;  // number of frames, number of atoms per frame
           double **coor; 
           } Trj;

struct Residue { char res_name[4];
                 char chain;
                 char res_ins;
                 int res_n;
                 int model;
                 char seg_name[5];
                 int beg, end ; /*begin and end atoms of sorted atom list*/
                 int n_alt_loc;  
                 int prev, next;
                 double CA[3], CG[3], CM[3]; 
                 int res_type;
               } Residue;

struct Chain { char ch_name;
               int beg, end;
               int model;
               char seg_name[5];
             } Chain;


struct Segment { char seg_name[5];
               int beg, end;
               int model;
               }  Segment;

struct Model { int model_n;
               int beg, end;
               }  Model;


void read_PDB_atoms(FILE *fp1, int *n_atoms, struct Atom **atoms, struct Trj *trj);
//void read_PDB_atoms(FILE *fp1, int *n_atoms, struct Atom *atoms);

void read_atom_pdb(char *buf, struct Atom *atom);

void hpsort(struct Atom *ra, int n);
void hpsort_trj(struct Atom *ra, int n, struct Trj *trj);
void hpsort_trj_order(struct Atom *ra, int n, struct Trj *trj);

void make_system(struct System *system, struct Trj *trj);

void copy_atom(struct Atom *atom_to, struct Atom atom_from);
int cmp_atoms(const void *p1, const void *p2);

FILE *file_open(char *fname, char *acc);

#define MAX_DIST_CA 4.3
void eliminate_alt_loc(struct System *system, struct Trj *trj);
double calc_order(double **coords1, double **coords2, int n_models);

/**************************************************
*  Each program has its options in flag_par 
*  (flags or values for parameters)
*
*  init_flag_par initializes flag_par
*  check_cmd_line reads options from command line
*  print_info_flag_par prints options       
*
***************************************************/
struct Flag_par {  
	       char file_in_pdb[120];
	       char file_out[120];
               char at1[5];
               char at2[5];
               int rat1;
               int rat2;
               int verbose;
	     } Flag_par;

void init_flag_par(struct Flag_par *flag_par);
void check_cmd_line(int argc, char *argv[], struct Flag_par *flag_par);
void print_info_flag_par(struct Flag_par flag_par);
double calc_order(double **coords1, double **coords2, int n_models);

/**** MAIN *****/

int main(int argc, char *argv[]) {
      FILE *fp1;
      char buf[1024], fn[120];
      struct Flag_par flag_par;
      struct System system,ref;
      struct Trj trj, trj_ref;
      double **coords1, **coords2, S2;
      int found_1, found_2;
      int i,j,k;

        /* Initialize options and parameters for this program */
        init_flag_par(&flag_par);      
        /* Read options from comman line */
        check_cmd_line(argc, argv, &flag_par);      
        /* Print information on options and parameters */
        print_info_flag_par(flag_par);      
	fp1 = fopen(flag_par.file_in_pdb,"r");
        read_PDB_atoms(fp1, &(system.n_atoms), &(system.atoms), &trj);
        fclose(fp1);

        /* create the residue, chain, segment structures */
        make_system(&system, &trj);  
printf("N. residui: %i -- N. atomi: %i\n", system.n_residues, system.n_atoms);

fp1 = file_open(flag_par.file_out,"w");
coords1 = calloc(trj.nf,sizeof(double *));
coords2 = calloc(trj.nf,sizeof(double *));
//printf("%s\n", flag_par.at1);
//printf("%s\n", flag_par.at2);
//exit(0);
for(j=0; j<system.n_residues; j++)
if(
   (flag_par.rat1 == 0 && flag_par.rat2 == 0) ||
   (flag_par.rat1 == -1 && flag_par.rat2 == 0 && j > 0) ||
   (flag_par.rat1 == 0 && flag_par.rat2 == 1 && j < system.n_residues - 1) 
  )
{
for(k=system.residues[j + flag_par.rat1].beg,found_1=0; k<=system.residues[j + flag_par.rat1].end; k++)
{
if( !strcmp(system.atoms[k].at_name,flag_par.at1)) 
{
found_1=1;
for(i=0; i< trj.nf; i++)
coords1[i] = trj.coor[k + i * trj.naxf];
}
}
for(k=system.residues[j + flag_par.rat2].beg,found_2=0; k<=system.residues[j + flag_par.rat2].end; k++)
{
if((!strcmp(system.atoms[k].at_name,flag_par.at2))) 
{
found_2=1;
for(i=0; i< trj.nf; i++)
coords2[i] = trj.coor[k + i * trj.naxf];
}
}
if(found_1 && found_2) 
{
S2 = calc_order(coords1,coords2, trj.nf);
fprintf(fp1,"%s %i %c %lf\n", system.residues[j].res_name, system.residues[j].res_n, system.residues[j].chain, S2);
}
}
fclose(fp1);
}
void check_cmd_line(int argc, char *argv[], struct Flag_par *flag_par)
{
	int i;
        char tmp[100];
        char extension[100];

	if(argc < 5) 
	{
	printf("Usage:\n"); 
	printf("./suppos file_pdb_in [-]atom_name_1 [+]atom_name_2 file_out  [Options]\n"); 
	printf("Options:\n"); 
	printf("-v (verbose mode)"); 
	printf("\n"); 
	exit(1);
	}
        
        strcpy((*flag_par).file_in_pdb, argv[1]);
        strcpy((*flag_par).at1, argv[2]);
        strcpy((*flag_par).at2, argv[3]);
        strcpy((*flag_par).file_out, argv[4]);

        if((*flag_par).at1[0] == '-') 
        {
         (*flag_par).rat1 = -1;
        for(i = 1; i< strlen((*flag_par).at1); i++)
         (*flag_par).at1[i-1] = (*flag_par).at1[i];
         (*flag_par).at1[i-1] = '\0';  
        }
        if((*flag_par).at2[0] == '+') 
        { 
         (*flag_par).rat2 = 1;
        for(i = 1; i< strlen((*flag_par).at2); i++)
         (*flag_par).at2[i-1] = (*flag_par).at2[i];
         (*flag_par).at2[i-1] = '\0';  
        }
printf("\'-\' before the name of the first atom refers to the atom in the preceding residue\n");
printf("\'+\' before the name of the second atom refers to the atom in the following residue\n");
if(
   ((*flag_par).rat1 == -1 && (*flag_par).rat2 == 0 ) ||
   ((*flag_par).rat1 == 0 && (*flag_par).rat2 == 1 ) ||
   ((*flag_par).rat1 == 0 && (*flag_par).rat2 == 0 )
  ) ; 
else
{
printf("Please, use either \'-\' or nothing before the name of the first atom\n either \'+\' or nothing before the name of the second atom\n...exiting...\n");
exit(0);
}
	for (i = 5; i < argc; i++) 
        {
		if (!strncmp(argv[i],"-v",3)) (*flag_par).verbose = 1;
		else 
                {
                 printf("I don't know option %s\n", argv[i]);
                 exit(2);
                }
        }
}

void init_flag_par(struct Flag_par *flag_par)
{
strcpy( (*flag_par).file_out,"");
(*flag_par).rat1=0;
(*flag_par).rat2=0;
(*flag_par).verbose=0;
}

void print_info_flag_par(struct Flag_par flag_par)
{
        printf("file pdb in: %s\n", flag_par.file_in_pdb);
        printf("file out: %s\n", flag_par.file_out);
        printf("\n");
}

double calc_order(double **coords1, double **coords2, int n_models)
{
int i;
double x,y,z,x2,y2,z2,xy,xz,yz,S2,t;
x=y=z=x2=y2=z2=xy=xz=yz=0;
for(i = 0; i< n_models; i++)
{
x=coords2[i][0] - coords1[i][0];
y=coords2[i][1] - coords1[i][1];
z=coords2[i][2] - coords1[i][2];
t=sqrt(x*x + y*y + z*z);
x = x/t;
y = y/t;
z = z/t;
x2 = x2 + x*x;
y2 = y2 + y*y;
z2 = z2 + z*z;
xy = xy + x*y;
xz = xz + x*z;
yz = yz + y*z;
}
x2 = x2 / (double) n_models;
y2 = y2 / (double) n_models;
z2 = z2 / (double) n_models;
xy = xy / (double) n_models;
xz = xz / (double) n_models;
yz = yz / (double) n_models;
S2 = (3.0/2.0)*(x2*x2 + y2*y2 + z2*z2 + 2 * xy *xy + 2 * xz * xz + 2* yz * yz) - 1.0/2.0;
return S2;
}

/* Legge i vari modelli dal file pdb e memoriza il numero di atomi*/
/* Questa funzione alloca la memoria necessaria e legge atomi e traiettoria */
void read_PDB_atoms(FILE *fp1, int *n_atoms, struct Atom *(*atoms), struct Trj *trj)
{
	char buf[120];
	int i=0, k, n_models,is_trj, n = 2;
        int mod_id=0;
        struct Atom tmp_atom;

// first check if it is a trajectory or a single structure
	while(fgets(buf,120,fp1) != NULL )
        {
    	if((!strncmp("ATOM",buf,4) || !strncmp("HETATM",buf,6)) && mod_id==0) i++;
	      if(!strncmp("ENDMDL",buf,6)) 
              mod_id++; 
        }
*n_atoms = i; 
n_models = mod_id;
if(mod_id <= 1) mod_id = 1;
  (*atoms) = calloc(*n_atoms , sizeof(struct Atom));
  if((*atoms) == NULL) 
    {
     printf("could not allocate memory for %i atoms... exiting...\n", *n_atoms);
     exit(0);
    }
  (*trj).nf = mod_id; 
  (*trj).naxf = *n_atoms; 
  trj->coor = calloc((*trj).nf * (*trj).naxf, sizeof(double *));
  if(trj->coor == NULL) 
    {
     printf("could not allocate memory for %i atoms... exiting...\n", (*trj).nf * (*trj).naxf);
     exit(0);
    }
   for(i = 0; i < (*trj).nf * (*trj).naxf; i++)
    {
    trj->coor[i] =  calloc(3, sizeof(double));
    if(trj->coor[i] == NULL) 
    {
     printf("could not allocate memory for %i-th atom coordinates... exiting...\n", i);
     exit(0);
    }
    }
        rewind(fp1);
        i = 0;
	while(fgets(buf,120,fp1) != NULL)
    	if(!strncmp("ATOM",buf,4) || !strncmp("HETATM",buf,6)) 
            {
                        if(i < *n_atoms)
                        {
			read_atom_pdb(buf, &((*atoms)[i]));
                        for(k=0;k<3;k++)
                        (*trj).coor[i][k] = (*atoms)[i].coor[k];
                        i++;
                        }
                        else if(mod_id > 1)
                        {
			read_atom_pdb(buf, &tmp_atom);
                        for(k=0;k<3;k++)
                        (*trj).coor[i][k] = tmp_atom.coor[k];
                        i++;
                        }
//            if(!(i%n)) {printf("%i atoms read\n",i); n = n*2;}
	    }
	    else
	      if(!strncmp("ENDMDL",buf,6)) 
              mod_id++; 
            for(i = 0; i < *n_atoms; i++) (*atoms)[i].i = i;
}

/* Legge un le informazioni di un atomo da una stringa buf e le salva 
nella strutture atom.*/
void read_atom_pdb(char *buf, struct Atom *atom)
{

    char at_rec[5];
    char tok[10];              
 
    strncpy(tok,buf,6);
    tok[6] = '\0';
    sscanf(tok,"%s", at_rec);

    if((strncmp("ATOM",at_rec,4)) * (strncmp("HETATM",at_rec,6)))
    {
    	printf("The ATOM line does not start with string ATOM... exiting...\n");
        exit(1);
    }

    strncpy(tok,buf + 6,5);
    tok[5] = '\0';
    sscanf(tok,"%i",&(atom->at_n));

    strncpy(tok,buf + 12,4);
    tok[4] = '\0';
    sscanf(tok,"%s", atom->at_name);
 
    strncpy(tok,buf + 16,1);
    tok[1] = '\0';
    if(sscanf(tok,"%c", &(atom->alt_loc)) == -1) atom->alt_loc=' ';
/*    else if ((atom->alt_loc=='A') || (atom->alt_loc=='1')) 
    atom->alt_loc=' '; */ 

	strncpy(tok,buf + 17,3);
    tok[3] = '\0'; 
    sscanf(tok,"%s", atom->res_name);

	strncpy(tok,buf + 21,1);
    tok[1] = '\0';
    if(sscanf(tok,"%c", &(atom->chain)) == EOF) atom->chain = ' ';

    strncpy(tok,buf + 22,4);
    tok[4] = '\0';
    sscanf(tok,"%i", &(atom->res_n));

	strncpy(tok,buf + 26,1);
    tok[1] = '\0';
    if (sscanf(tok,"%c", &(atom->res_ins)) == EOF) atom->res_ins=' ';

    strncpy(tok,buf + 30,8);
    tok[8] = '\0';
    sscanf(tok,"%lf", &(atom->coor[0]));

	strncpy(tok,buf + 38,8);
    tok[8] = '\0';
    sscanf(tok,"%lf", &(atom->coor[1]));

    strncpy(tok,buf + 46,8);
    tok[8] = '\0';
    sscanf(tok,"%lf", &(atom->coor[2]));
//    (*atom).coor = (*atom).ref;
    strncpy(tok,buf + 54,6);
    tok[6] = '\0';
    sscanf(tok,"%lf", &(atom->occ));

	strncpy(tok,buf + 60,6);
    tok[6] = '\0';
    sscanf(tok,"%lf", &(atom->temp));
/* Commentate per evitare problemi con i segmenti */
    if(strlen(buf) > 76)
    { 
      	strncpy(tok,buf + 72,4);
       	tok[4] = '\0';
		if (sscanf(tok,"%s", (atom->segid)) == EOF) 
			strcpy(atom->segid,"    ");
    }
    else strcpy(atom->segid,"    ");


    if(strlen(buf) > 78)
	{
       	strncpy(tok,buf + 76,2);
       	tok[2] = '\0';
       	if (sscanf(tok,"%s", (atom->element)) == EOF) 
			strcpy(atom->element,"UN");
	}

    if(strlen(buf) > 80)
	{
       	strncpy(tok,buf + 78,2);
       	tok[2] = '\0';
       	if (sscanf(tok,"%s", (atom->pdb_chrg)) == EOF) 
			strcpy(atom->pdb_chrg,"  ");
	}

}

void copy_atom(struct Atom *atom_to, struct Atom atom_from)
{
int i;
strcpy((*atom_to).at_name,atom_from.at_name);
(*atom_to).alt_loc = atom_from.alt_loc;
strcpy((*atom_to).res_name,atom_from.res_name);
(*atom_to).chain = atom_from.chain;
strcpy((*atom_to).element,atom_from.element);
(*atom_to).model = atom_from.model;
(*atom_to).at_n = atom_from.at_n;
(*atom_to).res_n = atom_from.res_n;
(*atom_to).res_ins = atom_from.res_ins;

for(i=0;i<3;i++)
{
//printf("copio %8.3lf\n", atom_from.coor[i]);
(*atom_to).coor[i] = atom_from.coor[i];
//printf("copiato %8.3lf\n", (*atom_to).coor[i]);
}
(*atom_to).occ = atom_from.occ;
(*atom_to).temp = atom_from.temp;
strcpy((*atom_to).segid,atom_from.segid);
strcpy((*atom_to).pdb_chrg,atom_from.pdb_chrg);
(*atom_to).mass = atom_from.mass;
strcpy((*atom_to).chem_group,atom_from.chem_group);
strcpy((*atom_to).group,atom_from.group);
(*atom_to).radius = atom_from.radius;
(*atom_to).charge = atom_from.charge;
(*atom_to).atom_type = atom_from.atom_type;
(*atom_to).i = atom_from.i;
}

/* funzione per controllo esistenza e apertura di un file*/
FILE *file_open(char *fname,char *acc) {
    FILE *fp;
    fp =fopen(fname,acc);
    if (fp==NULL)
        {
        fprintf(stderr,"unable to open file %s\n",fname);
        exit(1);
    }
    return(fp);
}


/* creazione ed inizializzazione strutture di sistema*/
void make_system(struct System *system, struct Trj *trj)
{
	int n_atoms = system->n_atoms;
	struct Atom *atoms = system->atoms;
	struct Residue *residues = system->residues; 
	int *p_n_residues = &(system->n_residues); 
	int *p_n_chains = &(system->n_chains); 
	struct Chain *chains = system->chains; 
	int *p_n_segments  = &(system->n_segments); 
	struct Segment *segments = system->segments;
	int *p_n_models = &(system->n_models); 
	struct Model *models = system->models; 

	int imodel=-1, isegment=-1, ichain=-1, iresidue=-1;
	char res_ins = '*';
	char chain = '*';
	char segid[5] = "*****";
	int model = -1;

	int i;
	int res_n = -MAX_N_RESIDUES;

//        if( system->n_atoms < 100000)
//        {
//        printf("sorting %i atoms using quicksort\n", system->n_atoms);  
//	qsort(system->atoms, system->n_atoms, sizeof(struct Atom), &cmp_atoms);
//        }
//        else
//      {

        printf("sorting %i atoms using heapsort\n", system->n_atoms);  
//        for(i=0;i<system->n_atoms;i++) system->atoms[i].i = i;
	hpsort_trj(system->atoms, system->n_atoms,trj);
//	hpsort(system->atoms, system->n_atoms);
//        }

        eliminate_alt_loc(system, trj);
        printf("after eliminate_alt_loc %i atoms left\n", system->n_atoms);  
	n_atoms = system->n_atoms;
// qui conta models, segments, chains, residue,
        model = -1;
	for (i=0; i<n_atoms; i++)
		if(atoms[i].model != model)
		{
//                        printf("atom %i: model: %i vs. %i\n", i, atoms[i].model, model);
			imodel++;
			isegment++;
			ichain++;
			iresidue++;
                        model=atoms[i].model;
                        strcpy(segid,atoms[i].segid);
                        chain=atoms[i].chain;
                        res_n=atoms[i].res_n;
                        res_ins=atoms[i].res_ins;
                }
		else if(strcmp(atoms[i].segid,segid))
			{
 //                       printf("atom %i: segid: %s vs. %s\n", i, atoms[i].segid, segid);
			isegment++;
			ichain++;
			iresidue++;
                        strcpy(segid,atoms[i].segid);
                        chain=atoms[i].chain;
                        res_n=atoms[i].res_n;
                        res_ins=atoms[i].res_ins;
                        }
		else if(atoms[i].chain != chain) 
                        {
			ichain++;
			iresidue++;
                        chain=atoms[i].chain;
                        res_n=atoms[i].res_n;
                        res_ins=atoms[i].res_ins;
                        }
		else if(atoms[i].res_n != res_n) 
                       {
			iresidue++;
                        res_n=atoms[i].res_n;
                        res_ins=atoms[i].res_ins;
                       }
		else if(atoms[i].res_ins != res_ins) 
                       {
			iresidue++;
                        res_ins=atoms[i].res_ins;
                       }
    imodel++;
    isegment++;
    ichain++;
    iresidue++;
    system->models = calloc(imodel, sizeof(Model));
    if(system->models==NULL) {printf("Could not allocate memory for System.models... Exiting...\n"); exit(0);}
    system->segments = calloc(isegment, sizeof(Segment));
    if(system->segments==NULL) {printf("Could not allocate memory for System.segments... Exiting...\n"); exit(0);}
    system->chains = calloc(ichain, sizeof(Chain));
    if(system->chains==NULL) {printf("Could not allocate memory for System.chains... Exiting...\n"); exit(0);}
    system->residues = calloc(iresidue, sizeof(Residue));
    if(system->residues==NULL) {printf("Could not allocate memory for System.residues... Exiting...\n"); exit(0);}
    system->valid_res = calloc(iresidue, sizeof(Valid_res));
    if(system->valid_res==NULL) {printf("Could not allocate memory for System.valid_res... Exiting...\n"); exit(0);}

//        printf("%i %i %i %i\n", imodel, isegment, ichain, iresidue);
	imodel=-1, isegment=-1, ichain=-1, iresidue=-1;
        model = -1;
	for (i=0; i<n_atoms; i++)
	{
		if(atoms[i].model != model)
		{
			imodel++;
			isegment++;
			ichain++;
			iresidue++;
//        printf("%i %i %i %i\n", imodel, isegment, ichain, iresidue);
 
			strcpy(system->residues[iresidue].res_name, atoms[i].res_name);
			system->residues[iresidue].res_n = atoms[i].res_n;
			system->residues[iresidue].res_ins = atoms[i].res_ins;
			system->residues[iresidue].chain = atoms[i].chain;
			system->residues[iresidue].model = atoms[i].model;
			strcpy(system->residues[iresidue].seg_name, atoms[i].segid);
			system->residues[iresidue].beg = i;
			if (iresidue != 0) system->residues[iresidue - 1].end = i-1;
        	chain = atoms[i].chain; 
        	res_n = atoms[i].res_n; 
        	res_ins = atoms[i].res_ins; 
  
/*		if (sscanf("%s", atoms[i].segid) == 0) 
            	strcpy(segid,"     "); 
        	else */
            	strcpy(segid, atoms[i].segid); 
			model=atoms[i].model;        

			system->chains[ichain].ch_name = chain;
			system->chains[ichain].beg = i;
                 	system->chains[ichain].model=atoms[i].model;
			strcpy(system->chains[ichain].seg_name, atoms[i].segid);
			if (ichain != 0) system->chains[ichain - 1].end = i-1;

			strcpy(system->segments[isegment].seg_name, segid);
        	system->segments[isegment].model=atoms[i].model;
        	system->segments[isegment].beg = i;
 			if (isegment != 0) system->segments[isegment - 1].end = i-1;

        	system->models[imodel].beg=i;
        	system->models[imodel].model_n=model;
        	if (imodel != 0) system->models[imodel - 1].end =i-1;
        //printf("Sono qui: model: %i prev: %i\n", atoms[i].model,model);
        //printf("1: %i %s\n", i, (*system).atoms[i].at_name);
		}
		else if(strcmp(atoms[i].segid,segid))
			{
				isegment++;
  				ichain++;
  				iresidue++;

				strcpy(system->residues[iresidue].res_name, atoms[i].res_name);
				system->residues[iresidue].res_n = atoms[i].res_n;
				system->residues[iresidue].res_ins = atoms[i].res_ins;
				system->residues[iresidue].chain = atoms[i].chain;
				system->residues[iresidue].model = atoms[i].model;
				strcpy(system->residues[iresidue].seg_name, atoms[i].segid);
				system->residues[iresidue].beg = i;
			/*	residues[iresidue].res_type = restyp(atoms[i].res_name); */
				if (iresidue != 0) system->residues[iresidue - 1].end = i-1;
			    chain = atoms[i].chain; 
    		    res_n = atoms[i].res_n; 
       		 	res_ins = atoms[i].res_ins; 
/*       	 		if (sscanf("%s", atoms[i].segid) == 0) 
           		 	strcpy(segid,"     "); 
        		else */
            		strcpy(segid, atoms[i].segid); 
        
				system->chains[ichain].ch_name = chain;
				system->chains[ichain].beg = i;
        		system->chains[ichain].model=atoms[i].model;
				strcpy(system->chains[ichain].seg_name, atoms[i].segid);
				if (ichain != 0) system->chains[ichain - 1].end = i-1;

	      		strcpy(system->segments[isegment].seg_name, segid);
    	    	system->segments[isegment].model=atoms[i].model;
       	 		system->segments[isegment].beg = i;
 				if (isegment != 0) system->segments[isegment - 1].end = i-1;
        //printf("Sono qui: segid: %s prev: %s\n", atoms[i].segid,segid);
        //printf("2: %i %s\n", i, (*system).atoms[i].at_name);

			}
			else if(atoms[i].chain != chain) 
				{
					ichain++;
					iresidue++;

					strcpy(system->residues[iresidue].res_name, atoms[i].res_name);
					system->residues[iresidue].res_n = atoms[i].res_n;
					system->residues[iresidue].res_ins = atoms[i].res_ins;
					system->residues[iresidue].chain = atoms[i].chain;
					system->residues[iresidue].model = atoms[i].model;
					strcpy(system->residues[iresidue].seg_name, atoms[i].segid);
					system->residues[iresidue].beg = i;
				/*	residues[iresidue].res_type = restyp(atoms[i].res_name); */
					if (iresidue != 0) system->residues[iresidue - 1].end = i-1;
        			chain = atoms[i].chain; 
       		 		res_n = atoms[i].res_n; 
        			res_ins = atoms[i].res_ins; 
					system->chains[ichain].ch_name = chain;
					system->chains[ichain].beg = i;
        			system->chains[ichain].model=atoms[i].model;
					strcpy(system->chains[ichain].seg_name, atoms[i].segid);
					if (ichain != 0) system->chains[ichain - 1].end = i-1;
        //printf("Sono qui: chain: %c prev: %c\n", atoms[i].chain,chain);
        //printf("3: %i %s\n", i, (*system).atoms[i].at_name);
				}
				else if(atoms[i].res_n != res_n) 
					{
					iresidue++;

					strcpy(system->residues[iresidue].res_name, atoms[i].res_name);
					res_n = atoms[i].res_n;
					res_ins = atoms[i].res_ins;
					system->residues[iresidue].res_n = res_n;
					system->residues[iresidue].res_ins = res_ins;
					system->residues[iresidue].chain = chain;
					system->residues[iresidue].model = atoms[i].model;
					strcpy(system->residues[iresidue].seg_name, atoms[i].segid);
					system->residues[iresidue].beg = i;
				/*	residues[iresidue].res_type = restyp(atoms[i].res_name); */
					if (iresidue != 0) system->residues[iresidue - 1].end = i-1;
        //printf("Sono qui: res_n: %i prev: %i\n", atoms[i].res_n,res_n);
        //printf("5: %i %s\n", i, (*system).atoms[i].at_name);
					}
					else if(atoms[i].res_ins != res_ins) 
					{
					iresidue++;
					strcpy(system->residues[iresidue].res_name, atoms[i].res_name);
					res_n = atoms[i].res_n;
					res_ins = atoms[i].res_ins;
					system->residues[iresidue].res_n = res_n;
					system->residues[iresidue].res_ins = res_ins;
					system->residues[iresidue].chain = chain;
					system->residues[iresidue].model = atoms[i].model;
					strcpy(system->residues[iresidue].seg_name, atoms[i].segid);
					system->residues[iresidue].beg = i;
					if (iresidue != 0) system->residues[iresidue - 1].end = i-1;
					}
	}
if(n_atoms != 0)
    {
    system->segments[isegment].end = i-1;
	system->chains[ichain].end = i-1;
	system->residues[iresidue].end = i-1;
    system->models[imodel].end = i-1;
    }
	*p_n_chains = ichain+1;
	*p_n_segments = isegment+1;
	*p_n_residues = iresidue+1;
	*p_n_models = imodel+1;

	printf("########################################################\n"); 
	printf("# SYSTEM:                                              #\n"); 
	printf("########################################################\n\n"); 
	printf("atoms = %8i, residues = %8i, chains = %8i, segments = %8i, models = %8i\n\n", n_atoms, *p_n_residues, *p_n_chains,*p_n_segments,*p_n_models  );
	printf("########################################################\n\n");
}
//void hpsort(struct Atom *ra, int n, int *cmp_atoms)
void hpsort(struct Atom *ra, int n)
{  
    int N, i, parent, child;  
    struct Atom rra;  
    N = n;
    i = n/2;
    for (;;) { 
        if (i > 0) { 
            i--;           
            copy_atom(&(rra),(ra[i]));
//             t = arr[i];    
        } else {     
            n--;           
            if (n == 0) return; 
            copy_atom(&(rra),(ra[n]));
//            t = arr[n];    
            copy_atom(&(ra[n]),(ra[0]));
//            arr[n] = arr[0]; 
        }  
  
        parent = i; 
        child = i*2 + 1; 
  
        while (child < n) {  
//            if (child + 1 < n  &&  arr[child + 1] > arr[child]) {  
            if (child + 1 < n  && (cmp_atoms(&(ra[child + 1]),&(ra[child])) > 0)) {
                child++; 
            }  
            if (cmp_atoms(&(ra[child]),&(rra)) > 0) {   
                copy_atom(&(ra[parent]),(ra[child]));
//                arr[parent] = arr[child]; 
                parent = child;  
                //child = parent*2-1; 
                child = parent*2+1;
            } else {  
                break; 
            }  
        }  
       copy_atom(&(ra[parent]),(rra));
//   arr[parent] = t;   
    }  
}
void hpsort_trj(struct Atom *ra, int n, struct Trj *trj) 
{  
    int N, i, parent, child, j, k;  
    struct Atom rra;  
    double **tmptrj;

    tmptrj = calloc(n,sizeof(double *));
    for(i = 0; i < n; i++) 
      tmptrj[i] = calloc(3, sizeof(double));

// per la traiettoria do ad ogni atomo un indice, ordino prima gli atomi
// e poi con gli indici ordino la traiettoria
//    for(i = 0; i < n; i++) ra[i].i = i;
   N = n;
    i = n/2;
    for (;;) { 
        if (i > 0) { 
            i--;           
            copy_atom(&(rra),(ra[i]));
//             t = arr[i];    
        } else {     
            n--;           
            if (n == 0) goto final_ops;
            copy_atom(&(rra),(ra[n]));
//            t = arr[n];    
            copy_atom(&(ra[n]),(ra[0]));
//            arr[n] = arr[0]; 
        }  
  
        parent = i; 
        child = i*2 + 1; 
  
        while (child < n) {  
//            if (child + 1 < n  &&  arr[child + 1] > arr[child]) {  
            if (child + 1 < n  && (cmp_atoms(&(ra[child + 1]),&(ra[child])) > 0)) {
                child++; 
            }  
            if (cmp_atoms(&(ra[child]),&(rra)) > 0) {   
                copy_atom(&(ra[parent]),(ra[child]));
//                arr[parent] = arr[child]; 
                parent = child;  
                //child = parent*2-1; 
                child = parent*2+1;
            } else {  
                break; 
            }  
        }  
       copy_atom(&(ra[parent]),(rra));
//   arr[parent] = t;   
    }  
    //qui la traiettoria
final_ops:
       n = (*trj).naxf;
//       for(j = 0; j < n; j++)
//        printf("%i %i\n", j, ra[j].i);
    for(i = 0; i < (*trj).nf; i++)   
       {
       for(j = 0; j < n; j++)
         for(k = 0; k < 3; k++)
            tmptrj[j][k] = (*trj).coor[i * (*trj).naxf + ra[j].i][k];
       for(j = 0; j < n; j++)
         for(k = 0; k < 3; k++)
            (*trj).coor[i * (*trj).naxf + j][k] =  tmptrj[j][k]; 
       }
}


void hpsort_trj_order(struct Atom *ra, int n, struct Trj *trj) 
{  
    int N, i, parent, child, j, k, *order;  
    struct Atom rra;  
    double **tmptrj;

    tmptrj = calloc(n ,sizeof(double *));
    for(i = 0; i < n; i++) 
      tmptrj[i] = calloc(3, sizeof(double));

    order = calloc(n,sizeof(int));
//printf("sono qui 2\n");
// per la traiettoria do ad ogni atomo un indice, ordino prima gli atomi
// e poi con gli indici ordino la traiettoria
//    for(i = 0; i < n; i++) ra[i].i = i;
       for(j = 0; j < n; j++)
        order[j] = ra[j].i;
//       for(j = 0; j < n; j++)
//        printf("%i %i\n", j, ra[j].i);
   N = n;
    i = n/2;
    for (;;) { 
        if (i > 0) { 
            i--;           
            copy_atom(&(rra),(ra[i]));
//             t = arr[i];    
        } else {     
            n--;           
            if (n == 0) goto final_ops;
            copy_atom(&(rra),(ra[n]));
//            t = arr[n];    
            copy_atom(&(ra[n]),(ra[0]));
//            arr[n] = arr[0]; 
        }  
  
        parent = i; 
        child = i*2 + 1; 
  
        while (child < n) {  
//            if (child + 1 < n  &&  arr[child + 1] > arr[child]) {  
            if (child + 1 < n  && (ra[child + 1].i > ra[child].i)) {
                child++; 
            }  
            if (ra[child].i > rra.i) {   
                copy_atom(&(ra[parent]),(ra[child]));
//                arr[parent] = arr[child]; 
                parent = child;  
                //child = parent*2-1; 
                child = parent*2+1;
            } else {  
                break; 
            }  
        }  
       copy_atom(&(ra[parent]),(rra));
//   arr[parent] = t;   
    }  
final_ops:
       n = (*trj).naxf;
    //qui la traiettoria
//       for(j = 0; j < n; j++)
//        printf("%i %i\n", j, ra[j].i);
    for(i = 0; i < (*trj).nf; i++)   
       {
       for(j = 0; j < n; j++)
         for(k = 0; k < 3; k++)
            tmptrj[order[j]][k] = (*trj).coor[i * (*trj).naxf + j][k];
       for(j = 0; j < n; j++)
         for(k = 0; k < 3; k++)
            (*trj).coor[i * (*trj).naxf + j][k] =  tmptrj[j][k]; 
       }
}

#define MAX_DIST_CA 4.3
void eliminate_alt_loc(struct System *system, struct Trj *trj)
{
int i,j,k,l;
char buf[256];
j = 1;
for (i=1; i< (*system).n_atoms; i++)   
{
/*
write_atom_pdb(buf,(*system).atoms[i]);
printf("%s\n",buf);
*/
if(
strcmp((*system).atoms[i-1].at_name,(*system).atoms[i].at_name) ||
strcmp((*system).atoms[i-1].res_name,(*system).atoms[i].res_name) ||
strcmp((*system).atoms[i-1].segid,(*system).atoms[i].segid) ||
((*system).atoms[i-1].chain != (*system).atoms[i].chain) ||
((*system).atoms[i-1].model != (*system).atoms[i].model) ||
((*system).atoms[i-1].res_n != (*system).atoms[i].res_n) ||
((*system).atoms[i-1].res_ins != (*system).atoms[i].res_ins)
)
{
//printf("sono qui %i/%i\n",i,(*system).n_atoms);
if(i != j)
{
copy_atom(&((*system).atoms[j]), ((*system).atoms[i]));
for(k = 0; k < (*trj).nf; k++)
for(l = 0; l<3; l++)
(*trj).coor[k * (*trj).naxf + j][l] = (*trj).coor[k * (*trj).naxf + i][l];
}
/* toglie la alt_loc */
(*system).atoms[j].alt_loc = ' ';
j++;
//if((i+1) == (*system).n_atoms-1) 
//copy_atom(&((*system).atoms[j]), ((*system).atoms[i+1]));
}
else 
{
// Commentate per evitare eccessivo output
//printf("Elimino l'atomo:\n");
//write_atom_pdb(buf,(*system).atoms[i]);
//printf("%s\n",buf);
//printf("perche' c'e' gia':\n");
//write_atom_pdb(buf,(*system).atoms[i-1]);
//printf("%s\n",buf);
}

}
if ((*system).n_atoms > 0)
{
(*system).n_atoms = j;
for(k = 0; k < (*trj).nf; k++)
for(j = 0; j< (*system).n_atoms; j++)
for(l = 0; l<3; l++)
(*trj).coor[k * (*system).n_atoms + j][l] = (*trj).coor[k * (*trj).naxf + j][l];
(*trj).naxf = (*system).n_atoms;
}
}

int cmp_atoms(const void *p1, const void *p2)
{
	struct Atom A_atom, B_atom;
	int check = 0 ; 

	A_atom = *((struct Atom *)p1); 
	B_atom = *((struct Atom *)p2); 


/**** confronto    if (numbers[min] <= numbers[mid]) *****/

	if( A_atom.model < B_atom.model) check = -1; 
	else if (A_atom.model == B_atom.model) 
    {
        check = 0;
    	if( strcmp(A_atom.segid,B_atom.segid) < 0) check = -1; 
    	else if (!strcmp(A_atom.segid,B_atom.segid)) 
        {
        	check = 0;
           	if( (int) A_atom.chain < (int) B_atom.chain ) check = -1;
           	else if( A_atom.chain ==  B_atom.chain )
            {
            	check = 0;
                if( A_atom.res_n <  B_atom.res_n ) check = -1;
                else if(A_atom.res_n ==  B_atom.res_n )
                {
                	check = 0;
                    if((int) A_atom.res_ins < (int) B_atom.res_ins ) check = -1;
                    else if( A_atom.res_ins ==  B_atom.res_ins )
                    {
                    	check = 0;
                        if( strcmp(A_atom.at_name,B_atom.at_name) < 0) check = -1; 
                        else if (!strcmp(A_atom.at_name,B_atom.at_name)) 
                        {
                        check = 0;
                        if( (int) A_atom.alt_loc < (int) B_atom.alt_loc ) check = -1;
                        else if( A_atom.alt_loc ==  B_atom.alt_loc ) check = 0;
                        else check = 1;
			}
                        else check = 1;
                    }
                        else check = 1;
                }
                        else check = 1;
            }
                        else check = 1;
        }
                        else check = 1;
	}
        else check = 1;
	return check;

}

