#include "include.h"
#include "include_dm.h"
#include "omega.h"
#include "IDetection.h"

extern int helicit_(int*, int*, int*, int*);
extern void init_kin_(double*, double*, double*, double*, double*, double*, double*, double*);
extern void squared_(char*,int*,char*,int*,double[],int*,int*);

#if RELIC==2
extern void squaredn_(char*,int*,char*,int*,double[],int*,int*);
#endif

/*--------------------------------------------------------------------*/

void defparticle(struct particle* p, char* name,  double mass, int type)
/* defines the mass, name and type of the particle p */
{
	p->name=name;
	p->mass=mass;
	p->type=type;
	return;
}

/*--------------------------------------------------------------------*/

void init_MSSMlist(struct parameters* param, struct particlelist* MSSMlist)
/* list of MSSM particles involved in indirect detection */
{
	int hscalar=0x001;
	int hfermion=0x10A;
	int hphoton=0x20A;
	int hgluon=hphoton;
	int hvector=0x20E;
	int hgravitino=0x31B;
	int hgraviton=0x41B;
	int htensor=0x41F;
	
	double relicmass=fabs(param->mass_neut[1]);
	double sqrtS=2.*relicmass;
	ModelReIni(sqrtS, param);
	double alphas=alphas_running(sqrtS,param->mass_top_pole,param->mass_b,param);
	double a=alphas/pi;
	double mtmt=mt_mt(param);

	double Mc=running_mass(param->mass_c,param->mass_c,2.*relicmass,param->mass_top_pole,param->mass_b_pole,param);
	double Mb=running_mass(param->mass_b,param->mass_b,2.*relicmass,param->mass_top_pole,param->mass_b_pole,param);
	int nf=5;
	Mb*=sqrt(1.+5.67*a+(35.94-1.36*nf)*a*a+(164.14-nf*(25.77-0.259*nf))*a*a*a);
	Mb=min(Mb,mb_pole(param));
	double Mt=running_mass(mtmt,mtmt,sqrtS,param->mass_top_pole,param->mass_b,param);
	nf=6;
	Mt*=sqrt(1.+5.67*a+(35.94-1.36*nf)*a*a+(164.14-nf*(25.77-0.259*nf))*a*a*a);
	Mt=min(Mt,param->mass_top_pole);
  
	MSSMlist->n=35;
	defparticle(&MSSMlist->list[0], "o1",  fabs(param->mass_neut[1]), hfermion);
	defparticle(&MSSMlist->list[1], "e",  param->mass_e, hfermion);
	defparticle(&MSSMlist->list[2],  "ebar", param->mass_e, hfermion);
	defparticle(&MSSMlist->list[3], "ne", param->mass_nue, hfermion);
	defparticle(&MSSMlist->list[4], "nebar",  param->mass_nue, hfermion);
	defparticle(&MSSMlist->list[5], "m", param->mass_mu, hfermion);
	defparticle(&MSSMlist->list[6], "mbar", param->mass_mu, hfermion);
	defparticle(&MSSMlist->list[7], "nm", param->mass_num, hfermion);
	defparticle(&MSSMlist->list[8], "nmbar", param->mass_num, hfermion);
	defparticle(&MSSMlist->list[9], "l",  param->mass_tau, hfermion);
	defparticle(&MSSMlist->list[10], "lbar", param->mass_tau, hfermion);
	defparticle(&MSSMlist->list[11], "nl",  param->mass_nut, hfermion);
	defparticle(&MSSMlist->list[12],  "nlbar", param->mass_nut, hfermion);
	defparticle(&MSSMlist->list[13], "u",  param->mass_u, hfermion);
	defparticle(&MSSMlist->list[14],  "ubar", param->mass_u, hfermion);
	defparticle(&MSSMlist->list[15], "d",  param->mass_d, hfermion);
	defparticle(&MSSMlist->list[16],  "dbar", param->mass_d, hfermion);
	defparticle(&MSSMlist->list[17], "s",  param->mass_s, hfermion);
	defparticle(&MSSMlist->list[18],  "sbar", param->mass_s, hfermion);
	defparticle(&MSSMlist->list[19], "c", Mc, hfermion);
	defparticle(&MSSMlist->list[20], "cbar", Mc, hfermion);
	defparticle(&MSSMlist->list[21], "b",  Mb, hfermion);
	defparticle(&MSSMlist->list[22],  "bbar", Mb, hfermion);
	defparticle(&MSSMlist->list[23], "t",  Mt, hfermion);
	defparticle(&MSSMlist->list[24],"tbar", Mt, hfermion);
	defparticle(&MSSMlist->list[25], "w", param->mass_W, hvector);
	defparticle(&MSSMlist->list[26],  "wbar", param->mass_W, hvector);
	defparticle(&MSSMlist->list[27], "z",  param->mass_Z, hvector);
	defparticle(&MSSMlist->list[28], "h",  param->mass_h0, hscalar);
	defparticle(&MSSMlist->list[29], "hc",  param->mass_H, hscalar);
	defparticle(&MSSMlist->list[30], "hcbar",  param->mass_H, hscalar);
	defparticle(&MSSMlist->list[31], "hh", param->mass_H0, hscalar);
	defparticle(&MSSMlist->list[32], "h3",  param->mass_A0, hscalar);
	defparticle(&MSSMlist->list[33], "a",  0, hvector);
	defparticle(&MSSMlist->list[34], "g",  0, hvector);
	
	return;
}

/*--------------------------------------------------------------------*/

void init_NMSSMlist(struct parameters* param, struct particlelist* NMSSMlist)
/* list of MSSM particles involved in indirect detection */
{
	int hscalar=0x001;
	int hfermion=0x10A;
	int hphoton=0x20A;
	int hgluon=hphoton;
	int hvector=0x20E;
	int hgravitino=0x31B;
	int hgraviton=0x41B;
	int htensor=0x41F;
	
	double relicmass=fabs(param->mass_neut[1]);
	double sqrtS=2.*relicmass;
	ModelReIni(sqrtS, param);
	double alphas=alphas_running(sqrtS,param->mass_top_pole,param->mass_b,param);
	double a=alphas/pi;

	double mtmt=mt_mt(param);

	double Mc=running_mass(param->mass_c,param->mass_c,2.*relicmass,param->mass_top_pole,param->mass_b_pole,param);
	double Mb=running_mass(param->mass_b,param->mass_b,2.*relicmass,param->mass_top_pole,param->mass_b_pole,param);
	int nf=5;
	Mb*=sqrt(1.+5.67*a+(35.94-1.36*nf)*a*a+(164.14-nf*(25.77-0.259*nf))*a*a*a);
	Mb=min(Mb,mb_pole(param));
	double Mt=running_mass(mtmt,mtmt,sqrtS,param->mass_top_pole,param->mass_b,param);
	nf=6;
	Mt*=sqrt(1.+5.67*a+(35.94-1.36*nf)*a*a+(164.14-nf*(25.77-0.259*nf))*a*a*a);
	Mt=min(Mt,param->mass_top_pole);
	
	NMSSMlist->n=37;
	defparticle(&NMSSMlist->list[0], "o1",  fabs(param->mass_neut[1]), hfermion);
	defparticle(&NMSSMlist->list[1], "e",  param->mass_e, hfermion);
	defparticle(&NMSSMlist->list[2],  "ebar", param->mass_e, hfermion);
	defparticle(&NMSSMlist->list[3], "ne", param->mass_nue, hfermion);
	defparticle(&NMSSMlist->list[4], "nebar",  param->mass_nue, hfermion);
	defparticle(&NMSSMlist->list[5], "m", param->mass_mu, hfermion);
	defparticle(&NMSSMlist->list[6], "mbar", param->mass_mu, hfermion);
	defparticle(&NMSSMlist->list[7], "nm", param->mass_num, hfermion);
	defparticle(&NMSSMlist->list[8], "nmbar", param->mass_num, hfermion);
	defparticle(&NMSSMlist->list[9], "l",  param->mass_tau, hfermion);
	defparticle(&NMSSMlist->list[10], "lbar", param->mass_tau, hfermion);
	defparticle(&NMSSMlist->list[11], "nl",  param->mass_nut, hfermion);
	defparticle(&NMSSMlist->list[12],  "nlbar", param->mass_nut, hfermion);
	defparticle(&NMSSMlist->list[13], "u",  param->mass_u, hfermion);
	defparticle(&NMSSMlist->list[14],  "ubar", param->mass_u, hfermion);
	defparticle(&NMSSMlist->list[15], "d",  param->mass_d, hfermion);
	defparticle(&NMSSMlist->list[16],  "dbar", param->mass_d, hfermion);
	defparticle(&NMSSMlist->list[17], "s",  param->mass_s, hfermion);
	defparticle(&NMSSMlist->list[18],  "sbar", param->mass_s, hfermion);
	defparticle(&NMSSMlist->list[19], "c", Mc, hfermion);
	defparticle(&NMSSMlist->list[20], "cbar", Mc, hfermion);
	defparticle(&NMSSMlist->list[21], "b",  Mb, hfermion);
	defparticle(&NMSSMlist->list[22],  "bbar", Mb, hfermion);
	defparticle(&NMSSMlist->list[23], "t",  Mt, hfermion);
	defparticle(&NMSSMlist->list[24],"tbar", Mt, hfermion);
	defparticle(&NMSSMlist->list[25], "w", param->mass_W, hvector);
	defparticle(&NMSSMlist->list[26],  "wbar", param->mass_W, hvector);
	defparticle(&NMSSMlist->list[27], "z",  param->mass_Z, hvector);
	defparticle(&NMSSMlist->list[28], "h1",  param->mass_h0, hscalar);
	defparticle(&NMSSMlist->list[29], "hc",  param->mass_H, hscalar);
	defparticle(&NMSSMlist->list[30], "hcbar",  param->mass_H, hscalar);
	defparticle(&NMSSMlist->list[31], "h2", param->mass_H0, hscalar);
	defparticle(&NMSSMlist->list[32], "h3",  param->mass_H03, hscalar);
	defparticle(&NMSSMlist->list[33], "ha",  param->mass_A0, hscalar);
	defparticle(&NMSSMlist->list[34], "hb",  param->mass_A02, hscalar);
	defparticle(&NMSSMlist->list[35], "a",  0, hvector);
	defparticle(&NMSSMlist->list[36], "g",  0, hvector);
	
	return;
}

/*--------------------------------------------------------------------*/

void defprocess(struct process *p, char *name, struct particle p1, struct particle p2, struct particle p3, struct particle p4)
/* defines a process p with its name and particles */
{
	p->name=name;
	p->particle1=p1;
	p->particle1=p2;
	p->particle1=p3;
	p->particle1=p4;
	
	return;
}

/*--------------------------------------------------------------------*/

int find_particle(struct particlelist* list, char* name, struct particle* p)
/* finds a particle p from its name in a given particle list */
{
	int npl=list->n;
	int ie=0;
	int found=0;
	while(ie<npl&&!found)
	{
			if(strcmp(name,list->list[ie].name)==0) found=1; 
			else ie++;
	}
	if(found) *p=list->list[ie];
	return found;
}

/*--------------------------------------------------------------------*/

void find_process(char* name, struct particlelist* pl, struct particle* relic, struct process* p, int type)
/* defines process p from its name, using a particle list pl and the dark matter particle relic */
{
	p->particle1=*relic;
	p->particle2=*relic;
	p->name=name;
		
	int i,j;
	int npl=pl->n;
	int np,anp;
	char* relicname=relic->name;
	int nrelic=strlen(relicname);
	int nname=strlen(name);
	int nname1;
	int ntype;
	if(type==0) ntype=2*nrelic; //annihilation
	else if (type==1) ntype=nrelic; //decay
	nname1=nname-ntype; 
	
	int found3,found4;
	found3=found4=0;
	
	char name1[nname1+1];
	for(i=0;i<nname1;i++) name1[i]=name[ntype+i];
	name1[nname1]='\0';
		
	i=0;
	
	while(!(found3&&found4)&&i<nname1-1)
	{
		char beg[i+2];	
		for(j=0;j<=i;j++) beg[j]=name1[j];
		beg[i+1]='\0';

		char end[nname1-i];	
		for(j=0;j<=nname1-i-1;j++) end[j]=name1[j+i+1];
		end[nname1-i-1]='\0';
				
		found3=find_particle(pl,beg,&(p->particle3));
		if(found3) found4=find_particle(pl,end,&(p->particle4));

		i++;
	}
		 
	if(found3&&found4)
	{
		if(!strcmp(p->particle3.name,p->particle4.name)) p->sf34=2; else p->sf34=1;
		int hel1=p->particle1.type;
		int hel3=p->particle3.type;
		int hel4=p->particle4.type;
		p->hel=helicit_(&hel1, &hel1, &hel3, &hel4);
	}
		
	return;
}

int init_IDparameters(struct parameters* param, struct IDparameters* idparam)
/* Initializes the structure idparam for indirect detection calculations using the model point defined in structure param */
{
	int scalar;
	double relicmass=findrelicmass(param, &scalar);
	if(relicmass!=fabs(param->mass_neut[1]))
	{
		printf("Warning: the neutralino1 is not the LSP. Cannot compute direct detection constraints.\n");
		return 0;
	}
	else if (search_mass(relicmass)<0)
	{
		printf("Neutralino mass out of range. Mass must be between 5 and 100000 GeV\n");
		return 0;
	}
	else
	{
		idparam->relicmass=fabs(param->mass_neut[1]);
		double maxenergy=1000;
		ModelIni(param, idparam->relicmass, maxenergy);
		if(param->mass_neut[5]==0) idparam->NMSSM=0;
		else idparam->NMSSM=1;

		if(!init_processes(idparam, param)) return 0;
		
		init_spectrum(idparam);
		
		return 1;
	}
}

/*--------------------------------------------------------------------*/

double IDsigmav(struct IDparameters* idparam)
/* calculates the velocity averaged annihilation cross section for every process in the list of processes proc */
{ 
	int nproc=idparam->o1o1.nproc;
	int flags=1;
	double sqrtS;
	double costh=0;
	double result_tmp[2];
	double sigma=0;
	char* name;
	int len, lenlib, sf34, hel;
	double mass1, mass3, mass4, p12, p34;
	int i;
	char* namelib="o1o1";
    
	for(i=0;i<nproc;i++)
	{ 
		name=idparam->o1o1.list[i].name;
		mass1=idparam->o1o1.list[i].particle1.mass;
		sqrtS=2.*mass1;
		mass3=idparam->o1o1.list[i].particle3.mass;
		mass4=idparam->o1o1.list[i].particle4.mass;
		hel=idparam->o1o1.list[i].hel;
		sf34=idparam->o1o1.list[i].sf34;
		len=strlen(idparam->o1o1.list[i].name);
		lenlib=strlen(namelib);
		
		if(sqrtS<mass3+mass4) idparam->o1o1.list[i].sigmav=0;
		else
		{			  
			init_kin_(&mass1,&mass1,&mass3,&mass4,&sqrtS,&costh,&p12,&p34);
		
			if(idparam->NMSSM==0) squared_(name,&len,namelib,&lenlib,result_tmp, &hel, &flags);
#if RELIC==1
			else
			{
				printf("Please set RELIC to 2 in Makefile and recompile the code, the NMSSM is currently disabled for the indirect detection!\n");
				return 0.;
			}
#endif
#if RELIC==2
			else squaredn_(name,&len,namelib,&lenlib,result_tmp, &hel, &flags);
#endif
			idparam->o1o1.list[i].sigmav=result_tmp[0]/(64.*pi*mass1*mass1)*sqrt(1-(mass3*mass3+mass4*mass4)/(2.*mass1*mass1) + (pow(mass3*mass3-mass4*mass4,2.))/(16*pow(mass1,4.)))*3.8937966e+8*2.997900e-26/2./sf34;
			sigma+=idparam->o1o1.list[i].sigmav;
		}
	}
	idparam->o1o1.dof=0;
	for(i=0;i<nproc;i++) if(idparam->o1o1.list[i].sigmav>0.01*sigma) idparam->o1o1.dof++;
	idparam->o1o1.sigmav=sigma;
	return sigma;
}

/*--------------------------------------------------------------------*/

int init_processes(struct IDparameters* idparam, struct parameters* param)
/* initializes neutralino annihilation processes and Higgs decay */
{
	if(param->width_h0*param->width_H0*param->width_A0*param->width_H==0.) return 0;
	
	double Brh0[11]={param->BRh0mumu,param->BRh0tautau,param->BRh0ss,param->BRh0cc, param->BRh0bb,param->BRh0tt,param->BRh0WW,param->BRh0ZZ, param->BRh0gaZ,param->BRh0gg,param->BRh0gaga};
 	double BrH0[13]={param->BRH0mumu,param->BRH0tautau,param->BRH0ss,param->BRH0cc, param->BRH0bb,param->BRH0tt,param->BRH0WW,param->BRH0ZZ, param->BRH0gaZ,param->BRH0gg,param->BRH0gaga, param->BRH0hZ, param->BRH0hh};
 	double BrA0[13]={param->BRA0mumu,param->BRA0tautau,param->BRA0ss,param->BRA0cc, param->BRA0bb,param->BRA0tt,param->BRA0WW,param->BRA0ZZ, param->BRA0gaZ,param->BRA0gg,param->BRA0gaga, param->BRA0hZ, param->BRA0hh};
 	double BrH03[13]={param->BRH03mumu,param->BRH03tautau,param->BRH03ss,param->BRH03cc, param->BRH03bb,param->BRH03tt,param->BRH03WW,param->BRH03ZZ, param->BRH03gaZ,param->BRH03gg,param->BRH03gaga, param->BRH03hZ, param->BRH03hh};
 	double BrA02[13]={param->BRA02mumu,param->BRA02tautau,param->BRA02ss,param->BRA02cc, param->BRA02bb,param->BRA02tt,param->BRA02WW,param->BRA02ZZ, param->BRA02gaZ,param->BRA02gg,param->BRA02gaga, param->BRA02hZ, param->BRA02hh};
 	double BrH[9]={param->BRHmunu,param->BRHtaunu,param->BRHub,param->BRHus, param->BRHcs,param->BRHcb,param->BRHtb,param->BRHWh, param->BRHWA};
 	
 	int i;
 	//MSSM
	if(param->mass_neut[5]==0)
	{
	
		struct particlelist MSSMlist;
		init_MSSMlist(param,&MSSMlist);
		
		//o1o1 annihilation
		struct particle relic=MSSMlist.list[0]; //neutralino1
		char* tab[25]={"o1o1ebare","o1o1nebarne","o1o1mbarm","o1o1nmbarnm","o1o1lbarl", "o1o1nlbarnl","o1o1ubaru", "o1o1dbard", "o1o1sbars", "o1o1cbarc", "o1o1bbarb", "o1o1tbart", "o1o1wbarw", "o1o1zz","o1o1hh", "o1o1hz","o1o1hcbarhc", "o1o1hcbarw","o1o1hhhh","o1o1hhh","o1o1hhz","o1o1h3h3","o1o1hh3","o1o1hhh3","o1o1h3z"};
		idparam->o1o1.nproc=25;
		for(i=0;i<25;i++) find_process(tab[i],&MSSMlist,&relic,&(idparam->o1o1.list[i]),0);
		
		
		//h0 decay
		struct particle part=MSSMlist.list[28]; //light higgs
		char* tabh0[11]={"hmbarm","hlbarl", "hsbars", "hcbarc", "hbbarb", "htbart", "hwbarw", "hzz", "haz","hgg", "haa"};
		idparam->h0.nproc=11;
	 	for(i=0;i<11;i++){
	  	find_process(tabh0[i],&MSSMlist,&part,&(idparam->h0.list[i]),1);
	  	idparam->h0.list[i].sigmav=Brh0[i];
/*	  	printf("hey %s 1%s 2%s 3 %s 4 %s\n", idparam->h0.list[i].name,idparam->h0.list[i].particle1.name,idparam->h0.list[i].particle2.name, idparam->h0.list[i].particle3.name, idparam->h0.list[i].particle4.name);*/
	  	}
	  	
	  	//H0 decay
	        part=MSSMlist.list[31]; 
	 	char* tabH0[13]={"hhmbarm","hhlbarl", "hhsbars", "hhcbarc", "hhbbarb", "hhtbart", "hhwbarw", "hhzz", "hhaz","hhgg", "hhaa", "hhhz", "hhhh"};
	 	idparam->H0.nproc=13;
	 	for(i=0;i<13;i++){
	  	find_process(tabH0[i],&MSSMlist,&part,&(idparam->H0.list[i]),1);
	  	idparam->H0.list[i].sigmav=BrH0[i];
	
	  	}
	  	
	  	//A0 decay
	    part=MSSMlist.list[32]; 
	 	char* tabA0[13]={"h3mbarm","h3lbarl", "h3sbars", "h3cbarc", "h3bbarb", "h3tbart", "h3wbarw", "h3zz", "h3az","h3gg", "h3aa", "h3hz", "h3hh"};
	 	idparam->A0.nproc=13;
	 	for(i=0;i<13;i++){
	  	find_process(tabA0[i],&MSSMlist,&part,&(idparam->A0.list[i]),1);
	  	idparam->A0.list[i].sigmav=BrA0[i];
	  	}
	  	
	  	//HC decay
	  	part=MSSMlist.list[29];
	  	char* tabHc[9]={"hcmnm", "hclnl", "hcub", "hcus", "hccs", "hccb", "hctb", "hcwh", "hcwh3"};
	  	idparam->hc.nproc=9;
	 	for(i=0;i<9;i++){
	        find_process(tabHc[i],&MSSMlist,&part,&(idparam->hc.list[i]),1);
	  	idparam->hc.list[i].sigmav=BrH[i];
	  	}
	 
	  }
	//NMSSM
	else
	{
		struct particlelist NMSSMlist;
		init_NMSSMlist(param,&NMSSMlist);
		struct particle relic=NMSSMlist.list[0]; //neutralino1
		char* tab[36]={"o1o1bbarb", "o1o1cbarc","o1o1dbard", "o1o1ebare", "o1o1h1h1", "o1o1h1h2", "o1o1h1h3", "o1o1h1ha", "o1o1h1hb", "o1o1h1z", "o1o1h2h2", "o1o1h2h3", "o1o1h2ha", "o1o1h2hb", "o1o1h2z", "o1o1h3h3", "o1o1h3ha", "o1o1h3hb", "o1o1h3z", "o1o1haha", "o1o1hahb","o1o1haz", "o1o1hbhb", "o1o1hbz", "o1o1hcbarhc", "o1o1hcbarw", "o1o1lbarl", "o1o1mbarm", "o1o1nebarne", "o1o1nlbarnl", "o1o1nmbarnm", "o1o1sbars", "o1o1tbart", "o1o1ubaru", "o1o1wbarw", "o1o1zz"};
		idparam->o1o1.nproc=36;
		for(i=0;i<36;i++) find_process(tab[i], &NMSSMlist,&relic, &(idparam->o1o1.list[i]),0);
		
		
		//h1
		struct particle part=NMSSMlist.list[28]; //light higgs
		char* tabh1[11]={"h1mbarm","h1lbarl", "h1sbars", "h1cbarc", "h1bbarb", "h1tbart", "h1wbarw", "h1zz", "h1az","h1gg", "h1aa"};
		idparam->h0.nproc=11;
	 	for(i=0;i<11;i++)
	 	{
			find_process(tabh1[i],&NMSSMlist,&part,&(idparam->h0.list[i]),1);
			idparam->h0.list[i].sigmav=Brh0[i];
	  	}
	  	
	  	
	  	//h2
		part=NMSSMlist.list[31]; //light higgs
		char* tabh2[13]={"h2mbarm","h2lbarl", "h2sbars", "h2cbarc", "h2bbarb", "h2tbart", "h2wbarw", "h2zz", "h2az","h2gg", "h2aa", "h2h1z", "h2h1h1"};
		idparam->H0.nproc=13;
	 	for(i=0;i<13;i++)
	 	{
			find_process(tabh2[i],&NMSSMlist,&part,&(idparam->H0.list[i]),1);
			idparam->H0.list[i].sigmav=BrH0[i];
	  	}
	  	
	  	//h3
		part=NMSSMlist.list[32]; //light higgs
		char* tabh3[13]={"h3mbarm","h3lbarl", "h3sbars", "h3cbarc", "h3bbarb", "h3tbart", "h3wbarw", "h3zz", "h3az","h3gg", "h3aa", "h3h1z", "h3h1h1"};
		idparam->h3.nproc=13;
	 	for(i=0;i<13;i++)
	 	{
			find_process(tabh3[i],&NMSSMlist,&part,&(idparam->h3.list[i]),1);
			idparam->h3.list[i].sigmav=BrH03[i];
	  	}
	  	
	  	//HC decay
	  	part=NMSSMlist.list[29];
	  	char* tabHc[9]={"hcmnm", "hclnl", "hcub", "hcus", "hccs", "hccb", "hctb", "hcwh1", "hcwha"};
	  	idparam->hc.nproc=9;
	 	for(i=0;i<9;i++)
	 	{
	        find_process(tabHc[i],&NMSSMlist,&part,&(idparam->hc.list[i]),1);
			idparam->hc.list[i].sigmav=BrH[i];
	  	}
	  	
	  	//A0 decay
	    part=NMSSMlist.list[33]; 
	 	char* tabA0[13]={"hambarm","halbarl", "hasbars", "hacbarc", "habbarb", "hatbart", "hawbarw", "hazz", "haaz","hagg", "haaa", "hah1z", "hah1h1"};
	 	idparam->A0.nproc=13;
	 	for(i=0;i<13;i++)
	 	{
			find_process(tabA0[i],&NMSSMlist,&part,&(idparam->A0.list[i]),1);
			idparam->A0.list[i].sigmav=BrA0[i];
	  	}
	  	
	  	//A02 decay
	    part=NMSSMlist.list[34]; 
	 	char* tabA02[13]={"hbmbarm","hblbarl", "hbsbars", "hbcbarc", "hbbbarb", "hbtbart", "hbwbarw", "hbzz", "hbaz","hbgg", "hbaa", "hbh1z", "hbh1h1"};
	 	idparam->hb.nproc=13;
	 	for(i=0;i<13;i++)
	 	{
			find_process(tabA02[i],&NMSSMlist,&part,&(idparam->hb.list[i]),1);
			idparam->hb.list[i].sigmav=BrA02[i];
	  	}
	}

	IDsigmav(idparam);
	
	return 1;
}
	
void copyprocesses(struct processes* proc1, struct processes* proc2)
/* copies the list of processes proc2 in proc1 */
{
	proc1->nproc=proc2->nproc;
	proc1->sigmav=proc2->sigmav;
	proc1->dof=proc2->dof;
	 
	int i;
	for(i=0;i<proc2->nproc;i++)
	{
		proc1->list[i].name=proc2->list[i].name;
		proc1->list[i].sigmav=proc2->list[i].sigmav;
		proc1->list[i].particle1=(proc2->list[i].particle1);
		proc1->list[i].particle2=(proc2->list[i].particle2);
		proc1->list[i].particle3=(proc2->list[i].particle3);
		proc1->list[i].particle4=(proc2->list[i].particle4);
	}

	return;
}

void higgsprocess(struct processes* proc, struct IDparameters* idparam, double sigmav, char* name)
/* initializes list of processes proc with the list of Higgs decays
The corresponding thermally averaged cross-sections are the Higgs branching ratios multipled by the cross section sigmav of the process which lead to the production of the Higgs */
{
 	int ih;
 	
	if(strcmp(name,"h")==0||strcmp(name,"h1")==0)
	{
		copyprocesses(proc,&(idparam->h0));

		for(ih=0;ih<proc->nproc;ih++)
		{
			proc->list[ih].sigmav=idparam->h0.list[ih].sigmav*sigmav;
		}
	}
	else if(strcmp(name,"hh")==0||strcmp(name,"h2")==0)
	{
		copyprocesses(proc,&(idparam->H0));
		for(ih=0;ih<proc->nproc;ih++)
		{
			proc->list[ih].sigmav=idparam->H0.list[ih].sigmav*sigmav;
		}
	}
	else if(strcmp(name,"h3")==0&&idparam->NMSSM)
	{
		copyprocesses(proc,&(idparam->h3));
		for(ih=0;ih<proc->nproc;ih++)
		{
		proc->list[ih].sigmav=idparam->h3.list[ih].sigmav*sigmav;
		}
	}
	else if(strcmp(name,"h3")==0||strcmp(name,"ha")==0)
	{
		copyprocesses(proc,&(idparam->A0));
		for(ih=0;ih<proc->nproc;ih++)
		{
			proc->list[ih].sigmav=idparam->A0.list[ih].sigmav*sigmav;
		}
	}
	else if(strcmp(name,"hb")==0)
	{
		copyprocesses(proc,&(idparam->hb));
		for(ih=0;ih<proc->nproc;ih++)
		{
			proc->list[ih].sigmav=idparam->hb.list[ih].sigmav*sigmav;
		}
	}
	else if(strcmp(name,"hc")==0||strcmp(name,"hcbar")==0)
	{
		copyprocesses(proc,&(idparam->hc));
		for(ih=0;ih<proc->nproc;ih++)
		{
			proc->list[ih].sigmav=idparam->hc.list[ih].sigmav*sigmav;
		}
	}
	else printf(" particle %s is not a Higgs\n", name);

	return;
}

void init_spectrumaux(struct processes* proc, struct IDparameters* idparam, double sqrtS)
/* initializes idparam->spectrum for a list of processes proc with center of mass energy sqrtS */
{
	int i,j,k;
	int channel;
	char names[2][10];
	double mass[2];
	double ene[2];
	int i_ene[2];

	for(i=0;i<proc->nproc;i++)
	{	
		if(proc->list[i].sigmav>=0.01*proc->sigmav)
		{
		
			sprintf(names[0],"%s",proc->list[i].particle3.name);
			sprintf(names[1],"%s",proc->list[i].particle4.name);
			mass[0]=proc->list[i].particle3.mass;
			mass[1]=proc->list[i].particle4.mass;
			ene[0]=0.5*(sqrtS +(mass[0]*mass[0]-mass[1]*mass[1])/sqrtS);
			ene[1]=sqrtS-ene[0];

			for(j=0;j<=1;j++)
			{ 
				//check if energy ene[j] is already present in idparam->spectrum.energy and find the corresponding index k
				//if not, increase the length of idparam->spectrum.energy

				k=0;
				while (ene[j]!=idparam->spectrum.energy.a[k][0]&&k<idparam->spectrum.energy.n1) k++;
				if(k==idparam->spectrum.energy.n1)
				{
					idparam->spectrum.energy.n1++;
					idparam->spectrum.sigmav.n1++;
					idparam->spectrum.energy.a[k][0]=ene[j];
					for(channel=0;channel<number_channels;channel++) idparam->spectrum.sigmav.a[k][channel]=0;
				}
				i_ene[j]=k;
				   
				//fill idparam->spectrum.sigmav according to the particle j   
				if(strcmp(names[j],"e")==0||strcmp(names[j],"ebar")==0)  idparam->spectrum.sigmav.a[i_ene[j]][0]+=0.5*proc->list[i].sigmav;
				else if(strcmp(names[j],"m")==0||strcmp(names[j],"mbar")==0)  idparam->spectrum.sigmav.a[i_ene[j]][1]+=0.5*proc->list[i].sigmav;
				else if(strcmp(names[j],"l")==0||strcmp(names[j],"lbar")==0)  idparam->spectrum.sigmav.a[i_ene[j]][2]+=0.5*proc->list[i].sigmav;
				else if(strcmp(names[j],"u")==0||strcmp(names[j],"ubar")==0)  idparam->spectrum.sigmav.a[i_ene[j]][3]+=0.5*proc->list[i].sigmav;
				else if(strcmp(names[j],"d")==0||strcmp(names[j],"dbar")==0)  idparam->spectrum.sigmav.a[i_ene[j]][3]+=0.5*proc->list[i].sigmav;
				else if(strcmp(names[j],"s")==0||strcmp(names[j],"sbar")==0)  idparam->spectrum.sigmav.a[i_ene[j]][3]+=0.5*proc->list[i].sigmav;
				else if(strcmp(names[j],"c")==0||strcmp(names[j],"cbar")==0)  idparam->spectrum.sigmav.a[i_ene[j]][4]+=0.5*proc->list[i].sigmav;
				else if(strcmp(names[j],"b")==0||strcmp(names[j],"bbar")==0)  idparam->spectrum.sigmav.a[i_ene[j]][5]+=0.5*proc->list[i].sigmav;
				else if(strcmp(names[j],"t")==0||strcmp(names[j],"tbar")==0)  idparam->spectrum.sigmav.a[i_ene[j]][6]+=0.5*proc->list[i].sigmav;
				else if(strcmp(names[j],"w")==0||strcmp(names[j],"wbar")==0)  idparam->spectrum.sigmav.a[i_ene[j]][7]+=0.5*proc->list[i].sigmav;
				else if(strcmp(names[j],"z")==0)  idparam->spectrum.sigmav.a[i_ene[j]][8]+=0.5*proc->list[i].sigmav;
				else if(strcmp(names[j],"g")==0)  idparam->spectrum.sigmav.a[i_ene[j]][9]+=0.5*proc->list[i].sigmav;
				else if(strcmp(names[j],"a")==0)  idparam->spectrum.sigmav.a[i_ene[j]][10]+=0.5*proc->list[i].sigmav;
				else if(strcmp(names[j],"ne")==0||strcmp(names[j],"nebar")==0)  idparam->spectrum.sigmav.a[i_ene[j]][11]+=0.5*proc->list[i].sigmav;
				else if(strcmp(names[j],"nm")==0||strcmp(names[j],"nmbar")==0)  idparam->spectrum.sigmav.a[i_ene[j]][12]+=0.5*proc->list[i].sigmav;
				else if(strcmp(names[j],"nl")==0||strcmp(names[j],"nlbar")==0)  idparam->spectrum.sigmav.a[i_ene[j]][13]+=0.5*proc->list[i].sigmav;
				else
				{ 
					//Recursive call to init_spectrumaux to handle Higgs decays
					struct processes ph;

					higgsprocess(&ph, idparam, 0.5*proc->list[i].sigmav,  names[j]);
					init_spectrumaux(&ph, idparam, ene[j]);
				} 
		
			}
		
		}
	
	}
	return;
}

void init_spectrum( struct IDparameters* idparam)
/* initializes idparam->spectrum
idparam->spectrum.sigmav.a[iene][iparticle] contains the thermally averaged cross-section at which SM particles are produced through DM annihilation. The particles hadronize with energy idparam->spectrum.energy.a[iene][0]
iparticle=0->electron, 1->muon, 2->tau, 3->light quarks, 4-> quark c, 5->quark b, 6->quark top, 7->W, 8->Z, 9->gluon, 10->photon, 11->neutrino e, 12-> neutrino mu, 13-> neutrino tau */
{
	double relicmass=idparam->relicmass;
	int i,j;
	double sqrtS=2.*relicmass;
	for(i=0;i<10;i++) idparam->spectrum.energy.a[i][0]=0;
	for(i=0;i<number_channels;i++)
	{
		idparam->spectrum.sigmav.a[0][i]=0;
	}
	idparam->spectrum.energy.a[0][0]=0.5*sqrtS;
	idparam->spectrum.energy.n1=1;
	idparam->spectrum.energy.n2=1;
	
	idparam->spectrum.sigmav.n1=1;
	idparam->spectrum.sigmav.n2=number_channels;
	idparam->spectrum.relicmass=relicmass;
	
	init_spectrumaux(&(idparam->o1o1),idparam, sqrtS);
	
	return;
}

/*--------------------------------------------------------------------*/

double annihilation_sigmav_calculator(char name[])
/* "container" function scanning the SLHA file "name" and the total annihilation cross section */
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;
	
	if(param.THDM_model>0||param.SM==1) return 0.;

	if(param.width_h0*param.width_H0*param.width_A0*param.width_H==0.) widthcalc(name,&param);
	if(param.width_h0*param.width_H0*param.width_A0*param.width_H==0.) return -1.;

	struct IDparameters idparam;

	if(init_IDparameters(&param,&idparam)) return idparam.o1o1.sigmav;
	else return -1.;
}

/*--------------------------------------------------------------------*/

int indirect_fermi_calculator(char name[], char* sigma)
/* "container" function scanning the SLHA file "name" and checking if the point is excluded by FERMI-LAT */
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return -1;
	
	if(param.THDM_model>0||param.SM==1) return -1;			
	
	if(param.width_h0*param.width_H0*param.width_A0*param.width_H==0.) widthcalc(name,&param);
	if(param.width_h0*param.width_H0*param.width_A0*param.width_H==0.) return -1.;	
	
	struct IDparameters idparam;
	
	if(init_IDparameters(&param,&idparam))
	{
		struct fermi fe;
		init_fermi(&fe);

		double delta_standard=deltalikelihood_fermi(&fe, &idparam, "standard");
		
		double dchi2lim;
		if (idparam.o1o1.dof==1) dchi2lim=4;
		else if (idparam.o1o1.dof==2) dchi2lim=6.18;
		else if (idparam.o1o1.dof==3) dchi2lim=8.02;
		else if (idparam.o1o1.dof==4) dchi2lim=9.72;
		else if (idparam.o1o1.dof==5) dchi2lim=11.31;
		else if (idparam.o1o1.dof==6) dchi2lim=12.85;
		
		int i=find_isigma(sigma);
		return (delta_standard<-pvalue_fromchi2(i,idparam.o1o1.dof));
	}
	else return -1;
}

/*--------------------------------------------------------------------*/

int indirect_ams02_calculator(char name[], char*sigma)
/* "container" function scanning the SLHA file "name" and checking if the point is excluded by AMS-02 */
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return -1;
	
	if(param.THDM_model>0||param.SM==1) return -1;
	
	if(param.width_h0*param.width_H0*param.width_A0*param.width_H==0.) widthcalc(name,&param);
	if(param.width_h0*param.width_H0*param.width_A0*param.width_H==0.) return -1.;
	
	struct IDparameters idparam;

	if(init_IDparameters(&param,&idparam))
	{
		double dchi2_standard=deltalikelihood_AMS(&idparam, "standard");
		
		double dchi2lim;
		
		int i=find_isigma(sigma);
		return (dchi2_standard<-pvalue_fromchi2(i,idparam.o1o1.dof));
	}
	else return -1;

}
