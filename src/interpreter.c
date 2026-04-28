#include "include.h"

/*---------------------------------------------------------------------*/

int check_nameobs(obsname* obs)
/* Checks if the observable name is valid or implemented into SuperIso. Returns 0 if invalid, 1 if valid */
{
	if(
		!(
			!strcmp(obs->type,"BR")||!strcmp(obs->type,"BRuntag")||!strcmp(obs->type,"dGamma/dq2")||!strcmp(obs->type,"R-1")||!strcmp(obs->type,"AI")||!strcmp(obs->type,"ACP")||!strcmp(obs->type,"AFB")||!strcmp(obs->type,"FL")||!strcmp(obs->type,"FT")||!strcmp(obs->type,"ATRe")||!strcmp(obs->type,"AT1")||!strcmp(obs->type,"AT2")||!strcmp(obs->type,"AT3")||!strcmp(obs->type,"AT4")||!strcmp(obs->type,"AT5")||!strcmp(obs->type,"HT1")||!strcmp(obs->type,"HT2")||!strcmp(obs->type,"HT3")||!strcmp(obs->type,"alpha")||!strcmp(obs->type,"AIm")||!strcmp(obs->type,"P1")||!strcmp(obs->type,"P2")||!strcmp(obs->type,"P3")||!strcmp(obs->type,"P6")||!strcmp(obs->type,"P4prime")||!strcmp(obs->type,"P5prime")||!strcmp(obs->type,"P6prime")||!strcmp(obs->type,"P8")||!strcmp(obs->type,"P8prime")||!strcmp(obs->type,"A3")||!strcmp(obs->type,"A4")||!strcmp(obs->type,"A5")||!strcmp(obs->type,"A7")||!strcmp(obs->type,"A8")||!strcmp(obs->type,"A9")||!strcmp(obs->type,"S3")||!strcmp(obs->type,"S4")||!strcmp(obs->type,"S5")||!strcmp(obs->type,"S7")||!strcmp(obs->type,"S8")||!strcmp(obs->type,"S9")||!strcmp(obs->type,"FH")||!strcmp(obs->type,"AlFB")||!strcmp(obs->type,"AhFB")||!strcmp(obs->type,"AlhFB")
			//SN: additional observables
			||!strcmp(obs->type,"ATReCP")||!strcmp(obs->type,"AT2CP")||!strcmp(obs->type,"P1CP")||!strcmp(obs->type,"P2CP")||!strcmp(obs->type,"P3CP")||!strcmp(obs->type,"P4primeCP")||!strcmp(obs->type,"P5primeCP")||!strcmp(obs->type,"P6primeCP")||!strcmp(obs->type,"P8primeCP")||!strcmp(obs->type,"ATIm")||!strcmp(obs->type,"ATImCP")||!strcmp(obs->type,"S6c")||!strcmp(obs->type,"A6s")||!strcmp(obs->type,"AFBfull")||!strcmp(obs->type,"A6c")||!strcmp(obs->type,"AFBCP")||!strcmp(obs->type,"AFBCPfull")||!strcmp(obs->type,"Q8m")||!strcmp(obs->type,"Q8p")||!strcmp(obs->type,"Q9")||!strcmp(obs->type,"S2s")||!strcmp(obs->type,"S1c")||!strcmp(obs->type,"S2c")||!strcmp(obs->type,"AFL")||!strcmp(obs->type,"A1c")||!strcmp(obs->type,"A2s")
			//
			/* Form Factors */
			||!strcmp(obs->type,"V")||!strcmp(obs->type,"A0")||!strcmp(obs->type,"A1")||!strcmp(obs->type,"A12")||!strcmp(obs->type,"T1")||!strcmp(obs->type,"T2")||!strcmp(obs->type,"T23")||!strcmp(obs->type,"A2")||!strcmp(obs->type,"fp")||!strcmp(obs->type,"f0")||!strcmp(obs->type,"fT")
			/* Result from Hadronic Fit */
			||!strcmp(obs->type,"hplusRe")||!strcmp(obs->type,"hplusIm")||!strcmp(obs->type,"hminusRe")||!strcmp(obs->type,"hminusIm")||!strcmp(obs->type,"hzeroRe")||!strcmp(obs->type,"hzeroIm")
			/*********************************************************/
			||!strcmp(obs->type,"C1")||!strcmp(obs->type,"C2")||!strcmp(obs->type,"C3")||!strcmp(obs->type,"C4")||!strcmp(obs->type,"C5")||!strcmp(obs->type,"C6")||!strcmp(obs->type,"C7")||!strcmp(obs->type,"C8")||!strcmp(obs->type,"C9")||!strcmp(obs->type,"C10")||!strcmp(obs->type,"CQ1")||!strcmp(obs->type,"CQ2")
			||!strcmp(obs->type,"Cprime1")||!strcmp(obs->type,"Cprime2")||!strcmp(obs->type,"Cprime3")||!strcmp(obs->type,"Cprime4")||!strcmp(obs->type,"Cprime5")||!strcmp(obs->type,"Cprime6")||!strcmp(obs->type,"Cprime7")||!strcmp(obs->type,"Cprime8")||!strcmp(obs->type,"Cprime9")||!strcmp(obs->type,"Cprime10")||!strcmp(obs->type,"CprimeQ1")||!strcmp(obs->type,"CprimeQ2")
			||!strcmp(obs->type,"ReC1")||!strcmp(obs->type,"ReC2")||!strcmp(obs->type,"ReC3")||!strcmp(obs->type,"ReC4")||!strcmp(obs->type,"ReC5")||!strcmp(obs->type,"ReC6")||!strcmp(obs->type,"ReC7")||!strcmp(obs->type,"ReC8")||!strcmp(obs->type,"ReC9")||!strcmp(obs->type,"ReC10")||!strcmp(obs->type,"ReCQ1")||!strcmp(obs->type,"ReCQ2")
			||!strcmp(obs->type,"ReCprime1")||!strcmp(obs->type,"ReCprime2")||!strcmp(obs->type,"ReCprime3")||!strcmp(obs->type,"ReCprime4")||!strcmp(obs->type,"ReCprime5")||!strcmp(obs->type,"ReCprime6")||!strcmp(obs->type,"ReCprime7")||!strcmp(obs->type,"ReCprime8")||!strcmp(obs->type,"ReCprime9")||!strcmp(obs->type,"ReCprime10")||!strcmp(obs->type,"ReCprimeQ1")||!strcmp(obs->type,"ReCprimeQ2")
			||!strcmp(obs->type,"ImC1")||!strcmp(obs->type,"ImC2")||!strcmp(obs->type,"ImC3")||!strcmp(obs->type,"ImC4")||!strcmp(obs->type,"ImC5")||!strcmp(obs->type,"ImC6")||!strcmp(obs->type,"ImC7")||!strcmp(obs->type,"ImC8")||!strcmp(obs->type,"ImC9")||!strcmp(obs->type,"ImC10")||!strcmp(obs->type,"ImCQ1")||!strcmp(obs->type,"ImCQ2")
			||!strcmp(obs->type,"ImCprime1")||!strcmp(obs->type,"ImCprime2")||!strcmp(obs->type,"ImCprime3")||!strcmp(obs->type,"ImCprime4")||!strcmp(obs->type,"ImCprime5")||!strcmp(obs->type,"ImCprime6")||!strcmp(obs->type,"ImCprime7")||!strcmp(obs->type,"ImCprime8")||!strcmp(obs->type,"ImCprime9")||!strcmp(obs->type,"ImCprime10")||!strcmp(obs->type,"ImCprimeQ1")||!strcmp(obs->type,"ImCprimeQ2")
			||!strcmp(obs->type,"Vud")||!strcmp(obs->type,"Vus")||!strcmp(obs->type,"Vub")||!strcmp(obs->type,"Vcd")||!strcmp(obs->type,"Vcs")||!strcmp(obs->type,"Vcb")||!strcmp(obs->type,"Vtd")||!strcmp(obs->type,"Vts")||!strcmp(obs->type,"Vtb")
			||!strcmp(obs->type,"ReVud")||!strcmp(obs->type,"ReVus")||!strcmp(obs->type,"ReVub")||!strcmp(obs->type,"ReVcd")||!strcmp(obs->type,"ReVcs")||!strcmp(obs->type,"ReVcb")||!strcmp(obs->type,"ReVtd")||!strcmp(obs->type,"ReVts")||!strcmp(obs->type,"ReVtb")
			||!strcmp(obs->type,"ImVud")||!strcmp(obs->type,"ImVus")||!strcmp(obs->type,"ImVub")||!strcmp(obs->type,"ImVcd")||!strcmp(obs->type,"ImVcs")||!strcmp(obs->type,"ImVcb")||!strcmp(obs->type,"ImVtd")||!strcmp(obs->type,"ImVts")||!strcmp(obs->type,"ImVtb")
		)
	)
	{
		printf("%s: Unknown observable type: %s\n", obs->name, obs->type);
		
		return 0;
	}
	
	if(
		!(
			!strcmp(obs->decay,"BXsgamma")||!strcmp(obs->decay,"BXdgamma")
			||!strcmp(obs->decay,"BXsmumu")||!strcmp(obs->decay,"BXsee")||!strcmp(obs->decay,"BXstautau")||!strcmp(obs->decay,"BXsll")
			||!strcmp(obs->decay,"BKstargamma")||!strcmp(obs->decay,"B0Kstar0gamma")
			||!strcmp(obs->decay,"Bsmumu")||!strcmp(obs->decay,"Bsee")||!strcmp(obs->decay,"Bstautau")||!strcmp(obs->decay,"Bsll")
			||!strcmp(obs->decay,"Bdmumu")||!strcmp(obs->decay,"Bdee")||!strcmp(obs->decay,"Bdtautau")||!strcmp(obs->decay,"Bdll")
			||!strcmp(obs->decay,"BKstarmumu")||!strcmp(obs->decay,"B0Kstar0mumu")
			||!strcmp(obs->decay,"BKstaree")||!strcmp(obs->decay,"B0Kstar0ee")||!strcmp(obs->decay,"BKstartautau")||!strcmp(obs->decay,"B0Kstar0tautau")||!strcmp(obs->decay,"BKstarll")||!strcmp(obs->decay,"B0Kstar0ll")
			||!strcmp(obs->decay,"BKmumu")||!strcmp(obs->decay,"B0K0mumu")
			||!strcmp(obs->decay,"BKee")||!strcmp(obs->decay,"B0K0ee")||!strcmp(obs->decay,"BKtautau")||!strcmp(obs->decay,"B0K0tautau")||!strcmp(obs->decay,"BKll")||!strcmp(obs->decay,"B0K0ll")
			||!strcmp(obs->decay,"Bsphimumu")
			||!strcmp(obs->decay,"Bsphiee")||!strcmp(obs->decay,"Bsphitautau")||!strcmp(obs->decay,"Bsphill")
			||!strcmp(obs->decay,"Bmunu")||!strcmp(obs->decay,"Benu")||!strcmp(obs->decay,"Btaunu")
			||!strcmp(obs->decay,"Dsmunu")||!strcmp(obs->decay,"Dsenu")||!strcmp(obs->decay,"Dstaunu")
			||!strcmp(obs->decay,"Dmunu")||!strcmp(obs->decay,"Denu")
			||!strcmp(obs->decay,"Kmunu/pimunu")
			||!strcmp(obs->decay,"BDmunu")||!strcmp(obs->decay,"BDenu")||!strcmp(obs->decay,"BDtaunu")
			||!strcmp(obs->decay,"BDstarmunu")||!strcmp(obs->decay,"BDstarenu")||!strcmp(obs->decay,"BDstartaunu")
			||!strcmp(obs->decay,"LambdabLambdamumu")||!strcmp(obs->decay,"LambdabLambdaee")||!strcmp(obs->decay,"LambdabLambdatautau")||!strcmp(obs->decay,"LambdabLambdall")
			||!strcmp(obs->decay,"KLpi0nunu")||!strcmp(obs->decay,"Kpinunu")
			||!strcmp(obs->decay,"KLmumu")||!strcmp(obs->decay,"KLee")||!strcmp(obs->decay,"KSmumu")
			||!strcmp(obs->decay,"KLpi0mumu")||!strcmp(obs->decay,"KLpi0ee")
			||!strcmp(obs->decay,"e")||!strcmp(obs->decay,"mu")||!strcmp(obs->decay,"tau") /* Wilson coefficients */
			||!strcmp(obs->decay,"CKM")
		)
	)
	{
		printf("%s: Unknown decay name: %s\n", obs->name, obs->decay);
		
		return 0;
	}
	
	if(
		(
			!strcmp(obs->decay,"BXsmumu")||!strcmp(obs->decay,"BXsee")||!strcmp(obs->decay,"BXstautau")||!strcmp(obs->decay,"BXsll")
			||!strcmp(obs->decay,"BKstarmumu")||!strcmp(obs->decay,"B0Kstar0mumu")
			||!strcmp(obs->decay,"BKstaree")||!strcmp(obs->decay,"B0Kstar0ee")||!strcmp(obs->decay,"BKstartautau")||!strcmp(obs->decay,"B0Kstar0tautau")||!strcmp(obs->decay,"BKstarll")||!strcmp(obs->decay,"B0Kstar0ll")
			||!strcmp(obs->decay,"BKmumu")||!strcmp(obs->decay,"B0K0mumu")
			||!strcmp(obs->decay,"BKee")||!strcmp(obs->decay,"B0K0ee")||!strcmp(obs->decay,"BKtautau")||!strcmp(obs->decay,"B0K0tautau")||!strcmp(obs->decay,"BKll")||!strcmp(obs->decay,"B0K0ll")
			||!strcmp(obs->decay,"Bsphimumu")
			||!strcmp(obs->decay,"Bsphiee")||!strcmp(obs->decay,"Bsphitautau")||!strcmp(obs->decay,"Bsphill")
			||!strcmp(obs->decay,"BDmunu")||!strcmp(obs->decay,"BDenu")||!strcmp(obs->decay,"BDtaunu")
			||!strcmp(obs->decay,"BDstarmunu")||!strcmp(obs->decay,"BDstarenu")||!strcmp(obs->decay,"BDstartaunu")
			||!strcmp(obs->decay,"LambdabLambdamumu")||!strcmp(obs->decay,"LambdabLambdaee")||!strcmp(obs->decay,"LambdabLambdatautau")||!strcmp(obs->decay,"LambdabLambdall")
			||!strcmp(obs->decay,"Dmumu")||!strcmp(obs->decay,"Dee")
		)&&(obs->low<=0.||obs->high<=0.||obs->low>obs->high)
	)
	{
		printf("%s: Wrong binning, please check: low=%g , high=%g\n", obs->name, obs->low, obs->high);
		
		return 0;
	}

	struct parameters param;
	
	Init_param(&param);
	
	slha_adjust(&param);

	param.SM=1;
	
	if(compute_nameobs_ref(obs,&param)<-1.e29)
	{
		printf("%s_%s: Not implemented or wrong name\n", obs->type, obs->decay);
		
		return 0;
	}

	return 1;
}

/*---------------------------------------------------------------------*/

void make_obslist(char names[][50], obsname obs[], int *nobs)
/* Makes a list of observables names "obs" from a string list "names". "nobs" is the original length of the string array */
{
	int ie,ke;
	int je=0;

 	for(ie=0;ie<*nobs;ie++)
	{
		int test=0;

		for(ke=0;ke<ie;ke++)
		{
			if(!strcmp(names[ie],names[ke]))
			{
				printf("Duplicate %s, keeping the first one...\n", names[ie]);

				for(int le=ke;le<*nobs-1;le++) sprintf(names[le],"%s",names[le+1]);
				*nobs=*nobs-1;
				ie--;

				test=1;

				ke=ie;
			}
		}

		if(test==0)
		{
			read_nameobs(names[ie],&obs[je]);
			if(check_nameobs(&obs[je])) je++;
		}
	}

	*nobs=je;

	return;
}

/*---------------------------------------------------------------------*/

int read_obs_list(char filename[], char names[][50])
/* Reads the names of observables from filename into the string array names */ 

{	
	char dummy[50];

	int nbobs=0;
	
	FILE *input;
	
	input=fopen(filename,"r");
	
	while(EOF!=fscanf(input,"%s",dummy))
	{
		if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(input,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
		else
		{
			strcpy(names[nbobs],dummy);
			nbobs++;
		}
	}
	
	fclose(input);
	
	return nbobs;
}

/*---------------------------------------------------------------------*/

double compute_nameobs_ref(obsname* obs, struct parameters* param)
/* Computes the value of observable "obs" using the set of parameters "param" and setting all standard deviation to zero (hence reference value) */
{
	return compute_nameobs(obs,0,param);
}

double compute_nameobs(obsname* obs, int ke, struct parameters* param) /* interpreter */
/* Computes the value of observable "obs" using the set of parameters "param", for the ke-th iteration of the theory error calculation */
{
	double mu_W,mu_b,mu_spec;
	
	if(param->log_mu_W_mass_W>0.) mu_W=pow(2.,param->log_mu_W_mass_W)*param->mass_W; else mu_W=param->mass_W;
	if(param->log_mu_b_mass_b>0.) mu_b=pow(2.,param->log_mu_b_mass_b)*param->mass_b_pole; else mu_b=param->mass_b_pole;
	double lambda_h=0.5;
	if(param->log_mu_spec_lambda_h_mass_b>0.) mu_spec=pow(2.,param->log_mu_spec_lambda_h_mass_b)*sqrt(lambda_h*param->mass_b); else mu_spec=sqrt(lambda_h*param->mass_b);
		
	int ie;
	
	static int keprev=0;
	
	int reset=(ke!=keprev||ke==0);
		
	keprev=ke;
	
	static int WCmu,WCe,WCtau,WCbase2;
	
	if(reset) WCmu=WCe=WCtau=WCbase2=0;

	static double binlowBKmumu,binlowBKee,binlowBKtautau;
	static double binlowBKstarmumu,binlowBKstaree,binlowBKstartautau;
	static double binlowBsphimumu,binlowBsphiee,binlowBsphitautau;
	static double binlowLbLmumu,binlowLbLee,binlowLbLtautau;

	static double binhighBKmumu,binhighBKee,binhighBKtautau;
	static double binhighBKstarmumu,binhighBKstaree,binhighBKstartautau;
	static double binhighBsphimumu,binhighBsphiee,binhighBsphitautau;
	static double binhighLbLmumu,binhighLbLee,binhighLbLtautau;
	
	if(reset) binlowBKmumu=binlowBKee=binlowBKtautau=binlowBKstarmumu=binlowBKstaree=binlowBKstartautau=binlowBsphimumu=binlowBsphiee=binlowBsphitautau=binlowLbLmumu=binlowLbLee=binlowLbLtautau=binhighBKmumu=binhighBKee=binhighBKtautau=binhighBKstarmumu=binhighBKstaree=binhighBKstartautau=binhighBsphimumu=binhighBsphiee=binhighBsphitautau=0.;

	static double obsBKstarmumu[Nobs_BKsll+1],obsBKmumu[Nobs_BKll+1],obsBsphimumu[Nobs_Bsphill+1],obsLbLmumu[Nobs_LbLll+1];
	static double obsBKstaree[Nobs_BKsll+1],obsBKee[Nobs_BKll+1],obsBsphiee[Nobs_Bsphill+1],obsLbLee[Nobs_LbLll+1];
	static double obsBKstartautau[Nobs_BKsll+1],obsBKtautau[Nobs_BKll+1],obsBsphitautau[Nobs_Bsphill+1],obsLbLtautau[Nobs_LbLll+1];


	static double binlowB0K0mumu,binlowB0K0ee,binlowB0K0tautau;
	static double binlowB0Kstar0mumu,binlowB0Kstar0ee,binlowB0Kstar0tautau;

	static double binhighB0K0mumu,binhighB0K0ee,binhighB0K0tautau;
	static double binhighB0Kstar0mumu,binhighB0Kstar0ee,binhighB0Kstar0tautau;
	
	if(reset) binlowB0K0mumu=binlowB0K0ee=binlowB0K0tautau=binlowB0Kstar0mumu=binlowB0Kstar0ee=binlowB0Kstar0tautau=binhighB0K0mumu=binhighB0K0ee=binhighB0K0tautau=binhighB0Kstar0mumu=binhighB0Kstar0ee=binhighB0Kstar0tautau=0.;
	
	static double obsB0Kstar0mumu[Nobs_BKsll+1],obsB0K0mumu[Nobs_BKll+1];
	static double obsB0Kstar0ee[Nobs_BKsll+1],obsB0K0ee[Nobs_BKll+1];
	static double obsB0Kstar0tautau[Nobs_BKsll+1],obsB0K0tautau[Nobs_BKll+1];



	static double complex C0b[11],C0spec[11],C1b[11],C0w[11],C1w[11],C2w[11],C2b[11],Cpb[11],CQpb[3],CQ0b[3],CQ1b[3];
	static double complex C0eb[11],C1eb[11],C0ew[11],C1ew[11],C2ew[11],C2eb[11],Cpeb[11],CQpeb[3],CQ0eb[3],CQ1eb[3];
	static double complex C0taub[11],C1taub[11],C0tauw[11],C1tauw[11],C2tauw[11],C2taub[11],Cptaub[11],CQptaub[3],CQ0taub[3],CQ1taub[3];

	static double complex C0w_base2[11],C1w_base2[11],C2w_base2[11],C0b_base2[11],C0spec_base2[11],C1b_base2[11],C1spec_base2[11],Cpb_base2[11],CQpb_base2[3];

	double alphas_mub;
	double complex Cmu[11],CQmu[3],Ce[11],CQe[3],Ctau[11],CQtau[3];

	if(!strcmp(obs->decay,"BXsgamma")||!strcmp(obs->decay,"BXdgamma")
	||!strcmp(obs->decay,"BXsmumu")||!strcmp(obs->decay,"BXsll")
	||!strcmp(obs->decay,"BKstargamma")||!strcmp(obs->decay,"B0Kstar0gamma")
	||!strcmp(obs->decay,"Bsmumu")
	||!strcmp(obs->decay,"Bdmumu")
	||!strcmp(obs->decay,"BKstarmumu")||!strcmp(obs->decay,"B0Kstar0mumu")
	||!strcmp(obs->decay,"BKstarll")||!strcmp(obs->decay,"B0Kstar0ll")
	||!strcmp(obs->decay,"BKmumu")||!strcmp(obs->decay,"B0K0mumu")
	||!strcmp(obs->decay,"BKll")||!strcmp(obs->decay,"B0K0ll")
	||!strcmp(obs->decay,"Bsphimumu")||!strcmp(obs->decay,"Bsphill")
	||!strcmp(obs->decay,"LambdabLambdamumu")||!strcmp(obs->decay,"LambdabLambdall")
	||!strcmp(obs->decay,"mu") // "mu" stands for WC with mu+mu- in final state
	)
	{
		if(WCmu==0)
		{
			CW_calculator(2,C0w,C1w,C2w,mu_W,param); /* 2 = muon */

			C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,param);
			Cprime_calculator(2,Cpb,CQpb,mu_W,mu_b,param);
			CQ_calculator(2,CQ0b,CQ1b,mu_W,mu_b,param);
			
			/* #pragma... is to suppress erroneous warnings for comparison of array with NULL */
			#pragma GCC diagnostic push
			#pragma GCC diagnostic ignored "-Wtautological-compare"
			if(param->deltaC!=NULL) for(ie=1;ie<=10;ie++) C0b[ie]+=param->deltaC[ie];
			if(param->deltaCp!=NULL) for(ie=1;ie<=10;ie++) Cpb[ie]+=param->deltaCp[ie];
			if(param->deltaCQ!=NULL) for(ie=1;ie<=2;ie++) CQ0b[ie]+=param->deltaCQ[ie];
			if(param->deltaCQp!=NULL) for(ie=1;ie<=2;ie++) CQpb[ie]+=param->deltaCQp[ie];
			#pragma GCC diagnostic pop
			WCmu=1;
		}
	}

	if(!strcmp(obs->decay,"BXsee")||!strcmp(obs->decay,"BXsll")
	||!strcmp(obs->decay,"Bsee")||!strcmp(obs->decay,"Bsll")
	||!strcmp(obs->decay,"Bdee")||!strcmp(obs->decay,"Bdll")
	||!strcmp(obs->decay,"BKstaree")||!strcmp(obs->decay,"B0Kstar0ee")
	||!strcmp(obs->decay,"BKstarll")||!strcmp(obs->decay,"B0Kstar0ll")
	||!strcmp(obs->decay,"BKee")||!strcmp(obs->decay,"B0K0ee")
	||!strcmp(obs->decay,"BKll")||!strcmp(obs->decay,"B0K0ll")
	||!strcmp(obs->decay,"Bsphiee")||!strcmp(obs->decay,"Bsphill")
	||!strcmp(obs->decay,"LambdabLambdaee")||!strcmp(obs->decay,"LambdabLambdall")
	||!strcmp(obs->decay,"e")) // "e" stands for WC with e+e- in final state
	{
		if(WCe==0)
		{
			CW_calculator(1,C0ew,C1ew,C2ew,mu_W,param); /* 1 = electron */

			C_calculator_base1(C0ew,C1ew,C2ew,mu_W,C0eb,C1eb,C2eb,mu_b,param);
			Cprime_calculator(1,Cpeb,CQpeb,mu_W,mu_b,param);
			
			CQ_calculator(1,CQ0eb,CQ1eb,mu_W,mu_b,param);
			
			/* #pragma... is to suppress erroneous warnings for comparison of array with NULL */
			#pragma GCC diagnostic push
			#pragma GCC diagnostic ignored "-Wtautological-compare"
			if(param->deltaC!=NULL) for(ie=1;ie<=10;ie++) C0eb[ie]+=param->deltaC[ie+10*(ie>8)];
			if(param->deltaCp!=NULL) for(ie=1;ie<=10;ie++) Cpeb[ie]+=param->deltaCp[ie+10*(ie>8)];
			if(param->deltaCQ!=NULL) for(ie=1;ie<=2;ie++) CQ0eb[ie]+=param->deltaCQ[ie+2];
			if(param->deltaCQp!=NULL) for(ie=1;ie<=2;ie++) CQpeb[ie]+=param->deltaCQp[ie+2];
			#pragma GCC diagnostic pop

			WCe=1;
		}
	}

	if(!strcmp(obs->decay,"BXstautau")
	||!strcmp(obs->decay,"Bstautau")
	||!strcmp(obs->decay,"Bdtautau")
	||!strcmp(obs->decay,"BKstartautau")||!strcmp(obs->decay,"B0Kstar0tautau")
	||!strcmp(obs->decay,"BKtautau")||!strcmp(obs->decay,"B0K0tautau")
	||!strcmp(obs->decay,"Bsphitautau")
	||!strcmp(obs->decay,"LambdabLambdatautau")
	||!strcmp(obs->decay,"tau")) // "tau" stands for WC with tau+tau- in final state
	{
		if(WCtau==0)
		{
			CW_calculator(3,C0tauw,C1tauw,C2tauw,mu_W,param); /* 3 = tau */

			C_calculator_base1(C0tauw,C1tauw,C2tauw,mu_W,C0taub,C1taub,C2taub,mu_b,param);
			Cprime_calculator(3,Cptaub,CQptaub,mu_W,mu_b,param);
			
			CQ_calculator(3,CQ0taub,CQ1taub,mu_W,mu_b,param);

			/* #pragma... is to suppress erroneous warnings for comparison of array with NULL */
			#pragma GCC diagnostic push
			#pragma GCC diagnostic ignored "-Wtautological-compare"
			if(param->deltaC!=NULL) for(ie=1;ie<=10;ie++) C0taub[ie]+=param->deltaC[ie+20*(ie>8)];
			if(param->deltaCp!=NULL) for(ie=1;ie<=10;ie++) Cptaub[ie]+=param->deltaCp[ie+20*(ie>8)];
			if(param->deltaCQ!=NULL) for(ie=1;ie<=2;ie++) CQ0taub[ie]+=param->deltaCQ[ie+4];
			if(param->deltaCQp!=NULL) for(ie=1;ie<=2;ie++) CQptaub[ie]+=param->deltaCQp[ie+4];
			#pragma GCC diagnostic pop

			WCtau=1;	
		}
	}

	if(!strcmp(obs->type,"AI")&&!strcmp(obs->decay,"BKstargamma"))
	{
		if(WCbase2==0)
		{
			CW_calculator(2,C0w_base2,C1w_base2,C2w_base2,mu_W,param); /* 2 = muon */

			C_calculator_base2(C0w_base2,C1w_base2,mu_W,C0b_base2,C1b_base2,mu_b,param);
			C_calculator_base2(C0w_base2,C1w_base2,mu_W,C0spec_base2,C1spec_base2,mu_spec,param);
			
			Cprime_calculator(2,Cpb_base2,CQpb_base2,mu_W,mu_b,param);		

			/* #pragma... is to suppress erroneous warnings for comparison of array with NULL */
			#pragma GCC diagnostic push
			#pragma GCC diagnostic ignored "-Wtautological-compare"
			if(param->deltaC!=NULL) for(ie=1;ie<=10;ie++) C0b_base2[ie]+=param->deltaC[ie];
			if(param->deltaCp!=NULL) for(ie=1;ie<=10;ie++) Cpb_base2[ie]+=param->deltaCp[ie];
			#pragma GCC diagnostic pop

			WCbase2=1;
		}
		
		double lambda_h=0.5;
		return delta0(C0b_base2,C0spec_base2,C1b_base2,C1spec_base2,Cpb_base2,param,mu_b,mu_spec,lambda_h);
	}

	if(!strcmp(obs->decay,"CKM")) /* CKM parameters */
	{
		if(!strcmp(obs->type,"Vud")||!strcmp(obs->type,"ReVud")) return creal(param->Vud);
		if(!strcmp(obs->type,"Vus")||!strcmp(obs->type,"ReVus")) return creal(param->Vus);
		if(!strcmp(obs->type,"Vub")||!strcmp(obs->type,"ReVub")) return creal(param->Vub);
		if(!strcmp(obs->type,"Vcd")||!strcmp(obs->type,"ReVcd")) return creal(param->Vcd);
		if(!strcmp(obs->type,"Vcs")||!strcmp(obs->type,"ReVcs")) return creal(param->Vcs);
		if(!strcmp(obs->type,"Vcb")||!strcmp(obs->type,"ReVcb")) return creal(param->Vcb);
		if(!strcmp(obs->type,"Vtd")||!strcmp(obs->type,"ReVtd")) return creal(param->Vtd);
		if(!strcmp(obs->type,"Vts")||!strcmp(obs->type,"ReVts")) return creal(param->Vts);
		if(!strcmp(obs->type,"Vtb")||!strcmp(obs->type,"ReVtb")) return creal(param->Vtb);

		if(!strcmp(obs->type,"ImVud")) return cimag(param->Vud);
		if(!strcmp(obs->type,"ImVus")) return cimag(param->Vus);
		if(!strcmp(obs->type,"ImVub")) return cimag(param->Vub);
		if(!strcmp(obs->type,"ImVcd")) return cimag(param->Vcd);
		if(!strcmp(obs->type,"ImVcs")) return cimag(param->Vcs);
		if(!strcmp(obs->type,"ImVcb")) return cimag(param->Vcb);
		if(!strcmp(obs->type,"ImVtd")) return cimag(param->Vtd);
		if(!strcmp(obs->type,"ImVts")) return cimag(param->Vts);
		if(!strcmp(obs->type,"ImVtb")) return cimag(param->Vtb);
	}

	if(!strcmp(obs->decay,"mu")) /* C_mu Wilson coefficients */
	{
		alphas_mub=alphas_running(mu_b,param->mass_top_pole,param->mass_b_pole,param);
		for(ie=1;ie<=10;ie++) Cmu[ie]=C0b[ie]+alphas_mub/4./pi*C1b[ie]+pow(alphas_mub/4./pi,2.)*C2b[ie];
		for(ie=1;ie<=2;ie++) CQmu[ie]=CQ0b[ie]+alphas_mub/4./pi*CQ1b[ie];
				
		if(!strcmp(obs->type,"C1")||!strcmp(obs->type,"ReC1")) return creal(Cmu[1]);
		if(!strcmp(obs->type,"C2")||!strcmp(obs->type,"ReC2")) return creal(Cmu[2]);
		if(!strcmp(obs->type,"C3")||!strcmp(obs->type,"ReC3")) return creal(Cmu[3]);
		if(!strcmp(obs->type,"C4")||!strcmp(obs->type,"ReC4")) return creal(Cmu[4]);
		if(!strcmp(obs->type,"C5")||!strcmp(obs->type,"ReC5")) return creal(Cmu[5]);
		if(!strcmp(obs->type,"C6")||!strcmp(obs->type,"ReC6")) return creal(Cmu[6]);
		if(!strcmp(obs->type,"C7")||!strcmp(obs->type,"ReC7")) return creal(Cmu[7]);
		if(!strcmp(obs->type,"C8")||!strcmp(obs->type,"ReC8")) return creal(Cmu[8]);
		if(!strcmp(obs->type,"C9")||!strcmp(obs->type,"ReC9")) return creal(Cmu[9]);
		if(!strcmp(obs->type,"C10")||!strcmp(obs->type,"ReC10")) return creal(Cmu[10]);

		if(!strcmp(obs->type,"Cprime1")||!strcmp(obs->type,"ReCprime1")) return creal(Cpb[1]);
		if(!strcmp(obs->type,"Cprime2")||!strcmp(obs->type,"ReCprime2")) return creal(Cpb[2]);
		if(!strcmp(obs->type,"Cprime3")||!strcmp(obs->type,"ReCprime3")) return creal(Cpb[3]);
		if(!strcmp(obs->type,"Cprime4")||!strcmp(obs->type,"ReCprime4")) return creal(Cpb[4]);
		if(!strcmp(obs->type,"Cprime5")||!strcmp(obs->type,"ReCprime5")) return creal(Cpb[5]);
		if(!strcmp(obs->type,"Cprime6")||!strcmp(obs->type,"ReCprime6")) return creal(Cpb[6]);
		if(!strcmp(obs->type,"Cprime7")||!strcmp(obs->type,"ReCprime7")) return creal(Cpb[7]);
		if(!strcmp(obs->type,"Cprime8")||!strcmp(obs->type,"ReCprime8")) return creal(Cpb[8]);
		if(!strcmp(obs->type,"Cprime9")||!strcmp(obs->type,"ReCprime9")) return creal(Cpb[9]);
		if(!strcmp(obs->type,"Cprime10")||!strcmp(obs->type,"ReCprime10")) return creal(Cpb[10]);
		
		if(!strcmp(obs->type,"CQ1")||!strcmp(obs->type,"ReCQ1")) return creal(CQmu[1]);
		if(!strcmp(obs->type,"CQ2")||!strcmp(obs->type,"ReCQ2")) return creal(CQmu[2]);

		if(!strcmp(obs->type,"CprimeQ1")||!strcmp(obs->type,"ReCQprime1")) return creal(CQpb[1]);
		if(!strcmp(obs->type,"CprimeQ2")||!strcmp(obs->type,"ReCQprime2")) return creal(CQpb[2]);
		
		if(!strcmp(obs->type,"ImC1")) return cimag(Cmu[1]);
		if(!strcmp(obs->type,"ImC2")) return cimag(Cmu[2]);
		if(!strcmp(obs->type,"ImC3")) return cimag(Cmu[3]);
		if(!strcmp(obs->type,"ImC4")) return cimag(Cmu[4]);
		if(!strcmp(obs->type,"ImC5")) return cimag(Cmu[5]);
		if(!strcmp(obs->type,"ImC6")) return cimag(Cmu[6]);
		if(!strcmp(obs->type,"ImC7")) return cimag(Cmu[7]);
		if(!strcmp(obs->type,"ImC8")) return cimag(Cmu[8]);
		if(!strcmp(obs->type,"ImC9")) return cimag(Cmu[9]);
		if(!strcmp(obs->type,"ImC10")) return cimag(Cmu[10]);

		if(!strcmp(obs->type,"ImCprime1")) return cimag(Cpb[1]);
		if(!strcmp(obs->type,"ImCprime2")) return cimag(Cpb[2]);
		if(!strcmp(obs->type,"ImCprime3")) return cimag(Cpb[3]);
		if(!strcmp(obs->type,"ImCprime4")) return cimag(Cpb[4]);
		if(!strcmp(obs->type,"ImCprime5")) return cimag(Cpb[5]);
		if(!strcmp(obs->type,"ImCprime6")) return cimag(Cpb[6]);
		if(!strcmp(obs->type,"ImCprime7")) return cimag(Cpb[7]);
		if(!strcmp(obs->type,"ImCprime8")) return cimag(Cpb[8]);
		if(!strcmp(obs->type,"ImCprime9")) return cimag(Cpb[9]);
		if(!strcmp(obs->type,"ImCprime10")) return cimag(Cpb[10]);
		
		if(!strcmp(obs->type,"ImCQ1")) return cimag(CQmu[1]);
		if(!strcmp(obs->type,"ImCQ2")) return cimag(CQmu[2]);

		if(!strcmp(obs->type,"ImCprimeQ1")) return cimag(CQpb[1]);
		if(!strcmp(obs->type,"ImCprimeQ2")) return cimag(CQpb[2]);
	}
	
	else if(!strcmp(obs->decay,"e")) /* C_e Wilson coefficients */
	{
		alphas_mub=alphas_running(mu_b,param->mass_top_pole,param->mass_b_pole,param);
		for(ie=1;ie<=10;ie++) Ce[ie]=C0eb[ie]+alphas_mub/4./pi*C1eb[ie]+pow(alphas_mub/4./pi,2.)*C2eb[ie];
		for(ie=1;ie<=2;ie++) CQe[ie]=CQ0eb[ie]+alphas_mub/4./pi*CQ1eb[ie];
		
		if(!strcmp(obs->type,"C1")||!strcmp(obs->type,"ReC1")) return creal(Ce[1]);
		if(!strcmp(obs->type,"C2")||!strcmp(obs->type,"ReC2")) return creal(Ce[2]);
		if(!strcmp(obs->type,"C3")||!strcmp(obs->type,"ReC3")) return creal(Ce[3]);
		if(!strcmp(obs->type,"C4")||!strcmp(obs->type,"ReC4")) return creal(Ce[4]);
		if(!strcmp(obs->type,"C5")||!strcmp(obs->type,"ReC5")) return creal(Ce[5]);
		if(!strcmp(obs->type,"C6")||!strcmp(obs->type,"ReC6")) return creal(Ce[6]);
		if(!strcmp(obs->type,"C7")||!strcmp(obs->type,"ReC7")) return creal(Ce[7]);
		if(!strcmp(obs->type,"C8")||!strcmp(obs->type,"ReC8")) return creal(Ce[8]);
		if(!strcmp(obs->type,"C9")||!strcmp(obs->type,"ReC9")) return creal(Ce[9]);
		if(!strcmp(obs->type,"C10")||!strcmp(obs->type,"ReC10")) return creal(Ce[10]);

		if(!strcmp(obs->type,"Cprime1")||!strcmp(obs->type,"ReCprime1")) return creal(Cpeb[1]);
		if(!strcmp(obs->type,"Cprime2")||!strcmp(obs->type,"ReCprime2")) return creal(Cpeb[2]);
		if(!strcmp(obs->type,"Cprime3")||!strcmp(obs->type,"ReCprime3")) return creal(Cpeb[3]);
		if(!strcmp(obs->type,"Cprime4")||!strcmp(obs->type,"ReCprime4")) return creal(Cpeb[4]);
		if(!strcmp(obs->type,"Cprime5")||!strcmp(obs->type,"ReCprime5")) return creal(Cpeb[5]);
		if(!strcmp(obs->type,"Cprime6")||!strcmp(obs->type,"ReCprime6")) return creal(Cpeb[6]);
		if(!strcmp(obs->type,"Cprime7")||!strcmp(obs->type,"ReCprime7")) return creal(Cpeb[7]);
		if(!strcmp(obs->type,"Cprime8")||!strcmp(obs->type,"ReCprime8")) return creal(Cpeb[8]);
		if(!strcmp(obs->type,"Cprime9")||!strcmp(obs->type,"ReCprime9")) return creal(Cpeb[9]);
		if(!strcmp(obs->type,"Cprime10")||!strcmp(obs->type,"ReCprime10")) return creal(Cpeb[10]);
		
		if(!strcmp(obs->type,"CQ1")||!strcmp(obs->type,"ReCQ1")) return creal(CQe[1]);
		if(!strcmp(obs->type,"CQ2")||!strcmp(obs->type,"ReCQ2")) return creal(CQe[2]);

		if(!strcmp(obs->type,"CprimeQ1")||!strcmp(obs->type,"ReCprimeQ1")) return creal(CQpeb[1]);
		if(!strcmp(obs->type,"CprimeQ2")||!strcmp(obs->type,"ReCprimeQ2")) return creal(CQpeb[2]);

		if(!strcmp(obs->type,"ImC1")) return cimag(Ce[1]);
		if(!strcmp(obs->type,"ImC2")) return cimag(Ce[2]);
		if(!strcmp(obs->type,"ImC3")) return cimag(Ce[3]);
		if(!strcmp(obs->type,"ImC4")) return cimag(Ce[4]);
		if(!strcmp(obs->type,"ImC5")) return cimag(Ce[5]);
		if(!strcmp(obs->type,"ImC6")) return cimag(Ce[6]);
		if(!strcmp(obs->type,"ImC7")) return cimag(Ce[7]);
		if(!strcmp(obs->type,"ImC8")) return cimag(Ce[8]);
		if(!strcmp(obs->type,"ImC9")) return cimag(Ce[9]);
		if(!strcmp(obs->type,"ImC10")) return cimag(Ce[10]);

		if(!strcmp(obs->type,"ImCprime1")) return cimag(Cpeb[1]);
		if(!strcmp(obs->type,"ImCprime2")) return cimag(Cpeb[2]);
		if(!strcmp(obs->type,"ImCprime3")) return cimag(Cpeb[3]);
		if(!strcmp(obs->type,"ImCprime4")) return cimag(Cpeb[4]);
		if(!strcmp(obs->type,"ImCprime5")) return cimag(Cpeb[5]);
		if(!strcmp(obs->type,"ImCprime6")) return cimag(Cpeb[6]);
		if(!strcmp(obs->type,"ImCprime7")) return cimag(Cpeb[7]);
		if(!strcmp(obs->type,"ImCprime8")) return cimag(Cpeb[8]);
		if(!strcmp(obs->type,"ImCprime9")) return cimag(Cpeb[9]);
		if(!strcmp(obs->type,"ImCprime10")) return cimag(Cpeb[10]);
		
		if(!strcmp(obs->type,"ImCQ1")) return cimag(CQe[1]);
		if(!strcmp(obs->type,"ImCQ2")) return cimag(CQe[2]);

		if(!strcmp(obs->type,"ImCprimeQ1")) return cimag(CQpeb[1]);
		if(!strcmp(obs->type,"ImCprimeQ2")) return cimag(CQpeb[2]);
	}

	else if(!strcmp(obs->decay,"tau")) /* C_tau Wilson coefficients */
	{
		alphas_mub=alphas_running(mu_b,param->mass_top_pole,param->mass_b_pole,param);
		for(ie=1;ie<=10;ie++) Ctau[ie]=C0taub[ie]+alphas_mub/4./pi*C1taub[ie]+pow(alphas_mub/4./pi,2.)*C2taub[ie];
		for(ie=1;ie<=2;ie++) CQtau[ie]=CQ0taub[ie]+alphas_mub/4./pi*CQ1taub[ie];
		
		if(!strcmp(obs->type,"C1")||!strcmp(obs->type,"ReC1")) return creal(Ctau[1]);
		if(!strcmp(obs->type,"C2")||!strcmp(obs->type,"ReC2")) return creal(Ctau[2]);
		if(!strcmp(obs->type,"C3")||!strcmp(obs->type,"ReC3")) return creal(Ctau[3]);
		if(!strcmp(obs->type,"C4")||!strcmp(obs->type,"ReC4")) return creal(Ctau[4]);
		if(!strcmp(obs->type,"C5")||!strcmp(obs->type,"ReC5")) return creal(Ctau[5]);
		if(!strcmp(obs->type,"C6")||!strcmp(obs->type,"ReC6")) return creal(Ctau[6]);
		if(!strcmp(obs->type,"C7")||!strcmp(obs->type,"ReC7")) return creal(Ctau[7]);
		if(!strcmp(obs->type,"C8")||!strcmp(obs->type,"ReC8")) return creal(Ctau[8]);
		if(!strcmp(obs->type,"C9")||!strcmp(obs->type,"ReC9")) return creal(Ctau[9]);
		if(!strcmp(obs->type,"C10")||!strcmp(obs->type,"ReC10")) return creal(Ctau[10]);

		if(!strcmp(obs->type,"Cprime1")||!strcmp(obs->type,"ReCprime1")) return creal(Cptaub[1]);
		if(!strcmp(obs->type,"Cprime2")||!strcmp(obs->type,"ReCprime2")) return creal(Cptaub[2]);
		if(!strcmp(obs->type,"Cprime3")||!strcmp(obs->type,"ReCprime3")) return creal(Cptaub[3]);
		if(!strcmp(obs->type,"Cprime4")||!strcmp(obs->type,"ReCprime4")) return creal(Cptaub[4]);
		if(!strcmp(obs->type,"Cprime5")||!strcmp(obs->type,"ReCprime5")) return creal(Cptaub[5]);
		if(!strcmp(obs->type,"Cprime6")||!strcmp(obs->type,"ReCprime6")) return creal(Cptaub[6]);
		if(!strcmp(obs->type,"Cprime7")||!strcmp(obs->type,"ReCprime7")) return creal(Cptaub[7]);
		if(!strcmp(obs->type,"Cprime8")||!strcmp(obs->type,"ReCprime8")) return creal(Cptaub[8]);
		if(!strcmp(obs->type,"Cprime9")||!strcmp(obs->type,"ReCprime9")) return creal(Cptaub[9]);
		if(!strcmp(obs->type,"Cprime10")||!strcmp(obs->type,"ReCprime10")) return creal(Cptaub[10]);
		
		if(!strcmp(obs->type,"CQ1")||!strcmp(obs->type,"ReCQ1")) return creal(CQtau[1]);
		if(!strcmp(obs->type,"CQ2")||!strcmp(obs->type,"ReCQ2")) return creal(CQtau[2]);

		if(!strcmp(obs->type,"CprimeQ1")||!strcmp(obs->type,"ReCprimeQ1")) return creal(CQptaub[1]);
		if(!strcmp(obs->type,"CprimeQ2")||!strcmp(obs->type,"ReCprimeQ2")) return creal(CQptaub[2]);
		
		if(!strcmp(obs->type,"ImC1")) return cimag(Ctau[1]);
		if(!strcmp(obs->type,"ImC2")) return cimag(Ctau[2]);
		if(!strcmp(obs->type,"ImC3")) return cimag(Ctau[3]);
		if(!strcmp(obs->type,"ImC4")) return cimag(Ctau[4]);
		if(!strcmp(obs->type,"ImC5")) return cimag(Ctau[5]);
		if(!strcmp(obs->type,"ImC6")) return cimag(Ctau[6]);
		if(!strcmp(obs->type,"ImC7")) return cimag(Ctau[7]);
		if(!strcmp(obs->type,"ImC8")) return cimag(Ctau[8]);
		if(!strcmp(obs->type,"ImC9")) return cimag(Ctau[9]);
		if(!strcmp(obs->type,"ImC10")) return cimag(Ctau[10]);

		if(!strcmp(obs->type,"ImCprime1")) return cimag(Cptaub[1]);
		if(!strcmp(obs->type,"ImCprime2")) return cimag(Cptaub[2]);
		if(!strcmp(obs->type,"ImCprime3")) return cimag(Cptaub[3]);
		if(!strcmp(obs->type,"ImCprime4")) return cimag(Cptaub[4]);
		if(!strcmp(obs->type,"ImCprime5")) return cimag(Cptaub[5]);
		if(!strcmp(obs->type,"ImCprime6")) return cimag(Cptaub[6]);
		if(!strcmp(obs->type,"ImCprime7")) return cimag(Cptaub[7]);
		if(!strcmp(obs->type,"ImCprime8")) return cimag(Cptaub[8]);
		if(!strcmp(obs->type,"ImCprime9")) return cimag(Cptaub[9]);
		if(!strcmp(obs->type,"ImCprime10")) return cimag(Cptaub[10]);
		
		if(!strcmp(obs->type,"ImCQ1")) return cimag(CQtau[1]);
		if(!strcmp(obs->type,"ImCQ2")) return cimag(CQtau[2]);

		if(!strcmp(obs->type,"ImCprimeQ1")) return cimag(CQptaub[1]);
		if(!strcmp(obs->type,"ImCprimeQ2")) return cimag(CQptaub[2]);
	}

	else if(!strcmp(obs->type,"BR")&&!strcmp(obs->decay,"BXsgamma")) return bsgamma(C0b,C1b,C2b,Cpb,mu_b,mu_W,param)*(1.+param->bsgamma_rand);

	else if(!strcmp(obs->type,"BR")&&!strcmp(obs->decay,"B0Kstar0gamma")) return BR_BKstargamma(0,C0b,C1b,C2b,Cpb,param,mu_b);

	else if(!strcmp(obs->type,"BR")&&!strcmp(obs->decay,"BKstargamma")) return BR_BKstargamma(1,C0b,C1b,C2b,Cpb,param,mu_b);
	
	else if(!strcmp(obs->decay,"Bsmumu"))
	{
		if(!strcmp(obs->type,"BR")) return Bsmumu(C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
		else if(!strcmp(obs->decay,"Bsmumu")) return Bsmumu_untag(C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
	}

	else if(!strcmp(obs->decay,"Bsee"))
	{
		if(!strcmp(obs->type,"BR")) return Bsll(1,C0eb,C1eb,C2eb,CQ0eb,CQ1eb,Cpeb,CQpeb,param,mu_b);
		else if(!strcmp(obs->type,"BRuntag")) return Bsll_untag(1,C0eb,C1eb,C2eb,CQ0eb,CQ1eb,Cpeb,CQpeb,param,mu_b);
	}

	else if(!strcmp(obs->decay,"Bstautau"))
	{
		if(!strcmp(obs->type,"BR")) return Bsll(3,C0taub,C1taub,C2taub,CQ0taub,CQ1taub,Cptaub,CQptaub,param,mu_b);
		else if(!strcmp(obs->type,"BRuntag")) return Bsll_untag(1,C0taub,C1taub,C2taub,CQ0taub,CQ1taub,Cptaub,CQptaub,param,mu_b);
	}
	
	else if(!strcmp(obs->type,"BR")&&!strcmp(obs->decay,"Bdmumu")) return Bdmumu(C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
	else if(!strcmp(obs->type,"BR")&&!strcmp(obs->decay,"Bdee")) return Bdll(1,C0eb,C1eb,C2eb,CQ0eb,CQ1eb,Cpeb,CQpeb,param,mu_b);
	else if(!strcmp(obs->type,"BR")&&!strcmp(obs->decay,"Bdtautau")) return Bdll(3,C0taub,C1taub,C2taub,CQ0taub,CQ1taub,Cptaub,CQptaub,param,mu_b);

	else if(!strcmp(obs->decay,"BXsmumu"))
	{
		if(!strcmp(obs->type,"BR")) return BRBXsll(2,obs->low,obs->high,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b)*(1.+param->BRBXsmumu_lowq2_rand*(obs->high<=8.)+param->BRBXsmumu_highq2_rand*(obs->low>=12.)+param->BRBXsmumu_full_rand*(obs->low<8.&&obs->high>12.));
		else if(!strcmp(obs->type,"AFB")) return A_BXsll_bin(2,obs->low,obs->high,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
	}

	else if(!strcmp(obs->decay,"BXsee"))
	{
		if(!strcmp(obs->type,"BR")) return BRBXsll(1,obs->low,obs->high,C0eb,C1eb,C2eb,CQ0eb,CQ1eb,Cpeb,CQpeb,param,mu_b)*(1.+param->BRBXsee_lowq2_rand*(obs->high<=8.)+param->BRBXsee_highq2_rand*(obs->low>=12.)+param->BRBXsee_full_rand*(obs->low<8.&&obs->high>12.));
		else if(!strcmp(obs->type,"AFB")) return A_BXsll_bin(1,obs->low,obs->high,C0eb,C1eb,C2eb,CQ0eb,CQ1eb,Cpeb,CQpeb,param,mu_b);
	}

	else if(!strcmp(obs->decay,"BXstautau"))
	{
		if(!strcmp(obs->type,"BR")) return BRBXsll(3,obs->low,obs->high,C0taub,C1taub,C2taub,CQ0taub,CQ1taub,Cptaub,CQptaub,param,mu_b)*(1.+param->BRBXstautau_lowq2_rand*(obs->high<=8.)+param->BRBXstautau_highq2_rand*(obs->low>=12.)+param->BRBXstautau_full_rand*(obs->low<8.&&obs->high>12.));
		else if(!strcmp(obs->type,"AFB")) return A_BXsll_bin(3,obs->low,obs->high,C0taub,C1taub,C2taub,CQ0taub,CQ1taub,Cptaub,CQptaub,param,mu_b);
	}

/* B -> K* l+ l- */
	else if(!strcmp(obs->decay,"B0Kstar0mumu"))
	{
		if(reset||(fabs(binlowB0Kstar0mumu-obs->low)>1.e-5)||(fabs(binhighB0Kstar0mumu-obs->high)>1.e-5))
		{
			obsB0Kstar0mumu[0]=BRBKstarll(2,0,obs->low,obs->high,obsB0Kstar0mumu,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b)/(obs->high-obs->low);
		
			binlowB0Kstar0mumu=obs->low;
			binhighB0Kstar0mumu=obs->high;
		}
		
		if(!strcmp(obs->type,"dGamma/dq2")) return obsB0Kstar0mumu[0];
		else if(!strcmp(obs->type,"BR")) return obsB0Kstar0mumu[0]*(obs->high-obs->low);
		else if(!strcmp(obs->type,"AFB")) return obsB0Kstar0mumu[1];
		else if(!strcmp(obs->type,"FL")) return obsB0Kstar0mumu[2];
		else if(!strcmp(obs->type,"FT")) return obsB0Kstar0mumu[3];
		else if(!strcmp(obs->type,"AT1")) return obsB0Kstar0mumu[4];
		else if(!strcmp(obs->type,"AT2")||!strcmp(obs->type,"P1")) return obsB0Kstar0mumu[5];
		else if(!strcmp(obs->type,"AT3")) return obsB0Kstar0mumu[6];
		else if(!strcmp(obs->type,"AT4")) return obsB0Kstar0mumu[7];
		else if(!strcmp(obs->type,"AT5")) return obsB0Kstar0mumu[8];
		else if(!strcmp(obs->type,"HT1")) return obsB0Kstar0mumu[9];
		else if(!strcmp(obs->type,"HT2")) return obsB0Kstar0mumu[10];
		else if(!strcmp(obs->type,"HT3")) return obsB0Kstar0mumu[11];
		else if(!strcmp(obs->type,"alpha")) return obsB0Kstar0mumu[12];
		else if(!strcmp(obs->type,"AIm")) return obsB0Kstar0mumu[13];
		else if(!strcmp(obs->type,"P2")) return obsB0Kstar0mumu[14];
		else if(!strcmp(obs->type,"P3")) return obsB0Kstar0mumu[15];
		else if(!strcmp(obs->type,"P6")) return obsB0Kstar0mumu[16];
		else if(!strcmp(obs->type,"P4prime")) return obsB0Kstar0mumu[17];
		else if(!strcmp(obs->type,"P5prime")) return obsB0Kstar0mumu[18];
		else if(!strcmp(obs->type,"P6prime")) return obsB0Kstar0mumu[19];
		else if(!strcmp(obs->type,"P8")) return obsB0Kstar0mumu[20];
		else if(!strcmp(obs->type,"P8prime")) return obsB0Kstar0mumu[21];
		else if(!strcmp(obs->type,"A7")) return obsB0Kstar0mumu[22];
		else if(!strcmp(obs->type,"A8")) return obsB0Kstar0mumu[23];
		else if(!strcmp(obs->type,"A9")) return obsB0Kstar0mumu[24];
		else if(!strcmp(obs->type,"S3")) return obsB0Kstar0mumu[25];
		else if(!strcmp(obs->type,"S4")) return obsB0Kstar0mumu[26];
		else if(!strcmp(obs->type,"S5")) return obsB0Kstar0mumu[27];
		else if(!strcmp(obs->type,"S7")) return obsB0Kstar0mumu[28];
		else if(!strcmp(obs->type,"S8")) return obsB0Kstar0mumu[29];
		else if(!strcmp(obs->type,"S9")) return obsB0Kstar0mumu[30];
		else if(!strcmp(obs->type,"A3")) return obsB0Kstar0mumu[31];
		else if(!strcmp(obs->type,"A4")) return obsB0Kstar0mumu[32];
		else if(!strcmp(obs->type,"A5")) return obsB0Kstar0mumu[33];
		else if(!strcmp(obs->type,"ATRe")) return obsB0Kstar0mumu[34];
		else if(!strcmp(obs->type,"ATReCP")) return obsB0Kstar0mumu[35];
		else if(!strcmp(obs->type,"AT2CP")||!strcmp(obs->type,"P1CP")) return obsB0Kstar0mumu[36];
		else if(!strcmp(obs->type,"P2CP")) return obsB0Kstar0mumu[37];
		else if(!strcmp(obs->type,"P3CP")) return obsB0Kstar0mumu[38];
		else if(!strcmp(obs->type,"P4primeCP")) return obsB0Kstar0mumu[39];
		else if(!strcmp(obs->type,"P5primeCP")) return obsB0Kstar0mumu[40];
		else if(!strcmp(obs->type,"P6primeCP")) return obsB0Kstar0mumu[41];
		else if(!strcmp(obs->type,"P8primeCP")) return obsB0Kstar0mumu[42];
		else if(!strcmp(obs->type,"ATIm")) return -2.*obsB0Kstar0mumu[15];
		else if(!strcmp(obs->type,"ATImCP")) return -2.*obsB0Kstar0mumu[38];
		else if(!strcmp(obs->type,"S6c")) return obsB0Kstar0mumu[43];
		else if(!strcmp(obs->type,"AFBfull")) return 3./8.*obsB0Kstar0mumu[43]+obsB0Kstar0mumu[1];
		else if(!strcmp(obs->type,"A6s")) return obsB0Kstar0mumu[44];
		else if(!strcmp(obs->type,"ACP")) return obsB0Kstar0mumu[45];
		else if(!strcmp(obs->type,"S1c")) return obsB0Kstar0mumu[46];
		else if(!strcmp(obs->type,"S2c")) return obsB0Kstar0mumu[47];
		else if(!strcmp(obs->type,"S2s")) return obsB0Kstar0mumu[48];
		else if(!strcmp(obs->type,"AFL")) return obsB0Kstar0mumu[49];
		else if(!strcmp(obs->type,"A1c")) return obsB0Kstar0mumu[50];
		else if(!strcmp(obs->type,"A2s")) return obsB0Kstar0mumu[51];
		else if(!strcmp(obs->type,"AFBCP")) return obsB0Kstar0mumu[52];
		else if(!strcmp(obs->type,"A6c")) return obsB0Kstar0mumu[53];
		else if(!strcmp(obs->type,"AFBCPfull")) return  3./8.*obsB0Kstar0mumu[53]+obsB0Kstar0mumu[52];
		/* Form Factors */
		else if(!strcmp(obs->type,"V"))   return obsB0Kstar0mumu[Nobs_BKsll_phys+1];
		else if(!strcmp(obs->type,"A0"))  return obsB0Kstar0mumu[Nobs_BKsll_phys+2];
		else if(!strcmp(obs->type,"A1"))  return obsB0Kstar0mumu[Nobs_BKsll_phys+3];
		else if(!strcmp(obs->type,"A12")) return obsB0Kstar0mumu[Nobs_BKsll_phys+4];
		else if(!strcmp(obs->type,"T1"))  return obsB0Kstar0mumu[Nobs_BKsll_phys+5];
		else if(!strcmp(obs->type,"T2"))  return obsB0Kstar0mumu[Nobs_BKsll_phys+6];
		else if(!strcmp(obs->type,"T23")) return obsB0Kstar0mumu[Nobs_BKsll_phys+7];
		else if(!strcmp(obs->type,"A2"))  return obsB0Kstar0mumu[Nobs_BKsll_phys+8];
		/* Results from Hadronic Fit */
		else if(!strcmp(obs->type,"hplusRe"))   return obsB0Kstar0mumu[Nobs_BKsll_phys+9];
		else if(!strcmp(obs->type,"hplusIm"))  return obsB0Kstar0mumu[Nobs_BKsll_phys+10];
		else if(!strcmp(obs->type,"hminusRe"))   return obsB0Kstar0mumu[Nobs_BKsll_phys+11];
		else if(!strcmp(obs->type,"hminusIm"))  return obsB0Kstar0mumu[Nobs_BKsll_phys+12];
		else if(!strcmp(obs->type,"hzeroRe"))   return obsB0Kstar0mumu[Nobs_BKsll_phys+13];
		else if(!strcmp(obs->type,"hzeroIm"))  return obsB0Kstar0mumu[Nobs_BKsll_phys+14];
		/*************************************************************/

	}

	else if(!strcmp(obs->decay,"BKstarmumu"))
	{
		if(reset||(fabs(binlowBKstarmumu-obs->low)>1.e-5)||(fabs(binhighBKstarmumu-obs->high)>1.e-5))
		{
			obsBKstarmumu[0]=BRBKstarll(2,1,obs->low,obs->high,obsBKstarmumu,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b)/(obs->high-obs->low);
		
			binlowBKstarmumu=obs->low;
			binhighBKstarmumu=obs->high;
		}
		
		if(!strcmp(obs->type,"dGamma/dq2")) return obsBKstarmumu[0];
		else if(!strcmp(obs->type,"BR")) return obsBKstarmumu[0]*(obs->high-obs->low);
		else if(!strcmp(obs->type,"AFB")) return obsBKstarmumu[1];
		else if(!strcmp(obs->type,"FL")) return obsBKstarmumu[2];
		else if(!strcmp(obs->type,"FT")) return obsBKstarmumu[3];
		else if(!strcmp(obs->type,"AT1")) return obsBKstarmumu[4];
		else if(!strcmp(obs->type,"AT2")||!strcmp(obs->type,"P1")) return obsBKstarmumu[5];
		else if(!strcmp(obs->type,"AT3")) return obsBKstarmumu[6];
		else if(!strcmp(obs->type,"AT4")) return obsBKstarmumu[7];
		else if(!strcmp(obs->type,"AT5")) return obsBKstarmumu[8];
		else if(!strcmp(obs->type,"HT1")) return obsBKstarmumu[9];
		else if(!strcmp(obs->type,"HT2")) return obsBKstarmumu[10];
		else if(!strcmp(obs->type,"HT3")) return obsBKstarmumu[11];
		else if(!strcmp(obs->type,"alpha")) return obsBKstarmumu[12];
		else if(!strcmp(obs->type,"AIm")) return obsBKstarmumu[13];
		else if(!strcmp(obs->type,"P2")) return obsBKstarmumu[14];
		else if(!strcmp(obs->type,"P3")) return obsBKstarmumu[15];
		else if(!strcmp(obs->type,"P6")) return obsBKstarmumu[16];
		else if(!strcmp(obs->type,"P4prime")) return obsBKstarmumu[17];
		else if(!strcmp(obs->type,"P5prime")) return obsBKstarmumu[18];
		else if(!strcmp(obs->type,"P6prime")) return obsBKstarmumu[19];
		else if(!strcmp(obs->type,"P8")) return obsBKstarmumu[20];
		else if(!strcmp(obs->type,"P8prime")) return obsBKstarmumu[21];
		else if(!strcmp(obs->type,"A7")) return obsBKstarmumu[22];
		else if(!strcmp(obs->type,"A8")) return obsBKstarmumu[23];
		else if(!strcmp(obs->type,"A9")) return obsBKstarmumu[24];
		else if(!strcmp(obs->type,"S3")) return obsBKstarmumu[25];
		else if(!strcmp(obs->type,"S4")) return obsBKstarmumu[26];
		else if(!strcmp(obs->type,"S5")) return obsBKstarmumu[27];
		else if(!strcmp(obs->type,"S7")) return obsBKstarmumu[28];
		else if(!strcmp(obs->type,"S8")) return obsBKstarmumu[29];
		else if(!strcmp(obs->type,"S9")) return obsBKstarmumu[30];
		else if(!strcmp(obs->type,"A3")) return obsBKstarmumu[31];
		else if(!strcmp(obs->type,"A4")) return obsBKstarmumu[32];
		else if(!strcmp(obs->type,"A5")) return obsBKstarmumu[33];
		else if(!strcmp(obs->type,"ATRe")) return obsBKstarmumu[34];
		else if(!strcmp(obs->type,"ATReCP")) return obsBKstarmumu[35];
		else if(!strcmp(obs->type,"AT2CP")||!strcmp(obs->type,"P1CP")) return obsBKstarmumu[36];
		else if(!strcmp(obs->type,"P2CP")) return obsBKstarmumu[37];
		else if(!strcmp(obs->type,"P3CP")) return obsBKstarmumu[38];
		else if(!strcmp(obs->type,"P4primeCP")) return obsBKstarmumu[39];
		else if(!strcmp(obs->type,"P5primeCP")) return obsBKstarmumu[40];
		else if(!strcmp(obs->type,"P6primeCP")) return obsBKstarmumu[41];
		else if(!strcmp(obs->type,"P8primeCP")) return obsBKstarmumu[42];
		else if(!strcmp(obs->type,"ATIm")) return -2.*obsBKstarmumu[15];
		else if(!strcmp(obs->type,"ATImCP")) return -2.*obsBKstarmumu[38];
		else if(!strcmp(obs->type,"S6c")) return obsBKstarmumu[43];
		else if(!strcmp(obs->type,"AFBfull")) return 3./8.*obsBKstarmumu[43]+obsBKstarmumu[1];
		else if(!strcmp(obs->type,"A6s")) return obsBKstarmumu[44];
		else if(!strcmp(obs->type,"ACP")) return obsBKstarmumu[45];
		else if(!strcmp(obs->type,"S1c")) return obsBKstarmumu[46];
		else if(!strcmp(obs->type,"S2c")) return obsBKstarmumu[47];
		else if(!strcmp(obs->type,"S2s")) return obsBKstarmumu[48];
		else if(!strcmp(obs->type,"AFL")) return obsBKstarmumu[49];
		else if(!strcmp(obs->type,"A1c")) return obsBKstarmumu[50];
		else if(!strcmp(obs->type,"A2s")) return obsBKstarmumu[51];
		else if(!strcmp(obs->type,"AFBCP")) return obsBKstarmumu[52];
		else if(!strcmp(obs->type,"A6c")) return obsBKstarmumu[53];
		else if(!strcmp(obs->type,"AFBCPfull")) return  3./8.*obsBKstarmumu[53]+obsBKstarmumu[52];
		else if(!strcmp(obs->type,"V"))   return obsBKstarmumu[Nobs_BKsll_phys+1];
		else if(!strcmp(obs->type,"A0"))  return obsBKstarmumu[Nobs_BKsll_phys+2];
		else if(!strcmp(obs->type,"A1"))  return obsBKstarmumu[Nobs_BKsll_phys+3];
		else if(!strcmp(obs->type,"A12")) return obsBKstarmumu[Nobs_BKsll_phys+4];
		else if(!strcmp(obs->type,"T1"))  return obsBKstarmumu[Nobs_BKsll_phys+5];
		else if(!strcmp(obs->type,"T2"))  return obsBKstarmumu[Nobs_BKsll_phys+6];
		else if(!strcmp(obs->type,"T23")) return obsBKstarmumu[Nobs_BKsll_phys+7];
		else if(!strcmp(obs->type,"A2"))  return obsBKstarmumu[Nobs_BKsll_phys+8];
		/*************************************************************/
	}
	
	else if(!strcmp(obs->decay,"B0Kstar0ee"))
	{
		if(reset||(fabs(binlowB0Kstar0ee-obs->low)>1.e-5)||(fabs(binhighB0Kstar0ee-obs->high)>1.e-5))
		{
			obsB0Kstar0ee[0]=BRBKstarll(1,0,obs->low,obs->high,obsB0Kstar0ee,C0eb,C1eb,C2eb,CQ0eb,CQ1eb,Cpeb,CQpeb,param,mu_b)/(obs->high-obs->low);
		
			binlowB0Kstar0ee=obs->low;
			binhighB0Kstar0ee=obs->high;
		}
		
		if(!strcmp(obs->type,"dGamma/dq2")) return obsB0Kstar0ee[0];
		else if(!strcmp(obs->type,"BR")) return obsB0Kstar0ee[0]*(obs->high-obs->low);
		else if(!strcmp(obs->type,"AFB")) return obsB0Kstar0ee[1];
		else if(!strcmp(obs->type,"FL")) return obsB0Kstar0ee[2];
		else if(!strcmp(obs->type,"FT")) return obsB0Kstar0ee[3];
		else if(!strcmp(obs->type,"AT1")) return obsB0Kstar0ee[4];
		else if(!strcmp(obs->type,"AT2")||!strcmp(obs->type,"P1")) return obsB0Kstar0ee[5];
		else if(!strcmp(obs->type,"AT3")) return obsB0Kstar0ee[6];
		else if(!strcmp(obs->type,"AT4")) return obsB0Kstar0ee[7];
		else if(!strcmp(obs->type,"AT5")) return obsB0Kstar0ee[8];
		else if(!strcmp(obs->type,"HT1")) return obsB0Kstar0ee[9];
		else if(!strcmp(obs->type,"HT2")) return obsB0Kstar0ee[10];
		else if(!strcmp(obs->type,"HT3")) return obsB0Kstar0ee[11];
		else if(!strcmp(obs->type,"alpha")) return obsB0Kstar0ee[12];
		else if(!strcmp(obs->type,"AIm")) return obsB0Kstar0ee[13];
		else if(!strcmp(obs->type,"P2")) return obsB0Kstar0ee[14];
		else if(!strcmp(obs->type,"P3")) return obsB0Kstar0ee[15];
		else if(!strcmp(obs->type,"P6")) return obsB0Kstar0ee[16];
		else if(!strcmp(obs->type,"P4prime")) return obsB0Kstar0ee[17];
		else if(!strcmp(obs->type,"P5prime")) return obsB0Kstar0ee[18];
		else if(!strcmp(obs->type,"P6prime")) return obsB0Kstar0ee[19];
		else if(!strcmp(obs->type,"P8")) return obsB0Kstar0ee[20];
		else if(!strcmp(obs->type,"P8prime")) return obsB0Kstar0ee[21];
		else if(!strcmp(obs->type,"A7")) return obsB0Kstar0ee[22];
		else if(!strcmp(obs->type,"A8")) return obsB0Kstar0ee[23];
		else if(!strcmp(obs->type,"A9")) return obsB0Kstar0ee[24];
		else if(!strcmp(obs->type,"S3")) return obsB0Kstar0ee[25];
		else if(!strcmp(obs->type,"S4")) return obsB0Kstar0ee[26];
		else if(!strcmp(obs->type,"S5")) return obsB0Kstar0ee[27];
		else if(!strcmp(obs->type,"S7")) return obsB0Kstar0ee[28];
		else if(!strcmp(obs->type,"S8")) return obsB0Kstar0ee[29];
		else if(!strcmp(obs->type,"S9")) return obsB0Kstar0ee[30];
		else if(!strcmp(obs->type,"A3")) return obsB0Kstar0ee[31];
		else if(!strcmp(obs->type,"A4")) return obsB0Kstar0ee[32];
		else if(!strcmp(obs->type,"A5")) return obsB0Kstar0ee[33];
		else if(!strcmp(obs->type,"ATRe")) return obsB0Kstar0ee[34];
		else if(!strcmp(obs->type,"ATReCP")) return obsB0Kstar0ee[35];
		else if(!strcmp(obs->type,"AT2CP")||!strcmp(obs->type,"P1CP")) return obsB0Kstar0ee[36];
		else if(!strcmp(obs->type,"P2CP")) return obsB0Kstar0ee[37];
		else if(!strcmp(obs->type,"P3CP")) return obsB0Kstar0ee[38];
		else if(!strcmp(obs->type,"P4primeCP")) return obsB0Kstar0ee[39];
		else if(!strcmp(obs->type,"P5primeCP")) return obsB0Kstar0ee[40];
		else if(!strcmp(obs->type,"P6primeCP")) return obsB0Kstar0ee[41];
		else if(!strcmp(obs->type,"P8primeCP")) return obsB0Kstar0ee[42];
		else if(!strcmp(obs->type,"ATIm")) return -2.*obsB0Kstar0ee[15];
		else if(!strcmp(obs->type,"ATImCP")) return -2.*obsB0Kstar0ee[38];
		else if(!strcmp(obs->type,"S6c")) return obsB0Kstar0ee[43];
		else if(!strcmp(obs->type,"AFBfull")) return 3./8.*obsB0Kstar0ee[43]+obsB0Kstar0ee[1];
		else if(!strcmp(obs->type,"A6s")) return obsB0Kstar0ee[44];
		else if(!strcmp(obs->type,"ACP")) return obsB0Kstar0ee[45];
		else if(!strcmp(obs->type,"S1c")) return obsB0Kstar0ee[46];
		else if(!strcmp(obs->type,"S2c")) return obsB0Kstar0ee[47];
		else if(!strcmp(obs->type,"S2s")) return obsB0Kstar0ee[48];
		else if(!strcmp(obs->type,"AFL")) return obsB0Kstar0ee[49];
		else if(!strcmp(obs->type,"A1c")) return obsB0Kstar0ee[50];
		else if(!strcmp(obs->type,"A2s")) return obsB0Kstar0ee[51];
		else if(!strcmp(obs->type,"AFBCP")) return obsB0Kstar0ee[52];
		else if(!strcmp(obs->type,"A6c")) return obsB0Kstar0ee[53];
		else if(!strcmp(obs->type,"AFBCPfull")) return  3./8.*obsB0Kstar0ee[53]+obsB0Kstar0ee[52];
		else if(!strcmp(obs->type,"V"))   return obsB0Kstar0ee[Nobs_BKsll_phys+1];
		else if(!strcmp(obs->type,"A0"))  return obsB0Kstar0ee[Nobs_BKsll_phys+2];
		else if(!strcmp(obs->type,"A1"))  return obsB0Kstar0ee[Nobs_BKsll_phys+3];
		else if(!strcmp(obs->type,"A12")) return obsB0Kstar0ee[Nobs_BKsll_phys+4];
		else if(!strcmp(obs->type,"T1"))  return obsB0Kstar0ee[Nobs_BKsll_phys+5];
		else if(!strcmp(obs->type,"T2"))  return obsB0Kstar0ee[Nobs_BKsll_phys+6];
		else if(!strcmp(obs->type,"T23")) return obsB0Kstar0ee[Nobs_BKsll_phys+7];
		else if(!strcmp(obs->type,"A2"))  return obsB0Kstar0ee[Nobs_BKsll_phys+8];
	}

	else if(!strcmp(obs->decay,"BKstaree"))
	{
		if(reset||(fabs(binlowBKstaree-obs->low)>1.e-5)||(fabs(binhighBKstaree-obs->high)>1.e-5))
		{
			obsBKstaree[0]=BRBKstarll(1,1,obs->low,obs->high,obsBKstaree,C0eb,C1eb,C2eb,CQ0eb,CQ1eb,Cpeb,CQpeb,param,mu_b)/(obs->high-obs->low);
		
			binlowBKstaree=obs->low;
			binhighBKstaree=obs->high;
		}
		
		if(!strcmp(obs->type,"dGamma/dq2")) return obsBKstaree[0];
		else if(!strcmp(obs->type,"BR")) return obsBKstaree[0]*(obs->high-obs->low);
		else if(!strcmp(obs->type,"AFB")) return obsBKstaree[1];
		else if(!strcmp(obs->type,"FL")) return obsBKstaree[2];
		else if(!strcmp(obs->type,"FT")) return obsBKstaree[3];
		else if(!strcmp(obs->type,"AT1")) return obsBKstaree[4];
		else if(!strcmp(obs->type,"AT2")||!strcmp(obs->type,"P1")) return obsBKstaree[5];
		else if(!strcmp(obs->type,"AT3")) return obsBKstaree[6];
		else if(!strcmp(obs->type,"AT4")) return obsBKstaree[7];
		else if(!strcmp(obs->type,"AT5")) return obsBKstaree[8];
		else if(!strcmp(obs->type,"HT1")) return obsBKstaree[9];
		else if(!strcmp(obs->type,"HT2")) return obsBKstaree[10];
		else if(!strcmp(obs->type,"HT3")) return obsBKstaree[11];
		else if(!strcmp(obs->type,"alpha")) return obsBKstaree[12];
		else if(!strcmp(obs->type,"AIm")) return obsBKstaree[13];
		else if(!strcmp(obs->type,"P2")) return obsBKstaree[14];
		else if(!strcmp(obs->type,"P3")) return obsBKstaree[15];
		else if(!strcmp(obs->type,"P6")) return obsBKstaree[16];
		else if(!strcmp(obs->type,"P4prime")) return obsBKstaree[17];
		else if(!strcmp(obs->type,"P5prime")) return obsBKstaree[18];
		else if(!strcmp(obs->type,"P6prime")) return obsBKstaree[19];
		else if(!strcmp(obs->type,"P8")) return obsBKstaree[20];
		else if(!strcmp(obs->type,"P8prime")) return obsBKstaree[21];
		else if(!strcmp(obs->type,"A7")) return obsBKstaree[22];
		else if(!strcmp(obs->type,"A8")) return obsBKstaree[23];
		else if(!strcmp(obs->type,"A9")) return obsBKstaree[24];
		else if(!strcmp(obs->type,"S3")) return obsBKstaree[25];
		else if(!strcmp(obs->type,"S4")) return obsBKstaree[26];
		else if(!strcmp(obs->type,"S5")) return obsBKstaree[27];
		else if(!strcmp(obs->type,"S7")) return obsBKstaree[28];
		else if(!strcmp(obs->type,"S8")) return obsBKstaree[29];
		else if(!strcmp(obs->type,"S9")) return obsBKstaree[30];
		else if(!strcmp(obs->type,"A3")) return obsBKstaree[31];
		else if(!strcmp(obs->type,"A4")) return obsBKstaree[32];
		else if(!strcmp(obs->type,"A5")) return obsBKstaree[33];
		else if(!strcmp(obs->type,"ATRe")) return obsBKstaree[34];
		else if(!strcmp(obs->type,"ATReCP")) return obsBKstaree[35];
		else if(!strcmp(obs->type,"AT2CP")||!strcmp(obs->type,"P1CP")) return obsBKstaree[36];
		else if(!strcmp(obs->type,"P2CP")) return obsBKstaree[37];
		else if(!strcmp(obs->type,"P3CP")) return obsBKstaree[38];
		else if(!strcmp(obs->type,"P4primeCP")) return obsBKstaree[39];
		else if(!strcmp(obs->type,"P5primeCP")) return obsBKstaree[40];
		else if(!strcmp(obs->type,"P6primeCP")) return obsBKstaree[41];
		else if(!strcmp(obs->type,"P8primeCP")) return obsBKstaree[42];
		else if(!strcmp(obs->type,"ATIm")) return -2.*obsBKstaree[15];
		else if(!strcmp(obs->type,"ATImCP")) return -2.*obsBKstaree[38];
		else if(!strcmp(obs->type,"S6c")) return obsBKstaree[43];
		else if(!strcmp(obs->type,"AFBfull")) return 3./8.*obsBKstaree[43]+obsBKstaree[1];
		else if(!strcmp(obs->type,"A6s")) return obsBKstaree[44];
		else if(!strcmp(obs->type,"ACP")) return obsBKstaree[45];
		else if(!strcmp(obs->type,"S1c")) return obsBKstaree[46];
		else if(!strcmp(obs->type,"S2c")) return obsBKstaree[47];
		else if(!strcmp(obs->type,"S2s")) return obsBKstaree[48];
		else if(!strcmp(obs->type,"AFL")) return obsBKstaree[49];
		else if(!strcmp(obs->type,"A1c")) return obsBKstaree[50];
		else if(!strcmp(obs->type,"A2s")) return obsBKstaree[51];
		else if(!strcmp(obs->type,"AFBCP")) return obsBKstaree[52];
		else if(!strcmp(obs->type,"A6c")) return obsBKstaree[53];
		else if(!strcmp(obs->type,"AFBCPfull")) return  3./8.*obsBKstaree[53]+obsBKstaree[52];
		else if(!strcmp(obs->type,"V"))   return obsBKstaree[Nobs_BKsll_phys+1];
		else if(!strcmp(obs->type,"A0"))  return obsBKstaree[Nobs_BKsll_phys+2];
		else if(!strcmp(obs->type,"A1"))  return obsBKstaree[Nobs_BKsll_phys+3];
		else if(!strcmp(obs->type,"A12")) return obsBKstaree[Nobs_BKsll_phys+4];
		else if(!strcmp(obs->type,"T1"))  return obsBKstaree[Nobs_BKsll_phys+5];
		else if(!strcmp(obs->type,"T2"))  return obsBKstaree[Nobs_BKsll_phys+6];
		else if(!strcmp(obs->type,"T23")) return obsBKstaree[Nobs_BKsll_phys+7];
		else if(!strcmp(obs->type,"A2"))  return obsBKstaree[Nobs_BKsll_phys+8];
	}

	else if(!strcmp(obs->decay,"B0Kstar0tautau"))
	{
		if(reset||(fabs(binlowB0Kstar0tautau-obs->low)>1.e-5)||(fabs(binhighB0Kstar0tautau-obs->high)>1.e-5))
		{
			obsB0Kstar0tautau[0]=BRBKstarll(3,0,obs->low,obs->high,obsB0Kstar0tautau,C0taub,C1taub,C2taub,CQ0taub,CQ1taub,Cptaub,CQptaub,param,mu_b)/(obs->high-obs->low);
		
			binlowB0Kstar0tautau=obs->low;
			binhighB0Kstar0tautau=obs->high;
		}
		
		if(!strcmp(obs->type,"dGamma/dq2")) return obsB0Kstar0tautau[0];
		else if(!strcmp(obs->type,"BR")) return obsB0Kstar0tautau[0]*(obs->high-obs->low);
		else if(!strcmp(obs->type,"AFB")) return obsB0Kstar0tautau[1];
		else if(!strcmp(obs->type,"FL")) return obsB0Kstar0tautau[2];
		else if(!strcmp(obs->type,"FT")) return obsB0Kstar0tautau[3];
		else if(!strcmp(obs->type,"AT1")) return obsB0Kstar0tautau[4];
		else if(!strcmp(obs->type,"AT2")||!strcmp(obs->type,"P1")) return obsB0Kstar0tautau[5];
		else if(!strcmp(obs->type,"AT3")) return obsB0Kstar0tautau[6];
		else if(!strcmp(obs->type,"AT4")) return obsB0Kstar0tautau[7];
		else if(!strcmp(obs->type,"AT5")) return obsB0Kstar0tautau[8];
		else if(!strcmp(obs->type,"HT1")) return obsB0Kstar0tautau[9];
		else if(!strcmp(obs->type,"HT2")) return obsB0Kstar0tautau[10];
		else if(!strcmp(obs->type,"HT3")) return obsB0Kstar0tautau[11];
		else if(!strcmp(obs->type,"alpha")) return obsB0Kstar0tautau[12];
		else if(!strcmp(obs->type,"AIm")) return obsB0Kstar0tautau[13];
		else if(!strcmp(obs->type,"P2")) return obsB0Kstar0tautau[14];
		else if(!strcmp(obs->type,"P3")) return obsB0Kstar0tautau[15];
		else if(!strcmp(obs->type,"P6")) return obsB0Kstar0tautau[16];
		else if(!strcmp(obs->type,"P4prime")) return obsB0Kstar0tautau[17];
		else if(!strcmp(obs->type,"P5prime")) return obsB0Kstar0tautau[18];
		else if(!strcmp(obs->type,"P6prime")) return obsB0Kstar0tautau[19];
		else if(!strcmp(obs->type,"P8")) return obsB0Kstar0tautau[20];
		else if(!strcmp(obs->type,"P8prime")) return obsB0Kstar0tautau[21];
		else if(!strcmp(obs->type,"A7")) return obsB0Kstar0tautau[22];
		else if(!strcmp(obs->type,"A8")) return obsB0Kstar0tautau[23];
		else if(!strcmp(obs->type,"A9")) return obsB0Kstar0tautau[24];
		else if(!strcmp(obs->type,"S3")) return obsB0Kstar0tautau[25];
		else if(!strcmp(obs->type,"S4")) return obsB0Kstar0tautau[26];
		else if(!strcmp(obs->type,"S5")) return obsB0Kstar0tautau[27];
		else if(!strcmp(obs->type,"S7")) return obsB0Kstar0tautau[28];
		else if(!strcmp(obs->type,"S8")) return obsB0Kstar0tautau[29];
		else if(!strcmp(obs->type,"S9")) return obsB0Kstar0tautau[30];
		else if(!strcmp(obs->type,"A3")) return obsB0Kstar0tautau[31];
		else if(!strcmp(obs->type,"A4")) return obsB0Kstar0tautau[32];
		else if(!strcmp(obs->type,"A5")) return obsB0Kstar0tautau[33];
		else if(!strcmp(obs->type,"ATRe")) return obsB0Kstar0tautau[34];
		else if(!strcmp(obs->type,"ATReCP")) return obsB0Kstar0tautau[35];
		else if(!strcmp(obs->type,"AT2CP")||!strcmp(obs->type,"P1CP")) return obsB0Kstar0tautau[36];
		else if(!strcmp(obs->type,"P2CP")) return obsB0Kstar0tautau[37];
		else if(!strcmp(obs->type,"P3CP")) return obsB0Kstar0tautau[38];
		else if(!strcmp(obs->type,"P4primeCP")) return obsB0Kstar0tautau[39];
		else if(!strcmp(obs->type,"P5primeCP")) return obsB0Kstar0tautau[40];
		else if(!strcmp(obs->type,"P6primeCP")) return obsB0Kstar0tautau[41];
		else if(!strcmp(obs->type,"P8primeCP")) return obsB0Kstar0tautau[42];
		else if(!strcmp(obs->type,"ATIm")) return -2.*obsB0Kstar0tautau[15];
		else if(!strcmp(obs->type,"ATImCP")) return -2.*obsB0Kstar0tautau[38];
		else if(!strcmp(obs->type,"S6c")) return obsB0Kstar0tautau[43];
		else if(!strcmp(obs->type,"AFBfull")) return 3./8.*obsB0Kstar0tautau[43]+obsB0Kstar0tautau[1];
		else if(!strcmp(obs->type,"A6s")) return obsB0Kstar0tautau[44];
		else if(!strcmp(obs->type,"ACP")) return obsB0Kstar0tautau[45];
		else if(!strcmp(obs->type,"S1c")) return obsB0Kstar0tautau[46];
		else if(!strcmp(obs->type,"S2c")) return obsB0Kstar0tautau[47];
		else if(!strcmp(obs->type,"S2s")) return obsB0Kstar0tautau[48];
		else if(!strcmp(obs->type,"AFL")) return obsB0Kstar0tautau[49];
		else if(!strcmp(obs->type,"A1c")) return obsB0Kstar0tautau[50];
		else if(!strcmp(obs->type,"A2s")) return obsB0Kstar0tautau[51];
		else if(!strcmp(obs->type,"AFBCP")) return obsB0Kstar0tautau[52];
		else if(!strcmp(obs->type,"A6c")) return obsB0Kstar0tautau[53];
		else if(!strcmp(obs->type,"AFBCPfull")) return  3./8.*obsB0Kstar0tautau[53]+obsB0Kstar0tautau[52];
		else if(!strcmp(obs->type,"V"))   return obsB0Kstar0tautau[Nobs_BKsll_phys+1];
		else if(!strcmp(obs->type,"A0"))  return obsB0Kstar0tautau[Nobs_BKsll_phys+2];
		else if(!strcmp(obs->type,"A1"))  return obsB0Kstar0tautau[Nobs_BKsll_phys+3];
		else if(!strcmp(obs->type,"A12")) return obsB0Kstar0tautau[Nobs_BKsll_phys+4];
		else if(!strcmp(obs->type,"T1"))  return obsB0Kstar0tautau[Nobs_BKsll_phys+5];
		else if(!strcmp(obs->type,"T2"))  return obsB0Kstar0tautau[Nobs_BKsll_phys+6];
		else if(!strcmp(obs->type,"T23")) return obsB0Kstar0tautau[Nobs_BKsll_phys+7];
		else if(!strcmp(obs->type,"A2"))  return obsB0Kstar0tautau[Nobs_BKsll_phys+8];
	}

	else if(!strcmp(obs->decay,"BKstartautau"))
	{
		if(reset||(fabs(binlowBKstartautau-obs->low)>1.e-5)||(fabs(binhighBKstartautau-obs->high)>1.e-5))
		{
			obsBKstartautau[0]=BRBKstarll(3,1,obs->low,obs->high,obsBKstartautau,C0taub,C1taub,C2taub,CQ0taub,CQ1taub,Cptaub,CQptaub,param,mu_b)/(obs->high-obs->low);
		
			binlowBKstartautau=obs->low;
			binhighBKstartautau=obs->high;
		}
		
		if(!strcmp(obs->type,"dGamma/dq2")) return obsBKstartautau[0];
		else if(!strcmp(obs->type,"BR")) return obsBKstartautau[0]*(obs->high-obs->low);
		else if(!strcmp(obs->type,"AFB")) return obsBKstartautau[1];
		else if(!strcmp(obs->type,"FL")) return obsBKstartautau[2];
		else if(!strcmp(obs->type,"FT")) return obsBKstartautau[3];
		else if(!strcmp(obs->type,"AT1")) return obsBKstartautau[4];
		else if(!strcmp(obs->type,"AT2")||!strcmp(obs->type,"P1")) return obsBKstartautau[5];
		else if(!strcmp(obs->type,"AT3")) return obsBKstartautau[6];
		else if(!strcmp(obs->type,"AT4")) return obsBKstartautau[7];
		else if(!strcmp(obs->type,"AT5")) return obsBKstartautau[8];
		else if(!strcmp(obs->type,"HT1")) return obsBKstartautau[9];
		else if(!strcmp(obs->type,"HT2")) return obsBKstartautau[10];
		else if(!strcmp(obs->type,"HT3")) return obsBKstartautau[11];
		else if(!strcmp(obs->type,"alpha")) return obsBKstartautau[12];
		else if(!strcmp(obs->type,"AIm")) return obsBKstartautau[13];
		else if(!strcmp(obs->type,"P2")) return obsBKstartautau[14];
		else if(!strcmp(obs->type,"P3")) return obsBKstartautau[15];
		else if(!strcmp(obs->type,"P6")) return obsBKstartautau[16];
		else if(!strcmp(obs->type,"P4prime")) return obsBKstartautau[17];
		else if(!strcmp(obs->type,"P5prime")) return obsBKstartautau[18];
		else if(!strcmp(obs->type,"P6prime")) return obsBKstartautau[19];
		else if(!strcmp(obs->type,"P8")) return obsBKstartautau[20];
		else if(!strcmp(obs->type,"P8prime")) return obsBKstartautau[21];
		else if(!strcmp(obs->type,"A7")) return obsBKstartautau[22];
		else if(!strcmp(obs->type,"A8")) return obsBKstartautau[23];
		else if(!strcmp(obs->type,"A9")) return obsBKstartautau[24];
		else if(!strcmp(obs->type,"S3")) return obsBKstartautau[25];
		else if(!strcmp(obs->type,"S4")) return obsBKstartautau[26];
		else if(!strcmp(obs->type,"S5")) return obsBKstartautau[27];
		else if(!strcmp(obs->type,"S7")) return obsBKstartautau[28];
		else if(!strcmp(obs->type,"S8")) return obsBKstartautau[29];
		else if(!strcmp(obs->type,"S9")) return obsBKstartautau[30];
		else if(!strcmp(obs->type,"A3")) return obsBKstartautau[31];
		else if(!strcmp(obs->type,"A4")) return obsBKstartautau[32];
		else if(!strcmp(obs->type,"A5")) return obsBKstartautau[33];
		else if(!strcmp(obs->type,"ATRe")) return obsBKstartautau[34];
		else if(!strcmp(obs->type,"ATReCP")) return obsBKstartautau[35];
		else if(!strcmp(obs->type,"AT2CP")||!strcmp(obs->type,"P1CP")) return obsBKstartautau[36];
		else if(!strcmp(obs->type,"P2CP")) return obsBKstartautau[37];
		else if(!strcmp(obs->type,"P3CP")) return obsBKstartautau[38];
		else if(!strcmp(obs->type,"P4primeCP")) return obsBKstartautau[39];
		else if(!strcmp(obs->type,"P5primeCP")) return obsBKstartautau[40];
		else if(!strcmp(obs->type,"P6primeCP")) return obsBKstartautau[41];
		else if(!strcmp(obs->type,"P8primeCP")) return obsBKstartautau[42];
		else if(!strcmp(obs->type,"ATIm")) return -2.*obsBKstartautau[15];
		else if(!strcmp(obs->type,"ATImCP")) return -2.*obsBKstartautau[38];
		else if(!strcmp(obs->type,"S6c")) return obsBKstartautau[43];
		else if(!strcmp(obs->type,"AFBfull")) return 3./8.*obsBKstartautau[43]+obsBKstartautau[1];
		else if(!strcmp(obs->type,"A6s")) return obsBKstartautau[44];
		else if(!strcmp(obs->type,"ACP")) return obsBKstartautau[45];
		else if(!strcmp(obs->type,"S1c")) return obsBKstartautau[46];
		else if(!strcmp(obs->type,"S2c")) return obsBKstartautau[47];
		else if(!strcmp(obs->type,"S2s")) return obsBKstartautau[48];
		else if(!strcmp(obs->type,"AFL")) return obsBKstartautau[49];
		else if(!strcmp(obs->type,"A1c")) return obsBKstartautau[50];
		else if(!strcmp(obs->type,"A2s")) return obsBKstartautau[51];
		else if(!strcmp(obs->type,"AFBCP")) return obsBKstartautau[52];
		else if(!strcmp(obs->type,"A6c")) return obsBKstartautau[53];
		else if(!strcmp(obs->type,"AFBCPfull")) return  3./8.*obsBKstartautau[53]+obsBKstartautau[52];
		else if(!strcmp(obs->type,"V"))   return obsBKstartautau[Nobs_BKsll_phys+1];
		else if(!strcmp(obs->type,"A0"))  return obsBKstartautau[Nobs_BKsll_phys+2];
		else if(!strcmp(obs->type,"A1"))  return obsBKstartautau[Nobs_BKsll_phys+3];
		else if(!strcmp(obs->type,"A12")) return obsBKstartautau[Nobs_BKsll_phys+4];
		else if(!strcmp(obs->type,"T1"))  return obsBKstartautau[Nobs_BKsll_phys+5];
		else if(!strcmp(obs->type,"T2"))  return obsBKstartautau[Nobs_BKsll_phys+6];
		else if(!strcmp(obs->type,"T23")) return obsBKstartautau[Nobs_BKsll_phys+7];
		else if(!strcmp(obs->type,"A2"))  return obsBKstartautau[Nobs_BKsll_phys+8];
	}

	else if(!strcmp(obs->type,"R-1")&&!strcmp(obs->decay,"B0Kstar0ll"))
	{
		binlowB0Kstar0mumu=binlowB0Kstar0ee=obs->low;
		binhighB0Kstar0mumu=binhighB0Kstar0ee=obs->high;
		
		obsB0Kstar0mumu[0]=BRBKstarll(2,0,obs->low,obs->high,obsB0Kstar0mumu,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b)/(obs->high-obs->low);
		obsB0Kstar0ee[0]=BRBKstarll(1,0,obs->low,obs->high,obsB0Kstar0ee,C0eb,C1eb,C2eb,CQ0eb,CQ1eb,Cpeb,CQpeb,param,mu_b)/(obs->high-obs->low);
		
		return obsB0Kstar0mumu[0]/obsB0Kstar0ee[0]-1.;
	}

	else if(!strcmp(obs->type,"R-1")&&!strcmp(obs->decay,"BKstarll"))
	{
		binlowBKstarmumu=binlowBKstaree=obs->low;
		binhighBKstarmumu=binhighBKstaree=obs->high;
		
		obsBKstarmumu[0]=BRBKstarll(2,1,obs->low,obs->high,obsBKstarmumu,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b)/(obs->high-obs->low);
		obsBKstaree[0]=BRBKstarll(1,1,obs->low,obs->high,obsBKstaree,C0eb,C1eb,C2eb,CQ0eb,CQ1eb,Cpeb,CQpeb,param,mu_b)/(obs->high-obs->low);
	
		return obsBKstarmumu[0]/obsBKstaree[0]-1.;
	}
	
/* B -> K l+ l- */

	else if(!strcmp(obs->decay,"B0K0mumu"))
	{
		if(reset||(fabs(binlowB0K0mumu-obs->low)>1.e-5)||(fabs(binhighB0K0mumu-obs->high)>1.e-5))
		{
			obsB0K0mumu[0]=BRBKll(2,0,obs->low,obs->high,obsB0K0mumu,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b)/(obs->high-obs->low);
		
			binlowB0K0mumu=obs->low;
			binhighB0K0mumu=obs->high;
		}
		
		if(!strcmp(obs->type,"dGamma/dq2")) return obsB0K0mumu[0];
		else if(!strcmp(obs->type,"BR")) return obsB0K0mumu[0]*(obs->high-obs->low);
		else if(!strcmp(obs->type,"AFB")) return obsB0K0mumu[1];
		else if(!strcmp(obs->type,"FH")) return obsB0K0mumu[2];
		/****** Adding non-physical quantities - 20/04/2025 ******/
		else if(!strcmp(obs->type,"fp"))   return obsB0K0mumu[Nobs_BKll_phys+1];
		else if(!strcmp(obs->type,"f0"))  return obsB0K0mumu[Nobs_BKll_phys+2];
		else if(!strcmp(obs->type,"fT"))  return obsB0K0mumu[Nobs_BKll_phys+3];
	}

	else if(!strcmp(obs->decay,"BKmumu"))
	{
		if(reset||(fabs(binlowBKmumu-obs->low)>1.e-5)||(fabs(binhighBKmumu-obs->high)>1.e-5))
		{
			obsBKmumu[0]=BRBKll(2,1,obs->low,obs->high,obsBKmumu,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b)/(obs->high-obs->low);
		
			binlowBKmumu=obs->low;
			binhighBKmumu=obs->high;
		}
		
		if(!strcmp(obs->type,"dGamma/dq2")) return obsBKmumu[0];
		else if(!strcmp(obs->type,"BR")) return obsBKmumu[0]*(obs->high-obs->low);
		else if(!strcmp(obs->type,"AFB")) return obsBKmumu[1];
		else if(!strcmp(obs->type,"FH")) return obsBKmumu[2];
		/****** Adding non-physical quantities - 20/04/2025 ******/
		else if(!strcmp(obs->type,"fp"))   return obsBKmumu[Nobs_BKll_phys+1];
		else if(!strcmp(obs->type,"f0"))  return obsBKmumu[Nobs_BKll_phys+2];
		else if(!strcmp(obs->type,"fT"))  return obsBKmumu[Nobs_BKll_phys+3];
	}
	
	else if(!strcmp(obs->decay,"B0K0ee"))
	{
		if(reset||(fabs(binlowB0K0ee-obs->low)>1.e-5)||(fabs(binhighB0K0ee-obs->high)>1.e-5))
		{
			obsB0K0ee[0]=BRBKll(1,0,obs->low,obs->high,obsB0K0ee,C0eb,C1eb,C2eb,CQ0eb,CQ1eb,Cpeb,CQpeb,param,mu_b)/(obs->high-obs->low);
		
			binlowB0K0ee=obs->low;
			binhighB0K0ee=obs->high;
		}
		
		if(!strcmp(obs->type,"dGamma/dq2")) return obsB0K0ee[0];
		else if(!strcmp(obs->type,"BR")) return obsB0K0ee[0]*(obs->high-obs->low);
		else if(!strcmp(obs->type,"AFB")) return obsB0K0ee[1];
		else if(!strcmp(obs->type,"FH")) return obsB0K0ee[2];
		/****** Adding non-physical quantities - 20/04/2025 ******/
		else if(!strcmp(obs->type,"fp"))   return obsB0K0ee[Nobs_BKll_phys+1];
		else if(!strcmp(obs->type,"f0"))  return obsB0K0ee[Nobs_BKll_phys+2];
		else if(!strcmp(obs->type,"fT"))  return obsB0K0ee[Nobs_BKll_phys+3];
	}

	else if(!strcmp(obs->decay,"BKee"))
	{
		if(reset||(fabs(binlowBKee-obs->low)>1.e-5)||(fabs(binhighBKee-obs->high)>1.e-5))
		{
			obsBKee[0]=BRBKll(1,1,obs->low,obs->high,obsBKee,C0eb,C1eb,C2eb,CQ0eb,CQ1eb,Cpeb,CQpeb,param,mu_b)/(obs->high-obs->low);
		
			binlowBKee=obs->low;
			binhighBKee=obs->high;
		}
		
		if(!strcmp(obs->type,"dGamma/dq2")) return obsBKee[0];
		else if(!strcmp(obs->type,"BR")) return obsBKee[0]*(obs->high-obs->low);
		else if(!strcmp(obs->type,"AFB")) return obsBKee[1];
		else if(!strcmp(obs->type,"FH")) return obsBKee[2];
		else if(!strcmp(obs->type,"fp"))   return obsBKee[Nobs_BKll_phys+1];
		else if(!strcmp(obs->type,"f0"))  return obsBKee[Nobs_BKll_phys+2];
		else if(!strcmp(obs->type,"fT"))  return obsBKee[Nobs_BKll_phys+3];
	}

	else if(!strcmp(obs->decay,"B0K0tautau"))
	{
		if(reset||(fabs(binlowB0K0tautau-obs->low)>1.e-5)||(fabs(binhighB0K0tautau-obs->high)>1.e-5))
		{
			obsB0K0tautau[0]=BRBKll(3,0,obs->low,obs->high,obsB0K0tautau,C0taub,C1taub,C2taub,CQ0taub,CQ1taub,Cptaub,CQptaub,param,mu_b)/(obs->high-obs->low);
		
			binlowB0K0tautau=obs->low;
			binhighB0K0tautau=obs->high;
		}
		
		if(!strcmp(obs->type,"dGamma/dq2")) return obsB0K0tautau[0];
		else if(!strcmp(obs->type,"BR")) return obsB0K0tautau[0]*(obs->high-obs->low);
		else if(!strcmp(obs->type,"AFB")) return obsB0K0tautau[1];
		else if(!strcmp(obs->type,"FH")) return obsB0K0tautau[2];
		else if(!strcmp(obs->type,"fp"))   return obsB0K0tautau[Nobs_BKll_phys+1];
		else if(!strcmp(obs->type,"f0"))  return obsB0K0tautau[Nobs_BKll_phys+2];
		else if(!strcmp(obs->type,"fT"))  return obsB0K0tautau[Nobs_BKll_phys+3];
	}

	else if(!strcmp(obs->decay,"BKtautau"))
	{
		if(reset||(fabs(binlowBKtautau-obs->low)>1.e-5)||(fabs(binhighBKtautau-obs->high)>1.e-5))
		{
			obsBKtautau[0]=BRBKll(3,1,obs->low,obs->high,obsBKtautau,C0taub,C1taub,C2taub,CQ0taub,CQ1taub,Cptaub,CQptaub,param,mu_b)/(obs->high-obs->low);
		
			binlowBKtautau=obs->low;
			binhighBKtautau=obs->high;
		}
		
		if(!strcmp(obs->type,"dGamma/dq2")) return obsBKtautau[0];
		else if(!strcmp(obs->type,"BR")) return obsBKtautau[0]*(obs->high-obs->low);
		else if(!strcmp(obs->type,"AFB")) return obsBKtautau[1];
		else if(!strcmp(obs->type,"FH")) return obsBKtautau[2];
		else if(!strcmp(obs->type,"fp"))   return obsBKtautau[Nobs_BKll_phys+1];
		else if(!strcmp(obs->type,"f0"))  return obsBKtautau[Nobs_BKll_phys+2];
		else if(!strcmp(obs->type,"fT"))  return obsBKtautau[Nobs_BKll_phys+3];
	}

	else if(!strcmp(obs->type,"R-1")&&!strcmp(obs->decay,"B0K0ll"))
	{
		binlowB0K0mumu=binlowB0K0ee=obs->low;
		binhighB0K0mumu=binhighB0K0ee=obs->high;
		
		obsB0K0mumu[0]=BRBKll(2,0,obs->low,obs->high,obsB0K0mumu,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b)/(obs->high-obs->low);
		obsB0K0ee[0]=BRBKll(1,0,obs->low,obs->high,obsB0K0ee,C0eb,C1eb,C2eb,CQ0eb,CQ1eb,Cpeb,CQpeb,param,mu_b)/(obs->high-obs->low);
		
		return obsB0K0mumu[0]/obsB0K0ee[0]-1.;
	}

	else if(!strcmp(obs->type,"R-1")&&!strcmp(obs->decay,"BKll"))
	{
		binlowBKmumu=binlowBKee=obs->low;
		binhighBKmumu=binhighBKee=obs->high;
		
		obsBKmumu[0]=BRBKll(2,1,obs->low,obs->high,obsBKmumu,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b)/(obs->high-obs->low);
		obsBKee[0]=BRBKll(1,1,obs->low,obs->high,obsBKee,C0eb,C1eb,C2eb,CQ0eb,CQ1eb,Cpeb,CQpeb,param,mu_b)/(obs->high-obs->low);

		return obsBKmumu[0]/obsBKee[0]-1;
	}

/* Bs -> phi l+ l- */

	else if(!strcmp(obs->decay,"Bsphimumu"))
	{
		if(reset||(fabs(binlowBsphimumu-obs->low)>1.e-5)||(fabs(binhighBsphimumu-obs->high)>1.e-5))
		{
			obsBsphimumu[0]=BRBsphill(2,obs->low,min(obs->high,pow(param->m_Bs-param->m_phi,2.)*0.999),obsBsphimumu,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b)/(min(obs->high,pow(param->m_Bs-param->m_phi,2.)*0.999)-obs->low);
				
			binlowBsphimumu=obs->low;
			binhighBsphimumu=obs->high;
		}
		
		if(!strcmp(obs->type,"dGamma/dq2")) return obsBsphimumu[0];
		else if(!strcmp(obs->type,"BR")) return obsBsphimumu[0]*(obs->high-obs->low);
		else if(!strcmp(obs->type,"FL")) return obsBsphimumu[1];
		else if(!strcmp(obs->type,"AT2")||!strcmp(obs->type,"P1")) return obsBsphimumu[2];
		else if(!strcmp(obs->type,"S3")) return obsBsphimumu[3];
		else if(!strcmp(obs->type,"S4")) return obsBsphimumu[4];
		else if(!strcmp(obs->type,"P4prime")) return obsBsphimumu[5];
		else if(!strcmp(obs->type,"S7")) return obsBsphimumu[6];
		else if(!strcmp(obs->type,"P6prime")) return obsBsphimumu[7];
		else if(!strcmp(obs->type,"ATReCP")) return obsBsphimumu[8];
		else if(!strcmp(obs->type,"ATImCP")) return obsBsphimumu[9];
		else if(!strcmp(obs->type,"A5")) return obsBsphimumu[10];
		else if(!strcmp(obs->type,"AFBCP")) return obsBsphimumu[11];
		else if(!strcmp(obs->type,"A8")) return obsBsphimumu[12];
		else if(!strcmp(obs->type,"A9")) return obsBsphimumu[13];
		else if(!strcmp(obs->type,"P2CP")) return obsBsphimumu[14];
		else if(!strcmp(obs->type,"P3CP")) return obsBsphimumu[15];
		else if(!strcmp(obs->type,"P5primeCP")) return obsBsphimumu[16];
		else if(!strcmp(obs->type,"P8primeCP")) return obsBsphimumu[17];
		else if(!strcmp(obs->type,"A6c")) return obsBsphimumu[18];
		else if(!strcmp(obs->type,"AFBCPfull")) return 3./8.*obsBsphimumu[18]+obsBsphimumu[11];
		else if(!strcmp(obs->type,"Q8m")) return obsBsphimumu[19];
		else if(!strcmp(obs->type,"Q8p")) return obsBsphimumu[20];
		else if(!strcmp(obs->type,"Q9")) return obsBsphimumu[21];
		else if(!strcmp(obs->type,"S2s")) return obsBsphimumu[22];
		else if(!strcmp(obs->type,"V"))   return obsBsphimumu[Nobs_Bsphill_phys+1];
		else if(!strcmp(obs->type,"A0"))  return obsBsphimumu[Nobs_Bsphill_phys+2];
		else if(!strcmp(obs->type,"A1"))  return obsBsphimumu[Nobs_Bsphill_phys+3];
		else if(!strcmp(obs->type,"A12")) return obsBsphimumu[Nobs_Bsphill_phys+4];
		else if(!strcmp(obs->type,"T1"))  return obsBsphimumu[Nobs_Bsphill_phys+5];
		else if(!strcmp(obs->type,"T2"))  return obsBsphimumu[Nobs_Bsphill_phys+6];
		else if(!strcmp(obs->type,"T23")) return obsBsphimumu[Nobs_Bsphill_phys+7];
		else if(!strcmp(obs->type,"A2"))  return obsBsphimumu[Nobs_Bsphill_phys+8];
		/*************************************************************/
	}

	else if(!strcmp(obs->decay,"Bsphiee"))
	{
		if(reset||(fabs(binlowBsphiee-obs->low)>1.e-5)||(fabs(binhighBsphiee-obs->high)>1.e-5))
		{
			obsBsphiee[0]=BRBsphill(1,obs->low,min(obs->high,pow(param->m_Bs-param->m_phi,2.)*0.999),obsBsphiee,C0eb,C1eb,C2eb,CQ0eb,CQ1eb,Cpeb,CQpeb,param,mu_b)/(min(obs->high,pow(param->m_Bs-param->m_phi,2.)*0.999)-obs->low);
		
			binlowBsphiee=obs->low;
			binhighBsphiee=obs->high;
		}
		
		if(!strcmp(obs->type,"dGamma/dq2")) return obsBsphiee[0];
		else if(!strcmp(obs->type,"BR")) return obsBsphiee[0]*(obs->high-obs->low);
		else if(!strcmp(obs->type,"FL")) return obsBsphiee[1];
		else if(!strcmp(obs->type,"AT2")||!strcmp(obs->type,"P1")) return obsBsphiee[2];
		else if(!strcmp(obs->type,"S3")) return obsBsphiee[3];
		else if(!strcmp(obs->type,"S4")) return obsBsphiee[4];
		else if(!strcmp(obs->type,"P4prime")) return obsBsphiee[5];
		else if(!strcmp(obs->type,"S7")) return obsBsphiee[6];
		/* SN: additional observables */
		else if(!strcmp(obs->type,"P6prime")) return obsBsphiee[7];
		else if(!strcmp(obs->type,"ATReCP")) return obsBsphiee[8];
		else if(!strcmp(obs->type,"ATImCP")) return obsBsphiee[9];
		else if(!strcmp(obs->type,"A5")) return obsBsphiee[10];
		else if(!strcmp(obs->type,"AFBCP")) return obsBsphiee[11];
		else if(!strcmp(obs->type,"A8")) return obsBsphiee[12];
		else if(!strcmp(obs->type,"A9")) return obsBsphiee[13];
		else if(!strcmp(obs->type,"P2CP")) return obsBsphiee[14];
		else if(!strcmp(obs->type,"P3CP")) return obsBsphiee[15];
		else if(!strcmp(obs->type,"P5primeCP")) return obsBsphiee[16];
		else if(!strcmp(obs->type,"P8primeCP")) return obsBsphiee[17];
		else if(!strcmp(obs->type,"A6c")) return obsBsphiee[18];
		else if(!strcmp(obs->type,"AFBCPfull")) return 3./8.*obsBsphiee[18]+obsBsphiee[11];
		else if(!strcmp(obs->type,"Q8m")) return obsBsphiee[19];
		else if(!strcmp(obs->type,"Q8p")) return obsBsphiee[20];
		else if(!strcmp(obs->type,"Q9")) return obsBsphiee[21];
		else if(!strcmp(obs->type,"S2s")) return obsBsphiee[22];
		else if(!strcmp(obs->type,"V"))   return obsBsphiee[Nobs_Bsphill_phys+1];
		else if(!strcmp(obs->type,"A0"))  return obsBsphiee[Nobs_Bsphill_phys+2];
		else if(!strcmp(obs->type,"A1"))  return obsBsphiee[Nobs_Bsphill_phys+3];
		else if(!strcmp(obs->type,"A12")) return obsBsphiee[Nobs_Bsphill_phys+4];
		else if(!strcmp(obs->type,"T1"))  return obsBsphiee[Nobs_Bsphill_phys+5];
		else if(!strcmp(obs->type,"T2"))  return obsBsphiee[Nobs_Bsphill_phys+6];
		else if(!strcmp(obs->type,"T23")) return obsBsphiee[Nobs_Bsphill_phys+7];
		else if(!strcmp(obs->type,"A2"))  return obsBsphiee[Nobs_Bsphill_phys+8];
		/*************************************************************/
	}

	else if(!strcmp(obs->decay,"Bsphitautau"))
	{
		if(reset||(fabs(binlowBsphitautau-obs->low)>1.e-5)||(fabs(binhighBsphitautau-obs->high)>1.e-5))
		{
			obsBsphitautau[0]=BRBsphill(3,obs->low,min(obs->high,pow(param->m_Bs-param->m_phi,2.)*0.999),obsBsphitautau,C0taub,C1taub,C2taub,CQ0taub,CQ1taub,Cptaub,CQptaub,param,mu_b)/(min(obs->high,pow(param->m_Bs-param->m_phi,2.)*0.999)-obs->low);
		
			binlowBsphitautau=obs->low;
			binhighBsphitautau=obs->high;
		}
		
		if(!strcmp(obs->type,"dGamma/dq2")) return obsBsphitautau[0];
		else if(!strcmp(obs->type,"BR")) return obsBsphitautau[0]*(obs->high-obs->low);
		else if(!strcmp(obs->type,"FL")) return obsBsphitautau[1];
		else if(!strcmp(obs->type,"AT2")||!strcmp(obs->type,"P1")) return obsBsphitautau[2];
		else if(!strcmp(obs->type,"S3")) return obsBsphitautau[3];
		else if(!strcmp(obs->type,"S4")) return obsBsphitautau[4];
		else if(!strcmp(obs->type,"P4prime")) return obsBsphitautau[5];
		else if(!strcmp(obs->type,"S7")) return obsBsphitautau[6];
		else if(!strcmp(obs->type,"P6prime")) return obsBsphitautau[7];
		else if(!strcmp(obs->type,"ATReCP")) return obsBsphitautau[8];
		else if(!strcmp(obs->type,"ATImCP")) return obsBsphitautau[9];
		else if(!strcmp(obs->type,"A5")) return obsBsphitautau[10];
		else if(!strcmp(obs->type,"AFBCP")) return obsBsphitautau[11];
		else if(!strcmp(obs->type,"A8")) return obsBsphitautau[12];
		else if(!strcmp(obs->type,"A9")) return obsBsphitautau[13];
		else if(!strcmp(obs->type,"P2CP")) return obsBsphitautau[14];
		else if(!strcmp(obs->type,"P3CP")) return obsBsphitautau[15];
		else if(!strcmp(obs->type,"P5primeCP")) return obsBsphitautau[16];
		else if(!strcmp(obs->type,"P8primeCP")) return obsBsphitautau[17];
		else if(!strcmp(obs->type,"A6c")) return obsBsphitautau[18];
		else if(!strcmp(obs->type,"AFBCPfull")) return 3./8.*obsBsphitautau[18]+obsBsphitautau[11];
		else if(!strcmp(obs->type,"Q8m")) return obsBsphitautau[19];
		else if(!strcmp(obs->type,"Q8p")) return obsBsphitautau[20];
		else if(!strcmp(obs->type,"Q9")) return obsBsphitautau[21];
		else if(!strcmp(obs->type,"S2s")) return obsBsphitautau[22];
		else if(!strcmp(obs->type,"V"))   return obsBsphitautau[Nobs_Bsphill_phys+1];
		else if(!strcmp(obs->type,"A0"))  return obsBsphitautau[Nobs_Bsphill_phys+2];
		else if(!strcmp(obs->type,"A1"))  return obsBsphitautau[Nobs_Bsphill_phys+3];
		else if(!strcmp(obs->type,"A12")) return obsBsphitautau[Nobs_Bsphill_phys+4];
		else if(!strcmp(obs->type,"T1"))  return obsBsphitautau[Nobs_Bsphill_phys+5];
		else if(!strcmp(obs->type,"T2"))  return obsBsphitautau[Nobs_Bsphill_phys+6];
		else if(!strcmp(obs->type,"T23")) return obsBsphitautau[Nobs_Bsphill_phys+7];
		else if(!strcmp(obs->type,"A2"))  return obsBsphitautau[Nobs_Bsphill_phys+8];
		/*************************************************************/
	}

	else if(!strcmp(obs->type,"R-1")&&!strcmp(obs->decay,"Bsphill"))
	{
		binlowBsphimumu=binlowBsphiee=obs->low;
		binhighBsphimumu=binhighBsphiee=obs->high;
		
		obsBsphimumu[0]=BRBsphill(2,obs->low,min(obs->high,pow(param->m_Bs-param->m_phi,2.)*0.999),obsBsphimumu,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b)/(obs->high-obs->low);
		obsBsphiee[0]=BRBsphill(1,obs->low,min(obs->high,pow(param->m_Bs-param->m_phi,2.)*0.999),obsBsphiee,C0eb,C1eb,C2eb,CQ0eb,CQ1eb,Cpeb,CQpeb,param,mu_b)/(obs->high-obs->low);
	
		return obsBsphimumu[0]/obsBsphiee[0]-1.;
	}

/* Lambda_b -> Lambda l+ l- */

	else if(!strcmp(obs->decay,"LambdabLambdamumu"))
	{
		if(reset||(fabs(binlowLbLmumu-obs->low)>1.e-5)||(fabs(binhighLbLmumu-obs->high)>1.e-5))
		{
			obsLbLmumu[0]=BRLbLll(2,obs->low,min(obs->high,pow(param->m_Lambdab-param->m_Lambda,2.)*0.999),obsLbLmumu,C0b,C1b,C2b,Cpb,param,mu_b)/(min(obs->high,pow(param->m_Lambdab-param->m_Lambda,2.)*0.999)-obs->low);
				
			binlowLbLmumu=obs->low;
			binhighLbLmumu=obs->high;
		}
		
		if(!strcmp(obs->type,"dGamma/dq2")) return obsLbLmumu[0];
		else if(!strcmp(obs->type,"BR")) return obsLbLmumu[0]*(obs->high-obs->low);
		else if(!strcmp(obs->type,"AlFB")) return obsLbLmumu[1];
		else if(!strcmp(obs->type,"AhFB")) return obsLbLmumu[2];
		else if(!strcmp(obs->type,"AlhFB")) return obsLbLmumu[3];
		else if(!strcmp(obs->type,"FL")) return obsLbLmumu[4];
		else if(!strcmp(obs->type,"FT")) return obsLbLmumu[5];
	}

	else if(!strcmp(obs->decay,"LambdabLambdaee"))
	{
		if(reset||(fabs(binlowLbLee-obs->low)>1.e-5)||(fabs(binhighLbLee-obs->high)>1.e-5))
		{
			obsLbLee[0]=BRLbLll(1,obs->low,min(obs->high,pow(param->m_Lambdab-param->m_Lambda,2.)*0.999),obsLbLee,C0eb,C1eb,C2eb,Cpeb,param,mu_b)/(min(obs->high,pow(param->m_Lambdab-param->m_Lambda,2.)*0.999)-obs->low);
				
			binlowLbLee=obs->low;
			binhighLbLee=obs->high;
		}
		
		if(!strcmp(obs->type,"dGamma/dq2")) return obsLbLee[0];
		else if(!strcmp(obs->type,"BR")) return obsLbLee[0]*(obs->high-obs->low);
		else if(!strcmp(obs->type,"AlFB")) return obsLbLee[1];
		else if(!strcmp(obs->type,"AhFB")) return obsLbLee[2];
		else if(!strcmp(obs->type,"AlhFB")) return obsLbLee[3];
		else if(!strcmp(obs->type,"FL")) return obsLbLee[4];
		else if(!strcmp(obs->type,"FT")) return obsLbLee[5];
	}

	else if(!strcmp(obs->decay,"LambdabLambdatautau"))
	{
		if(reset||(fabs(binlowLbLtautau-obs->low)>1.e-5)||(fabs(binhighLbLtautau-obs->high)>1.e-5))
		{
			obsLbLtautau[0]=BRLbLll(3,obs->low,min(obs->high,pow(param->m_Lambdab-param->m_Lambda,2.)*0.999),obsLbLtautau,C0taub,C1taub,C2taub,Cptaub,param,mu_b)/(min(obs->high,pow(param->m_Lambdab-param->m_Lambda,2.)*0.999)-obs->low);
				
			binlowLbLtautau=obs->low;
			binhighLbLtautau=obs->high;
		}
		
		if(!strcmp(obs->type,"dGamma/dq2")) return obsLbLtautau[0];
		else if(!strcmp(obs->type,"BR")) return obsLbLtautau[0]*(obs->high-obs->low);
		else if(!strcmp(obs->type,"AlFB")) return obsLbLtautau[1];
		else if(!strcmp(obs->type,"AhFB")) return obsLbLtautau[2];
		else if(!strcmp(obs->type,"AlhFB")) return obsLbLtautau[3];
		else if(!strcmp(obs->type,"FL")) return obsLbLtautau[4];
		else if(!strcmp(obs->type,"FT")) return obsLbLtautau[5];
	}

	else if(!strcmp(obs->type,"R-1")&&!strcmp(obs->decay,"LambdabLambdall"))
	{
		binlowLbLmumu=binlowLbLee=obs->low;
		binhighLbLmumu=binhighLbLee=obs->high;
		
		obsLbLmumu[0]=BRLbLll(2,obs->low,min(obs->high,pow(param->m_Lambdab-param->m_Lambda,2.)*0.999),obsLbLmumu,C0b,C1b,C2b,Cpb,param,mu_b)/(obs->high-obs->low);
		obsLbLee[0]=BRLbLll(1,obs->low,min(obs->high,pow(param->m_Lambdab-param->m_Lambda,2.)*0.999),obsLbLee,C0eb,C1eb,C2eb,Cpeb,param,mu_b)/(obs->high-obs->low);
	
		return obsLbLmumu[0]/obsLbLee[0]-1.;
	}
	
	else if(!strcmp(obs->decay,"Kpinunu"))
	{
		if(!strcmp(obs->type,"BR")) return BRKppipnunu(0.,0.,0.,0.,0.,0.,param);
	}

	else if(!strcmp(obs->decay,"KLpi0nunu"))
	{
		if(!strcmp(obs->type,"BR")) return BRKLpinunu(0.,0.,0.,0.,0.,0.,param);
	}

	else if(!strcmp(obs->decay,"KLmumu"))
	{
		if(!strcmp(obs->type,"BR")) return BRKLll(2,0,0.,0.,0.,0.,0.,0.,param);
	}

	else if(!strcmp(obs->decay,"KLee"))
	{
		if(!strcmp(obs->type,"BR")) return BRKLll(1,0,0.,0.,0.,0.,0.,0.,param);
	}

	else if(!strcmp(obs->decay,"KSmumu"))
	{
		if(!strcmp(obs->type,"BR")) return BRKSmumu(0,0.,0.,0.,0.,0.,0.,param);
	}

	else if(!strcmp(obs->decay,"KLpi0mumu"))
	{
		if(!strcmp(obs->type,"BR")) return BRKLpill(2,0.,0.,0.,0.,0.,0.,0.,0.,param);
		else if(!strcmp(obs->type,"AFB")) return AFBmu_KLpill(0.,0.,0.,0.,0.,0.,0.,0.,param);
	}

	else if(!strcmp(obs->decay,"KLpi0ee"))
	{
		if(!strcmp(obs->type,"BR")) return BRKLpill(1,0.,0.,0.,0.,0.,0.,0.,0.,param);
	}

	return -1.e30;
}

/*---------------------------------------------------------------------*/

/*---------------------------------------------------------------------*/

void get_predictions(char names[][50], int *nbobs, double** predictions, struct parameters* param)
/* Computes the predictions/reference values of observables in string array "names" */
{		
	struct nuisance nuisparam;
	read_nuisance(param->nuisance_values,&nuisparam);
		
	get_predictions_nuisance(names,nbobs,predictions,param,&nuisparam);

	return;
}

/*---------------------------------------------------------------------*/

void get_th_covariance(double ***covariance_th, char names[][50], int *nbobs, struct parameters* param)
/* Computes the theory covariance matrix for the "nbobs" observables with names given in "names", using the nuisance parameters contained in the file "param->nuisance_values" */
{				
	int n=NNUIS; /* number of nuisance parameters */

	struct nuisance nuisparam;
	read_nuisance(param->nuisance_values,&nuisparam);
	
	char namenuisance[NNUIS+1][50];
	observables(0,NULL,0,NULL,NULL,&nuisparam,namenuisance,param);

	int ie,je;
	double **corr; /* correlation matrix between the parameters */
	corr=(double **) malloc((n+1)*sizeof(double *));
	for(ie=0;ie<=n;ie++) corr[ie]=(double *) malloc((n+1)*sizeof(double));
	for(ie=0;ie<=n;ie++) for(je=0;je<=n;je++) corr[ie][je]=(ie==je);
	
	read_correlation(param->nuisance_corr,corr,namenuisance,n);

	get_th_covariance_nuisance(covariance_th,names,nbobs,param,&nuisparam,corr);

	for (ie = 0; ie < n+1; ie++) free(corr[ie]);
	free(corr);

	return;
}

