#include "src/include.h"
/*--------------------------------------------------------*/
/* Calculation of the observables in the Standard Model   */
/*--------------------------------------------------------*/

int main()
{ 
	struct parameters param;
		
	Init_param(&param);
	
	slha_adjust(&param);
	param.SM=1;
	
	printf("\n");
	
	printf("SuperIso " VERSION " - %d\n\n", YEAR);
	printf("Standard Model predictions\n\n");
	printf("Observable\t\t\tValue\n\n");

	double complex C0b[11],C0spec[11],C1b[11],C1spec[11],C0w[11],C1w[11],C2w[11],C2b[11],Cpb[11],CQpb[3],CQ0b[3],CQ1b[3];
	CQ0b[1]=CQ0b[2]=CQ1b[1]=CQ1b[2]=CQpb[1]=CQpb[2]=0.;
	double obs[Nobs_BKsll+1];

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole/2.;
	CW_calculator(2,C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	Cprime_calculator(2,Cpb,CQpb,mu_W,mu_b,&param);
	printf("BR(b->s gamma)\t\t\t%.3e\n",bsgamma(C0b,C1b,C2b,Cpb,mu_b,mu_W,&param));
	
	double lambda_h=0.5;
	double mu_spec=sqrt(lambda_h*param.mass_b);		
	C_calculator_base2(C0w,C1w,mu_W,C0b,C1b,mu_b,&param);
	C_calculator_base2(C0w,C1w,mu_W,C0spec,C1spec,mu_spec,&param);
	Cprime_calculator(2,Cpb,CQpb,mu_W,mu_b,&param);
	printf("delta0(B->K* gamma)\t\t%.3e\n\n",delta0(C0b,C0spec,C1b,C1spec,Cpb,&param,mu_b,mu_spec,lambda_h));
	
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	Cprime_calculator(2,Cpb,CQpb,mu_W,mu_b,&param);
	
	printf("BR(B0->K* gamma)\t\t%.3e\n",BR_BKstargamma(0,C0b,C1b,C2b,Cpb,&param,mu_b));
	printf("BR(B+->K* gamma)\t\t%.3e\n\n",BR_BKstargamma(1,C0b,C1b,C2b,Cpb,&param,mu_b));
	
	printf("BR(Bs->mu mu)_untag\t\t%.3e\n",Bsmumu_untag(C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b));
	printf("BR(Bd->mu mu)\t\t\t%.3e\n\n",Bdmumu(C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b));

 	printf("BR(B->tau nu)\t\t\t%.3e\n",Btaunu(&param));
    printf("R(B->tau nu)\t\t\t%.3e\n\n",RBtaunu(&param));
    printf("BR(B->D tau nu)\t\t\t%.3e\n",BRBDlnu_full(3,1,obs,&param));
    printf("AFB(B->D tau nu)\t\t%.3e\n",obs[1]);
    printf("Ptau(B->D tau nu)\t\t%.3e\n",obs[2]);
	printf("BR(B->D tau nu)/BR(B->D e nu)\t%.3e\n\n",BDtaunu_BDenu(&param));
    printf("BR(B->D* tau nu)\t\t%.3e\n",BRBDstarlnu_full(3,1,obs,&param));
    printf("AFB(B->D* tau nu)\t\t%.3e\n",obs[1]);
    printf("Ptau(B->D* tau nu)\t\t%.3e\n",obs[2]);
    printf("PD*(B->D* tau nu)\t\t%.3e\n",obs[3]);
	printf("BR(B->D* tau nu)/BR(B->D* e nu)\t%.3e\n\n",BDstartaunu_BDstarenu(&param));
	printf("BR(Ds->tau nu)\t\t\t%.3e\n",Dstaunu(&param));
	printf("BR(Ds->mu nu)\t\t\t%.3e\n",Dsmunu(&param));
	printf("BR(D->mu nu)\t\t\t%.3e\n",Dmunu(&param));
	printf("BR(K->mu nu)/BR(pi->mu nu)\t%.3e\n",Kmunu_pimunu(&param));
	printf("Rmu23(K->mu nu)\t\t\t%.3e\n\n",Rmu23(&param));


	double complex CLe, CRe, CLmu, CRmu, CLtau, CRtau;
	CLe = CRe = CLmu = CRmu = CLtau = CRtau = 0.;

    /*** K->pinunu ***/
	printf("Br(K+->pi+ nu nu):\t\t%.2e\n",BRKppipnunu(CLe, CRe, CLmu, CRmu, CLtau, CRtau, &param));
	printf("Br(KL->pi0 nu nu):\t\t%.2e\n\n",BRKLpinunu(CLe, CRe, CLmu, CRmu, CLtau, CRtau, &param));

    /*** K->ll ***/
	int SL;
	SL = 0; //SL=0 => full, SL=1 => only Short Distance, SL=2 => only Long Distance
	double complex C10, C10p, CQ1, CQ1p, CQ2, CQ2p;
	C10 = C10p = CQ1 = CQ1p = CQ2 = CQ2p = 0.;
	/* sign of LD contribution for KL->mumu is unknown (see eq.2.16 1711.11030).
	 * In leshouches.c the default value is positive (param.LDsign_KLmumu = 1)*/
	param.LDsign_KLmumu = 1;
	printf("Br(KL->mu mu)[full;LD%+d]:\t%.2e\n",param.LDsign_KLmumu, BRKLll(2, SL, C10, C10p, CQ1, CQ1p, CQ2, CQ2p, &param));
	param.LDsign_KLmumu = -1; //changing the sign of long-distance contribution to negative
	printf("Br(KL->mu mu)[full;LD%+d]:\t%.2e\n\n",param.LDsign_KLmumu, BRKLll(2,SL, C10, C10p, CQ1, CQ1p, CQ2, CQ2p, &param));
	param.LDsign_KLmumu = 1;
	printf("Br(KL->e e)[full;LD%+d]:\t\t%.2e\n",param.LDsign_KLmumu, BRKLll(1, SL, C10, C10p, CQ1, CQ1p, CQ2, CQ2p, &param));
	param.LDsign_KLmumu = -1; //changing the sign of long-distance contribution to negative
	printf("Br(KL->e e)[full;LD%+d]:\t\t%.2e\n\n",param.LDsign_KLmumu, BRKLll(1,SL, C10, C10p, CQ1, CQ1p, CQ2, CQ2p, &param));
	printf("Br(KS->mu mu)_full:\t\t%.2e\n\n",BRKSmumu(SL, C10, C10p, CQ1, CQ1p, CQ2, CQ2p, &param));

    /*** KL -> pi l+ l- ***/
	double complex C9,C9p;
	C9=C9p=C10=C10p=0.;
	/* Constructive interference is theoretically preferred and by default assumed in SuperIso*/
	param.sign_KLpill = 1; //by default is +1 (avoiding confusion in case in leshouches.c it has been changed)
	printf("Br(KL->pi0ee)_const-interf.:\t%.2e\n",BRKLpill(1, C9, C9p, C10, C10p, CQ1, CQ1p, CQ2, CQ2p, &param));

	param.sign_KLpill = -1; //Changing to destructive interference
	printf("Br(KL->pi0ee)_dest-interf.:\t%.2e\n\n",BRKLpill(1, C9, C9p, C10, C10p, CQ1, CQ1p, CQ2, CQ2p, &param));

	param.sign_KLpill = 1; //Changing to constructive interference
	printf("Br(KL->pi0mumu)_const-interf.:\t%.2e\n",BRKLpill(2, C9, C9p, C10, C10p, CQ1, CQ1p, CQ2, CQ2p, &param));

	param.sign_KLpill = -1; //Changing to destructive interference: -1
	printf("Br(KL->pi0mumu)_dest-interf.:\t%.2e\n\n",BRKLpill(2, C9, C9p, C10, C10p, CQ1, CQ1p, CQ2, CQ2p, &param));




	return 1;
}
