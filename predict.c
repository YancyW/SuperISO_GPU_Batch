#include "src/include.h"

/*------------------------------------------------------------------------*/
/* Calculation of prediction +/- error, with given NP Wilson coefficients */
/*------------------------------------------------------------------------*/

/* Change the number of digits displayed with PRECISION below  */
#define PRECISION 5
#define MAX_OBS_NAME_LEN 50

void display_instructions();
int check_nonphys_obs(obsname* obs);
int check_non_bTOs_decay(obsname* obs);
const char* get_coefficient_name(const char* type, int index);
void parse_line(char* line, struct parameters* param, int* has_np_contributions);
int get_exponent(double number);
void format_with_matched_exponent(double main_value, double uncertainty, char *formatted_str, int precision);
void pretty_print(char *str, double a, double b, int precision, int padding_size);
void save_predictions_to_csv(const char* filename, char names[][MAX_OBS_NAME_LEN], double* predictions, double** covariance_th, int nobs);


/****************************************************************/
void display_instructions() {
    printf("\n");
    printf("***************************************************************************************\n");
    printf("*************** Theory Prediction for Flavour Observables *****************************\n");
    printf("***************************************************************************************\n");
    printf("*  \033[1;36mPrediction Information:\033[0m                                                            *\n");
    printf("*  - The observable prediction is displayed with its uncertainty (±)                  *\n");
    printf("*  - Predictions can account for New Physics effects via inputted                     *\n");
    printf("*    Wilson coefficients at the hadronic scale:                                       *\n");
    printf("*     + SM-like: C7, C8, C9, C10 & primed version (C7p, C8p, C9p, C10p)               *\n");
    printf("*     + Scalar & pseudoscalar: CQ1, CQ2 & primed versions (CQ1p,CQ2p)                 *\n");
    printf("*     + Optional: Lepton Flavour Universality Violating (LFUV)                        *\n");
    printf("*  - For the basis, see the manual at:                                                *\n");
    printf("*    \033[4;34mhttp://superiso.in2p3.fr/#manual\033[0m                                                 *\n");
    printf("*                                                                                     *\n");
    printf("* \033[1;36mAvailable Observables:\033[0m                                                              *\n");
    printf("*  - Enter one observable via terminal or (optionally) several via file               *\n");
    printf("*  - Observable should be in the form: observable_decay[_q2min_q2max]                 *\n");
    printf("*  - [q2-dependent] observables should include range (q2min to q2max)                 *\n");
    printf("* Example: \033[1;32mP5prime_B0Kstar0mumu_1_6.1\033[0m                                                 *\n");
    printf("* For more examples, see the list at chi2_input/myobs.in                              *\n");
    printf("*                                                                                     *\n");
    printf("* \033[1;36mProgram Execution Options:\033[0m                                                          *\n");
    printf("* 1. Interactive mode (no arguments):                                                 *\n");
    printf("*    \033[1;32m./predict.x\033[0m                                                                      *\n");
    printf("*    - Prompts for observable and real-valued Wilson coefficients                     *\n");
    printf("*    - To include complex Wilson coeffs, pass chi2_input/mycoeffs.in as argument      *\n");
    printf("*                                                                                     *\n");
    printf("* 2. Batch mode with arguments:                                                       *\n");
    printf("*    \033[1;32m./predict.x <obsfile> [outputname] [wilsonfile] [nuisfile] [nuiscorrfile]\033[0m        *\n");
    printf("*    - obsfile: list of observables (e.g., chi2_input/myobs.in)                       *\n");
    printf("*    - outputname: base name for output CSV file (predictions saved if provided)      *\n");
    printf("*    - wilsonfile: NP coefficients (e.g., chi2_input/mycoeffs.in)                     *\n");
    printf("*    - nuisfile: inputs and errors (default: chi2_input/nuisance.in)                  *\n");
    printf("*    - corrfile: error correlations (default: chi2_input/nuisance_corr.in)            *\n");
    printf("*                                                                                     *\n");
    printf("* Example 1: Prediction of observable list with default nuisance parameters           *\n");
    printf("*           \033[1;32m./predict.x chi2_input/myobs.in results\033[0m                                   *\n");
    printf("* Example 2: Prediction with custom Wilson coefficients and nuisance files            *\n");
    printf("*           \033[1;32m./predict.x chi2_input/myobs.in results chi2_input/mycoeffs.in \033[0m           *\n");
    printf("*           \033[1;32m      custom_nuis.in custom_corr.in\033[0m                                       *\n");
    printf("***************************************************************************************\n\n");
}

#include <sys/stat.h>  // Needed for mkdir()
#include <sys/types.h> // Needed for stat()


int main(int argc,char** argv)
{
    int has_np_contributions = 0;  // Flag to track non-zero contributions
    struct parameters param;
    int ie = 0;
    int nobs;
    obsname obs;

    char obsfile[500], outputname[500], wilsonfile[500], nuisfile[500], nuiscorrfile[500];
    char names[NBOBSMAX][MAX_OBS_NAME_LEN];

    // Initialize parameters
    Init_param(&param);
    slha_adjust(&param);
    param.SM=1;

    // Set default file names
    sprintf(nuisfile,"chi2_input/nuisance.in");
    sprintf(nuiscorrfile,"chi2_input/nuisance_corr.in");
    outputname[0] = '\0';  // Initialize outputname as empty

    // Parse command line arguments
    if(argc >= 2) {
        sscanf(argv[1],"%s",obsfile);
    }
    if(argc >= 3) {
        sscanf(argv[2],"%s",outputname);
    }
    if(argc >= 4) {
        sscanf(argv[3],"%s",wilsonfile);
    }
    if(argc >= 5) {
        sscanf(argv[4],"%s",nuisfile);
    }
    if(argc >= 6) {
        sscanf(argv[5],"%s",nuiscorrfile);
    }

    sprintf(param.nuisance_values, "%s", nuisfile);
    sprintf(param.nuisance_corr, "%s", nuiscorrfile);

    if(argc >= 5) printf("Reading inputs and uncertainties from %s\n",nuisfile);
    if(argc >= 6) printf("Reading uncertainty correlations from %s\n\n",nuiscorrfile);

    if(argc >= 2) {
        printf("Reading observable list from %s\n",obsfile);
        nobs=read_obs_list(obsfile,names);
    }

    // Handle Wilson coefficients - interactive mode if wilsonfile not provided
    if(argc >= 4) {
        // Read from file if wilsonfile is provided
        char line[256];
        FILE* file = fopen(wilsonfile, "r");
        if (!file) {
            perror("Error opening file");
            exit(1);
        }

        printf("Reading Wilson coefficients from %s\n",wilsonfile);

        while (fgets(line, sizeof(line), file)) {
            parse_line(line, &param, &has_np_contributions);
        }

        fclose(file);

        if (has_np_contributions) {
            printf("Predictions are given assuming above New Physics contributions.\n");
        } else {
            printf("No New Physics contributions considered, SM predictions are:\n");
        }
    } else if (argc < 2) {
        // Interactive mode for single observable
        char inputted_obs[MAX_OBS_NAME_LEN];
        char show_instructions;

        printf("Do you want to see the instructions? (y/n): ");
        scanf(" %c", &show_instructions);

        if (show_instructions == 'y' || show_instructions == 'Y') {
            display_instructions();
        }

        printf("\033[1;32mEnter the observable: \033[0m");
        scanf("%49s", inputted_obs);

        if (read_nameobs(inputted_obs, &obs) == 0) {
            printf("Failed to read observable name.\n");
            return 1;
        }
        if (check_nameobs(&obs) == 0) {
            printf("Invalid observable.\n");
            return 1;
        }

        nobs = 1;
        strcat(names[0], inputted_obs);
    }

    // Interactive Wilson coefficient input if no wilsonfile provided, obs physical and  b->s
    if (argc < 4 && !(check_nonphys_obs(&obs)) && !(check_non_bTOs_decay(&obs))) {
        char response_NP;
        printf("Do you want to consider New Physics contributions? (y/n): ");
        scanf(" %c", &response_NP);

		if (response_NP == 'y' || response_NP == 'Y') {
			printf("Do you want NP in radiative coefficients C7, C8 (y/n): ");
			char response_radiative;
			scanf(" %c", &response_radiative);
			if (response_radiative == 'y' || response_radiative == 'Y') {
						double dC7, dC8;
						printf("Enter the values for dC7, dC8 (separated by commas): ");
						scanf("%lf, %lf", &dC7, &dC8);
						param.deltaC[7] = dC7;
						param.deltaC[8] = dC8;
			}

			// Ask whether LFUV Wilson coefficients wanted
			char response_lfuv;
			printf("Do you want to input LFUV Wilson coefficients? (y/n): ");
			scanf(" %c", &response_lfuv);  // Space before %c to consume any leftover newline characters


			// Ask about SM-type Wilson coefficients
			char response_lh;
			printf("Do you want to input SM-like Wilson coefficients C9, C10? (y/n): ");
			scanf(" %c", &response_lh);

			if (response_lh == 'y' || response_lh == 'Y') {

				if (response_lfuv == 'y' || response_lfuv == 'Y') {
					// LFUV Wilson Coefficients
					char response_lepton;
					printf("Do you want non-zero NP in the electron mode? (y/n): ");
					scanf(" %c", &response_lepton);
					if (response_lepton == 'y' || response_lepton == 'Y') {
						double dC9e, dC10e;
						printf("Enter the values for dC9e, dC10e (separated by commas): ");
						scanf("%lf, %lf", &dC9e, &dC10e);

						param.deltaC[10+9] = dC9e;
						param.deltaC[10+10] = dC10e;
					}

					printf("Do you want non-zero NP in the muon mode? (y/n): ");
					scanf(" %c", &response_lepton);
					if (response_lepton == 'y' || response_lepton == 'Y') {
						double dC9mu, dC10mu;
						printf("Enter the values for dC9mu, dC10mu (separated by commas): ");
						scanf("%lf, %lf", &dC9mu, &dC10mu);

						param.deltaC[9] = dC9mu;
						param.deltaC[10] = dC10mu;
					}

					printf("Do you want non-zero NP in the tau mode? (y/n): ");
					scanf(" %c", &response_lepton);
					if (response_lepton == 'y' || response_lepton == 'Y') {
						double dC9tau, dC10tau;
						printf("Enter the values for dC9tau, dC10tau (separated by commas): ");
						scanf("%lf, %lf", &dC9tau, &dC10tau);

						param.deltaC[20+9] = dC9tau;
						param.deltaC[20+10] = dC10tau;
					}
				} else {
					// Flavor-universal Wilson coefficients
					double dC9, dC10;
					printf("Enter the values for dC9, dC10 (separated by commas): ");
					scanf("%lf, %lf", &dC9, &dC10);

					param.deltaC[9] = dC9;
					param.deltaC[10] = dC10;
					param.deltaC[10+9] = dC9;
					param.deltaC[10+10] = dC10;
					param.deltaC[20+9] = dC9;
					param.deltaC[20+10] = dC10;
				}
			}

			// Ask about right-handed Wilson coefficients
			char response_rh;
			printf("Do you want to input primed Wilson coefficients C7', C8', C9', C10'? (y/n): ");
			scanf(" %c", &response_rh);

			if (response_rh == 'y' || response_rh == 'Y') {
				printf("Do you want NP in primed radiative coefficients C7', C8' (y/n): ");
				char response_radiative_prime;
				scanf(" %c", &response_radiative_prime);
				if (response_radiative_prime == 'y' || response_radiative_prime == 'Y') {
							double dC7p, dC8p;
							printf("Enter the values for dC7, dC8 (separated by commas): ");
							scanf("%lf, %lf", &dC7p, &dC8p);
							param.deltaC[7] = dC7p;
							param.deltaC[8] = dC8p;
				}
				if (response_lfuv == 'y' || response_lfuv == 'Y') {
					// LFUV right-handed Wilson Coefficients
					char response_lepton;

					// Electron mode
					printf("Do you want non-zero NP in the electron mode (right-handed)? (y/n): ");
					scanf(" %c", &response_lepton);
					if (response_lepton == 'y' || response_lepton == 'Y') {
						double dC9p_e, dC10p_e;
						printf("Enter the values for dC9p_e, dC10p_e (separated by commas): ");
						scanf("%lf, %lf", &dC9p_e, &dC10p_e);

						param.deltaCp[10 + 9] = dC9p_e;
						param.deltaCp[10 + 10] = dC10p_e;
					}

					// Muon mode
					printf("Do you want non-zero NP in the muon mode (right-handed)? (y/n): ");
					scanf(" %c", &response_lepton);
					if (response_lepton == 'y' || response_lepton == 'Y') {
						double dC9p_mu, dC10p_mu;
						printf("Enter the values for dC9p_mu, dC10p_mu (separated by commas): ");
						scanf("%lf, %lf", &dC9p_mu, &dC10p_mu);

						param.deltaCp[9] = dC9p_mu;
						param.deltaCp[10] = dC10p_mu;
					}

					// Tau mode
					printf("Do you want non-zero NP in the tau mode (right-handed)? (y/n): ");
					scanf(" %c", &response_lepton);
					if (response_lepton == 'y' || response_lepton == 'Y') {
						double dC9p_tau, dC10p_tau;
						printf("Enter the values for dC9p_tau, dC10p_tau (separated by commas): ");
						scanf("%lf, %lf", &dC9p_tau, &dC10p_tau);

						param.deltaCp[20 + 9] = dC9p_tau;
						param.deltaCp[20 + 10] = dC10p_tau;
					}
				} else {
					// Flavor-universal right-handed Wilson coefficients
					double dC9p, dC10p;
					printf("Enter the values for dC9p, dC10p (separated by commas): ");
					scanf("%lf, %lf", &dC9p, &dC10p);

					param.deltaCp[9] = dC9p;
					param.deltaCp[10] = dC10p;
					param.deltaCp[10 + 9] = dC9p;
					param.deltaCp[10 + 10] = dC10p;
					param.deltaCp[20 + 9] = dC9p;
					param.deltaCp[20 + 10] = dC10p;
				}
			}

			// Ask about (pseudo-)scalar Wilson coefficients
			char response_scalar;
			printf("Do you want to consider scalar Wilson coefficients? (y/n): ");
			scanf(" %c", &response_scalar);  // Space before %c to consume any leftover newline characters
			if (response_scalar == 'y' || response_scalar == 'Y') {
				// LFUV scalar Wilson Coefficients
				if (response_lfuv == 'y' || response_lfuv == 'Y') {

					char response_scalar_mode;

					// Electron
					printf("Do you want non-zero NP in the electron mode? (y/n): ");
					scanf(" %c", &response_scalar_mode);
					if (response_scalar_mode == 'y' || response_scalar_mode == 'Y') {
					double dCQ1e, dCQ2e;
					printf("Enter the values for dCQ1e, dCQ2e (separated by commas): ");
					scanf("%lf, %lf", &dCQ1e, &dCQ2e);

					param.deltaCQ[1 + 2] = dCQ1e;
					param.deltaCQ[2 + 2] = dCQ2e;
					}

					// Muon
					printf("Do you want non-zero NP in the muon mode? (y/n): ");
					scanf(" %c", &response_scalar_mode);
					if (response_scalar_mode == 'y' || response_scalar_mode == 'Y') {
					double dCQ1mu, dCQ2mu;
					printf("Enter the values for dCQ1mu, dCQ2mu (separated by commas): ");
					scanf("%lf, %lf", &dCQ1mu, &dCQ2mu);

					param.deltaCQ[1] = dCQ1mu;
					param.deltaCQ[2] = dCQ2mu;
					}

					// Tau
					printf("Do you want non-zero NP in the tau mode? (y/n): ");
					scanf(" %c", &response_scalar_mode);
					if (response_scalar_mode == 'y' || response_scalar_mode == 'Y') {
					double dCQ1tau, dCQ2tau;
					printf("Enter the values for dCQ1tau, dCQ2tau (separated by commas): ");
					scanf("%lf, %lf", &dCQ1tau, &dCQ2tau);

					param.deltaCQ[1 + 4] = dCQ1tau;
					param.deltaCQ[2 + 4] = dCQ2tau;
					}
				} else {
					// Flavor-universal Wilson coefficients
					double dCQ1, dCQ2;
					printf("Enter the values for dCQ1, dCQ2 (separated by commas): ");
					scanf("%lf, %lf", &dCQ1, &dCQ2);

					param.deltaCQ[1] = dCQ1;
					param.deltaCQ[2] = dCQ2;
					param.deltaCQ[3] = dCQ1;
					param.deltaCQ[4] = dCQ2;
					param.deltaCQ[5] = dCQ1;
					param.deltaCQ[6] = dCQ2;
				}
			}


			// Ask about primed (pseudo-)scalar Wilson coefficients
			char response_scalar_rh;
			printf("Do you want to consider primed scalar Wilson coefficients? (y/n): ");
			scanf(" %c", &response_scalar_rh);  // Space before %c to consume any leftover newline characters
			if (response_scalar_rh == 'y' || response_scalar_rh == 'Y') {
				// LFUV Wilson primed scalar Coefficients
				if (response_lfuv == 'y' || response_lfuv == 'Y') {

					char response_scalar_mode;

					// Electron
					printf("Do you want non-zero NP in the electron mode? (y/n): ");
					scanf(" %c", &response_scalar_mode);
					if (response_scalar_mode == 'y' || response_scalar_mode == 'Y') {
					double dCQp1e, dCQp2e;
					printf("Enter the values for dCQp1e, dCQp2e (separated by commas): ");
					scanf("%lf, %lf", &dCQp1e, &dCQp2e);

					param.deltaCQp[1 + 2] = dCQp1e;
					param.deltaCQp[2 + 2] = dCQp2e;
					}

					// Muon
					printf("Do you want non-zero NP in the muon mode? (y/n): ");
					scanf(" %c", &response_scalar_mode);
					if (response_scalar_mode == 'y' || response_scalar_mode == 'Y') {
					double dCQp1mu, dCQp2mu;
					printf("Enter the values for dCQp1mu, dCQp2mu (separated by commas): ");
					scanf("%lf, %lf", &dCQp1mu, &dCQp2mu);

					param.deltaCQp[1] = dCQp1mu;
					param.deltaCQp[2] = dCQp2mu;
					}

					// Tau
					printf("Do you want non-zero NP in the tau mode? (y/n): ");
					scanf(" %c", &response_scalar_mode);
					if (response_scalar_mode == 'y' || response_scalar_mode == 'Y') {
					double dCQp1tau, dCQp2tau;
					printf("Enter the values for dCQp1tau, dCQp2tau (separated by commas): ");
					scanf("%lf, %lf", &dCQp1tau, &dCQp2tau);

					param.deltaCQp[1 + 4] = dCQp1tau;
					param.deltaCQp[2 + 4] = dCQp2tau;
					}
				} else {
					// Flavor-universal Wilson coefficients
					double dCQp1, dCQp2;
					printf("Enter the values for dCQp1, dCQp2 (separated by commas): ");
					scanf("%lf, %lf", &dCQp1, &dCQp2);

					param.deltaCQp[1] = dCQp1;
					param.deltaCQp[2] = dCQp2;
					param.deltaCQp[3] = dCQp1;
					param.deltaCQp[4] = dCQp2;
					param.deltaCQp[5] = dCQp1;
					param.deltaCQp[6] = dCQp2;
				}
			}
		}
		/* *************************************************************** */
    }

    // Interactive Wilson coefficient input if no wilsonfile provided, obs physical and  s->d
    if (argc < 4 && !(check_nonphys_obs(&obs)) && (check_non_bTOs_decay(&obs))) {
        char response_NP;
        printf("Do you want to consider New Physics contributions? (y/n): ");
        scanf(" %c", &response_NP);

		if (response_NP == 'y' || response_NP == 'Y') {
			// Ask whether LFUV Wilson coefficients wanted
			char response_lfuv;
			printf("Do you want to input LFUV Wilson coefficients? (y/n): ");
			scanf(" %c", &response_lfuv);  // Space before %c to consume any leftover newline characters


			// Ask about Wilson coefficients for neutrinos
			char response_CLR;
			printf("Do you want to input neutrino Wilson coefficient CL? (y/n): ");
			scanf(" %c", &response_CLR);

			if (response_CLR == 'y' || response_CLR == 'Y') {

				if (response_lfuv == 'y' || response_lfuv == 'Y') {
					// LFUV Wilson Coefficients
					char response_lepton;
					printf("Do you want non-zero NP in the electron mode? (y/n): ");
					scanf(" %c", &response_lepton);
					if (response_lepton == 'y' || response_lepton == 'Y') {
						double dCLe, dCRe;
						printf("Enter the value for dCLe, dCRe (separated by commas): ");
						scanf("%lf, %lf", &dCLe, &dCRe);

						param.KCLe = dCLe;
						param.KCRe = dCRe;
					}

					printf("Do you want non-zero NP in the muon mode? (y/n): ");
					scanf(" %c", &response_lepton);
					if (response_lepton == 'y' || response_lepton == 'Y') {
						double dCLmu, dCRmu;
						printf("Enter the value for dCLmu, dCRmu (separated by commas): ");
						scanf("%lf, %lf", &dCLmu, &dCRmu);

						param.KCLmu = dCLmu;
						param.KCRmu = dCRmu;
					}

					printf("Do you want non-zero NP in the tau mode? (y/n): ");
					scanf(" %c", &response_lepton);
					if (response_lepton == 'y' || response_lepton == 'Y') {
						double dCLtau, dCRtau;
						printf("Enter the value for dCLtau, dCRtau (separated by commas): ");
						scanf("%lf, %lf", &dCLtau, &dCRtau);

						param.KCLtau = dCLtau;
						param.KCRtau = dCRtau;
					}
				} else {
					// Flavor-universal Wilson coefficients
					double dCL, dCR;
					printf("Enter the values for dCL, dCR (separated by commas): ");
					scanf("%lf, %lf", &dCL, &dCR);

					param.KCLe = param.KCLmu = param.KCLtau = dCL;
					param.KCRe = param.KCRmu = param.KCRtau = dCR;
				}
			}
			// Ask about V-A Wilson coefficients
			char response_lh;
			printf("Do you want to input SM-like Wilson coefficients C9, C10? (y/n): ");
			scanf(" %c", &response_lh);

			if (response_lh == 'y' || response_lh == 'Y') {

				if (response_lfuv == 'y' || response_lfuv == 'Y') {
					// LFUV Wilson Coefficients
					char response_lepton;
					printf("Do you want non-zero NP in the electron mode? (y/n): ");
					scanf(" %c", &response_lepton);
					if (response_lepton == 'y' || response_lepton == 'Y') {
						double dC9e, dC10e;
						printf("Enter the values for dC9e, dC10e (separated by commas): ");
						scanf("%lf, %lf", &dC9e, &dC10e);

						param.KC9e = dC9e;
						param.KC10e = dC10e;
					}

					printf("Do you want non-zero NP in the muon mode? (y/n): ");
					scanf(" %c", &response_lepton);
					if (response_lepton == 'y' || response_lepton == 'Y') {
						double dC9mu, dC10mu;
						printf("Enter the values for dC9mu, dC10mu (separated by commas): ");
						scanf("%lf, %lf", &dC9mu, &dC10mu);

						param.KC9mu = dC9mu;
						param.KC10mu = dC10mu;
					}

				} else {
					// Flavor-universal Wilson coefficients
					double dC9, dC10;
					printf("Enter the values for dC9, dC10 (separated by commas): ");
					scanf("%lf, %lf", &dC9, &dC10);

					param.KC9e = param.KC9mu = dC9;
					param.KC10e = param.KC10mu = dC10;
				}
			}

			// Ask about right-handed Wilson coefficients
			char response_rh;
			printf("Do you want to input primed Wilson coefficients C9', C10'? (y/n): ");
			scanf(" %c", &response_rh);

			if (response_rh == 'y' || response_rh == 'Y') {
				if (response_lfuv == 'y' || response_lfuv == 'Y') {
					// LFUV right-handed Wilson Coefficients
					char response_lepton;

					// Electron mode
					printf("Do you want non-zero NP in the electron mode (right-handed)? (y/n): ");
					scanf(" %c", &response_lepton);
					if (response_lepton == 'y' || response_lepton == 'Y') {
						double dC9p_e, dC10p_e;
						printf("Enter the values for dC9p_e, dC10p_e (separated by commas): ");
						scanf("%lf, %lf", &dC9p_e, &dC10p_e);

						param.KC9pe = dC9p_e;
						param.KC10pe = dC10p_e;
					}

					// Muon mode
					printf("Do you want non-zero NP in the muon mode (right-handed)? (y/n): ");
					scanf(" %c", &response_lepton);
					if (response_lepton == 'y' || response_lepton == 'Y') {
						double dC9p_mu, dC10p_mu;
						printf("Enter the values for dC9p_mu, dC10p_mu (separated by commas): ");
						scanf("%lf, %lf", &dC9p_mu, &dC10p_mu);

						param.KC9pmu = dC9p_mu;
						param.KC10pmu = dC10p_mu;
					}

				} else {
					// Flavor-universal right-handed Wilson coefficients
					double dC9p, dC10p;
					printf("Enter the values for dC9p, dC10p (separated by commas): ");
					scanf("%lf, %lf", &dC9p, &dC10p);

					param.KC9pe = param.KC9pmu = dC9p;
					param.KC10pe = param.KC10pmu = dC10p;
				}
			}

			// Ask about (pseudo-)scalar Wilson coefficients
			char response_scalar;
			printf("Do you want to consider scalar Wilson coefficients? (y/n): ");
			scanf(" %c", &response_scalar);  // Space before %c to consume any leftover newline characters
			if (response_scalar == 'y' || response_scalar == 'Y') {
				// LFUV scalar Wilson Coefficients
				if (response_lfuv == 'y' || response_lfuv == 'Y') {

					char response_scalar_mode;

					// Electron
					printf("Do you want non-zero NP in the electron mode? (y/n): ");
					scanf(" %c", &response_scalar_mode);
					if (response_scalar_mode == 'y' || response_scalar_mode == 'Y') {
					double dCQ1e, dCQ2e;
					printf("Enter the values for dCQ1e, dCQ2e (separated by commas): ");
					scanf("%lf, %lf", &dCQ1e, &dCQ2e);

					param.KCQ1e = dCQ1e;
					param.KCQ2e = dCQ2e;
					}

					// Muon
					printf("Do you want non-zero NP in the muon mode? (y/n): ");
					scanf(" %c", &response_scalar_mode);
					if (response_scalar_mode == 'y' || response_scalar_mode == 'Y') {
					double dCQ1mu, dCQ2mu;
					printf("Enter the values for dCQ1mu, dCQ2mu (separated by commas): ");
					scanf("%lf, %lf", &dCQ1mu, &dCQ2mu);

					param.KCQ1mu = dCQ1mu;
					param.KCQ2mu = dCQ2mu;
					}

				} else {
					// Flavor-universal Wilson coefficients
					double dCQ1, dCQ2;
					printf("Enter the values for dCQ1, dCQ2 (separated by commas): ");
					scanf("%lf, %lf", &dCQ1, &dCQ2);

					param.KCQ1e = param.KCQ1mu = dCQ1;
					param.KCQ2e = param.KCQ2mu = dCQ2;
				}
			}


			// Ask about primed (pseudo-)scalar Wilson coefficients
			char response_scalar_rh;
			printf("Do you want to consider primed scalar Wilson coefficients? (y/n): ");
			scanf(" %c", &response_scalar_rh);  // Space before %c to consume any leftover newline characters
			if (response_scalar_rh == 'y' || response_scalar_rh == 'Y') {
				// LFUV Wilson primed scalar Coefficients
				if (response_lfuv == 'y' || response_lfuv == 'Y') {

					char response_scalar_mode;

					// Electron
					printf("Do you want non-zero NP in the electron mode? (y/n): ");
					scanf(" %c", &response_scalar_mode);
					if (response_scalar_mode == 'y' || response_scalar_mode == 'Y') {
					double dCQp1e, dCQp2e;
					printf("Enter the values for dCQp1e, dCQp2e (separated by commas): ");
					scanf("%lf, %lf", &dCQp1e, &dCQp2e);

					param.KCQ1pe = dCQp1e;
					param.KCQ2pe = dCQp2e;
					}

					// Muon
					printf("Do you want non-zero NP in the muon mode? (y/n): ");
					scanf(" %c", &response_scalar_mode);
					if (response_scalar_mode == 'y' || response_scalar_mode == 'Y') {
					double dCQp1mu, dCQp2mu;
					printf("Enter the values for dCQp1mu, dCQp2mu (separated by commas): ");
					scanf("%lf, %lf", &dCQp1mu, &dCQp2mu);

					param.KCQ1pmu = dCQp1mu;
					param.KCQ2pmu = dCQp2mu;
					}

				} else {
					// Flavor-universal Wilson coefficients
					double dCQp1, dCQp2;
					printf("Enter the values for dCQp1, dCQp2 (separated by commas): ");
					scanf("%lf, %lf", &dCQp1, &dCQp2);

					param.KCQ1pe = param.KCQ1pmu = dCQp1;
					param.KCQ2pe = param.KCQ2pmu = dCQp2;
				}
			}
		}
		/* *************************************************************** */
    }

    // Get prediction
    double *predictions;
    get_predictions(names, &nobs, &predictions, &param);

    double **covariance_th;
    get_th_covariance(&covariance_th, names, &nobs, &param);

    printf("\n");
    if (nobs==1) {
        pretty_print(names[ie], predictions[ie], sqrt(covariance_th[ie][ie]), PRECISION, 0);
    } else {
        for(ie=0;ie<nobs;ie++) {
            pretty_print(names[ie], predictions[ie], sqrt(covariance_th[ie][ie]), PRECISION, 30);
        }
    }
    printf("\n");

    // Save to CSV if outputname provided and more than one argument given
    if (argc >= 3 && outputname[0] != '\0') {
        save_predictions_to_csv(outputname, names, predictions, covariance_th, nobs);
    }

    // Free variables
    for (ie = 0; ie < nobs; ie++) {
        free(covariance_th[ie]);
    }
    free(covariance_th);
    free(predictions);

    return 0;
}









int check_nonphys_obs(obsname* obs)
/* Checks if the observable is a non-physical one. Returns 1 if True, 0 if False */
{
	if(
		(!strcmp(obs->type,"C1")||!strcmp(obs->type,"C2")||!strcmp(obs->type,"C3")||!strcmp(obs->type,"C4")||!strcmp(obs->type,"C5")||!strcmp(obs->type,"C6")||!strcmp(obs->type,"C7")||!strcmp(obs->type,"C8")||!strcmp(obs->type,"C9")||!strcmp(obs->type,"C10")||!strcmp(obs->type,"CQ1")||!strcmp(obs->type,"CQ2")||!strcmp(obs->type,"Cprime1")||!strcmp(obs->type,"Cprime2")||!strcmp(obs->type,"Cprime3")||!strcmp(obs->type,"Cprime4")||!strcmp(obs->type,"Cprime5")||!strcmp(obs->type,"Cprime6")||!strcmp(obs->type,"Cprime7")||!strcmp(obs->type,"Cprime8")||!strcmp(obs->type,"Cprime9")||!strcmp(obs->type,"Cprime10")||!strcmp(obs->type,"CprimeQ1")||!strcmp(obs->type,"CprimeQ2")||!strcmp(obs->type,"ReC1")||!strcmp(obs->type,"ReC2")||!strcmp(obs->type,"ReC3")||!strcmp(obs->type,"ReC4")||!strcmp(obs->type,"ReC5")||!strcmp(obs->type,"ReC6")||!strcmp(obs->type,"ReC7")||!strcmp(obs->type,"ReC8")||!strcmp(obs->type,"ReC9")||!strcmp(obs->type,"ReC10")||!strcmp(obs->type,"ReCQ1")||!strcmp(obs->type,"ReCQ2")||!strcmp(obs->type,"ReCprime1")||!strcmp(obs->type,"ReCprime2")||!strcmp(obs->type,"ReCprime3")||!strcmp(obs->type,"ReCprime4")||!strcmp(obs->type,"ReCprime5")||!strcmp(obs->type,"ReCprime6")||!strcmp(obs->type,"ReCprime7")||!strcmp(obs->type,"ReCprime8")||!strcmp(obs->type,"ReCprime9")||!strcmp(obs->type,"ReCprime10")||!strcmp(obs->type,"ReCprimeQ1")||!strcmp(obs->type,"ReCprimeQ2")||!strcmp(obs->type,"ImC1")||!strcmp(obs->type,"ImC2")||!strcmp(obs->type,"ImC3")||!strcmp(obs->type,"ImC4")||!strcmp(obs->type,"ImC5")||!strcmp(obs->type,"ImC6")||!strcmp(obs->type,"ImC7")||!strcmp(obs->type,"ImC8")||!strcmp(obs->type,"ImC9")||!strcmp(obs->type,"ImC10")||!strcmp(obs->type,"ImCQ1")||!strcmp(obs->type,"ImCQ2")||!strcmp(obs->type,"ImCprime1")||!strcmp(obs->type,"ImCprime2")||!strcmp(obs->type,"ImCprime3")||!strcmp(obs->type,"ImCprime4")||!strcmp(obs->type,"ImCprime5")||!strcmp(obs->type,"ImCprime6")||!strcmp(obs->type,"ImCprime7")||!strcmp(obs->type,"ImCprime8")||!strcmp(obs->type,"ImCprime9")||!strcmp(obs->type,"ImCprime10")||!strcmp(obs->type,"ImCprimeQ1")||!strcmp(obs->type,"ImCprimeQ2")||!strcmp(obs->type,"Vud")||!strcmp(obs->type,"Vus")||!strcmp(obs->type,"Vub")||!strcmp(obs->type,"Vcd")||!strcmp(obs->type,"Vcs")||!strcmp(obs->type,"Vcb")||!strcmp(obs->type,"Vtd")||!strcmp(obs->type,"Vts")||!strcmp(obs->type,"Vtb")||!strcmp(obs->type,"ReVud")||!strcmp(obs->type,"ReVus")||!strcmp(obs->type,"ReVub")||!strcmp(obs->type,"ReVcd")||!strcmp(obs->type,"ReVcs")||!strcmp(obs->type,"ReVcb")||!strcmp(obs->type,"ReVtd")||!strcmp(obs->type,"ReVts")||!strcmp(obs->type,"ReVtb")||!strcmp(obs->type,"ImVud")||!strcmp(obs->type,"ImVus")||!strcmp(obs->type,"ImVub")||!strcmp(obs->type,"ImVcd")||!strcmp(obs->type,"ImVcs")||!strcmp(obs->type,"ImVcb")||!strcmp(obs->type,"ImVtd")||!strcmp(obs->type,"ImVts")||!strcmp(obs->type,"ImVtb")
		)
	) return 1;// Non-physical

	return 0;// Is not among the non-physical list
}

int check_non_bTOs_decay(obsname* obs)
/* Checks if the decay is not b->s. Returns 1 if True, 0 if False */
{
	if(
		(!strcmp(obs->decay,"KLpi0nunu")||!strcmp(obs->decay,"Kpinunu")||!strcmp(obs->decay,"KLmumu")||!strcmp(obs->decay,"KLee")||!strcmp(obs->decay,"KSmumu")||!strcmp(obs->decay,"KLpi0ee")||!strcmp(obs->decay,"KLpi0mumu")
		)
	) return 1;// Is not among b->s decays

	return 0;// Is among b->s decays
}


// Function to map coefficient indices to readable names
const char* get_coefficient_name(const char* type, int index) {
    static const char* deltaC_names[] = {
        "", "", "", "", "", "", "", "deltaC7_mu", "", "deltaC9_mu", "deltaC10_mu",
        "", "", "", "", "", "", "deltaC7_e", "", "deltaC9_e", "deltaC10_e",
        "", "", "", "", "", "", "deltaC7_tau", "", "deltaC9_tau", "deltaC10_tau"
    };

    static const char* deltaCp_names[] = {
        "", "", "", "", "", "", "", "deltaC7prime_mu", "", "deltaC9prime_mu", "deltaC10prime_mu",
        "", "", "", "", "", "", "deltaC7prime_e", "", "deltaC9prime_e", "deltaC10prime_e",
        "", "", "", "", "", "", "deltaC7prime_tau", "", "deltaC9prime_tau", "deltaC10prime_tau"
    };

    static const char* deltaCQ_names[] = {
        "", "deltaCQ1_mu", "deltaCQ2_mu", "deltaCQ1_e", "deltaCQ2_e", "deltaCQ1_tau", "deltaCQ2_tau"
    };

    static const char* deltaCQp_names[] = {
        "", "deltaCQ1prime_mu", "deltaCQ2prime_mu", "deltaCQ1prime_e", "deltaCQ2prime_e",
        "deltaCQ1prime_tau", "deltaCQ2prime_tau"
    };

    if (strcmp(type, "deltaC") == 0) {
        if (index >= 0 && index < sizeof(deltaC_names) / sizeof(deltaC_names[0])) {
            return deltaC_names[index];
        }
    } else if (strcmp(type, "deltaCp") == 0) {
        if (index >= 0 && index < sizeof(deltaCp_names) / sizeof(deltaCp_names[0])) {
            return deltaCp_names[index];
        }
    } else if (strcmp(type, "deltaCQ") == 0) {
        if (index >= 0 && index < sizeof(deltaCQ_names) / sizeof(deltaCQ_names[0])) {
            return deltaCQ_names[index];
        }
    } else if (strcmp(type, "deltaCQp") == 0) {
        if (index >= 0 && index < sizeof(deltaCQp_names) / sizeof(deltaCQp_names[0])) {
            return deltaCQp_names[index];
        }
    }

    // If no match is found
    return "Unknown coefficient";
}



// Function to parse Wilson coefficients assignments
void parse_line(char* line, struct parameters* param, int* has_np_contributions) {
    static int in_comment = 0; // Track if inside a comment block

    while (1) {
        if (in_comment) {
            char* end = strstr(line, "*/");
            if (end) {
                in_comment = 0; // Exit comment block
                // Shift the remaining line after */
                memmove(line, end + 2, strlen(end + 2) + 1);
            } else {
                // Entire line is still within a comment block
                return;
            }
        }

        char* start = strstr(line, "/*");
        if (start) {
            char* end = strstr(start, "*/");
            if (end) {
                // Inline comment: truncate the comment and continue
                memmove(start, end + 2, strlen(end + 2) + 1);
            } else {
                in_comment = 1; // Enter comment block
                *start = '\0';  // Truncate at the start of the comment
            }
        }

        // Ignore empty or fully commented lines
        if (line[0] == '\0' || line[0] == '\n' || (line[0] == '/' && line[1] == '/')) {
            return;
        }

        // Parse the valid part of the line before the next comment
		int index = 0;
		double real = 0.0, imag = 0.0;
		double value;
		int SM_like_max_index = 30;
		int SM_like_prime_max_index = 30;
		int scalar_max_index = 6;
		int scalar_prime_max_index = 6;

        // Parse b->s coefficients first
		if (sscanf(line, "deltaC[%d] = %lf + I*(%lf)", &index, &real, &imag) == 3) {
			if (index > 0 && index <= SM_like_max_index) {
				param->deltaC[index] = real + I * imag;
				if (real != 0.0 || imag != 0.0) {
					*has_np_contributions = 1; // Mark that a non-zero value exists
					printf("Set %-18s = %f + I*(%f)\n", get_coefficient_name("deltaC", index), creal(param->deltaC[index]), cimag(param->deltaC[index]));
				}
			}
		} else if (sscanf(line, "deltaCp[%d] = %lf + I*(%lf)", &index, &real, &imag) == 3) {
			if (index > 0 && index <= SM_like_prime_max_index) {
				param->deltaCp[index] = real + I * imag;
				if (real != 0.0 || imag != 0.0) {
					*has_np_contributions = 1; // Mark that a non-zero value exists
					printf("Set %-18s = %f + I*(%f)\n", get_coefficient_name("deltaCp", index), creal(param->deltaCp[index]), cimag(param->deltaCp[index]));
				}
			}
		} else if (sscanf(line, "deltaCQ[%d] = %lf + I*(%lf)", &index, &real, &imag) == 3) {
			if (index > 0 && index <= scalar_max_index) {
				param->deltaCQ[index] = real + I * imag;
				if (real != 0.0 || imag != 0.0) {
					*has_np_contributions = 1; // Mark that a non-zero value exists
					printf("Set %-18s = %f + I*(%f)\n", get_coefficient_name("deltaCQ", index), creal(param->deltaCQ[index]), cimag(param->deltaCQ[index]));
				}
			}
		} else if (sscanf(line, "deltaCQp[%d] = %lf + I*(%lf)", &index, &real, &imag) == 3) {
			if (index >= 0 && index <= scalar_prime_max_index) {
				param->deltaCQp[index] = real + I * imag;
				if (real != 0.0 || imag != 0.0) {
					*has_np_contributions = 1; // Mark that a non-zero value exists
					printf("Set %-18s = %f + I*(%f)\n", get_coefficient_name("deltaCQp", index), creal(param->deltaCQp[index]), cimag(param->deltaCQp[index]));
				}
			}
		}
		// Parse s->d coefficients
		else if (sscanf(line, "KC9e = %lf + I*(%lf)", &real, &imag) == 2) {
            param->KC9e = real + I * imag;
            if (real != 0.0 || imag != 0.0) {
                *has_np_contributions = 1;
                printf("Set %-18s = %f + I*(%f)\n", "KC9e", creal(param->KC9e), cimag(param->KC9e));
            }
        } else if (sscanf(line, "KC9mu = %lf + I*(%lf)", &real, &imag) == 2) {
            param->KC9mu = real + I * imag;
            if (real != 0.0 || imag != 0.0) {
                *has_np_contributions = 1;
                printf("Set %-18s = %f + I*(%f)\n", "KC9mu", creal(param->KC9mu), cimag(param->KC9mu));
            }
        } else if (sscanf(line, "KC10e = %lf + I*(%lf)", &real, &imag) == 2) {
            param->KC10e = real + I * imag;
            if (real != 0.0 || imag != 0.0) {
                *has_np_contributions = 1;
                printf("Set %-18s = %f + I*(%f)\n", "KC10e", creal(param->KC10e), cimag(param->KC10e));
            }
        } else if (sscanf(line, "KC10mu = %lf + I*(%lf)", &real, &imag) == 2) {
            param->KC10mu = real + I * imag;
            if (real != 0.0 || imag != 0.0) {
                *has_np_contributions = 1;
                printf("Set %-18s = %f + I*(%f)\n", "KC10mu", creal(param->KC10mu), cimag(param->KC10mu));
            }
        } else if (sscanf(line, "KC9pe = %lf + I*(%lf)", &real, &imag) == 2) {
            param->KC9pe = real + I * imag;
            if (real != 0.0 || imag != 0.0) {
                *has_np_contributions = 1;
                printf("Set %-18s = %f + I*(%f)\n", "KC9pe", creal(param->KC9pe), cimag(param->KC9pe));
            }
        } else if (sscanf(line, "KC9pmu = %lf + I*(%lf)", &real, &imag) == 2) {
            param->KC9pmu = real + I * imag;
            if (real != 0.0 || imag != 0.0) {
                *has_np_contributions = 1;
                printf("Set %-18s = %f + I*(%f)\n", "KC9pmu", creal(param->KC9pmu), cimag(param->KC9pmu));
            }
        } else if (sscanf(line, "KC10pe = %lf + I*(%lf)", &real, &imag) == 2) {
            param->KC10pe = real + I * imag;
            if (real != 0.0 || imag != 0.0) {
                *has_np_contributions = 1;
                printf("Set %-18s = %f + I*(%f)\n", "KC10pe", creal(param->KC10pe), cimag(param->KC10pe));
            }
        } else if (sscanf(line, "KC10pmu = %lf + I*(%lf)", &real, &imag) == 2) {
            param->KC10pmu = real + I * imag;
            if (real != 0.0 || imag != 0.0) {
                *has_np_contributions = 1;
                printf("Set %-18s = %f + I*(%f)\n", "KC10pmu", creal(param->KC10pmu), cimag(param->KC10pmu));
            }
        } else if (sscanf(line, "KCQ1e = %lf + I*(%lf)", &real, &imag) == 2) {
            param->KCQ1e = real + I * imag;
            if (real != 0.0 || imag != 0.0) {
                *has_np_contributions = 1;
                printf("Set %-18s = %f + I*(%f)\n", "KCQ1e", creal(param->KCQ1e), cimag(param->KCQ1e));
            }
        } else if (sscanf(line, "KCQ1mu = %lf + I*(%lf)", &real, &imag) == 2) {
            param->KCQ1mu = real + I * imag;
            if (real != 0.0 || imag != 0.0) {
                *has_np_contributions = 1;
                printf("Set %-18s = %f + I*(%f)\n", "KCQ1mu", creal(param->KCQ1mu), cimag(param->KCQ1mu));
            }
        } else if (sscanf(line, "KCQ2e = %lf + I*(%lf)", &real, &imag) == 2) {
            param->KCQ2e = real + I * imag;
            if (real != 0.0 || imag != 0.0) {
                *has_np_contributions = 1;
                printf("Set %-18s = %f + I*(%f)\n", "KCQ2e", creal(param->KCQ2e), cimag(param->KCQ2e));
            }
        } else if (sscanf(line, "KCQ2mu = %lf + I*(%lf)", &real, &imag) == 2) {
            param->KCQ2mu = real + I * imag;
            if (real != 0.0 || imag != 0.0) {
                *has_np_contributions = 1;
                printf("Set %-18s = %f + I*(%f)\n", "KCQ2mu", creal(param->KCQ2mu), cimag(param->KCQ2mu));
            }
        } else if (sscanf(line, "KCQ1pe = %lf + I*(%lf)", &real, &imag) == 2) {
            param->KCQ1pe = real + I * imag;
            if (real != 0.0 || imag != 0.0) {
                *has_np_contributions = 1;
                printf("Set %-18s = %f + I*(%f)\n", "KCQ1pe", creal(param->KCQ1pe), cimag(param->KCQ1pe));
            }
        } else if (sscanf(line, "KCQ1pmu = %lf + I*(%lf)", &real, &imag) == 2) {
            param->KCQ1pmu = real + I * imag;
            if (real != 0.0 || imag != 0.0) {
                *has_np_contributions = 1;
                printf("Set %-18s = %f + I*(%f)\n", "KCQ1pmu", creal(param->KCQ1pmu), cimag(param->KCQ1pmu));
            }
        } else if (sscanf(line, "KCQ2pe = %lf + I*(%lf)", &real, &imag) == 2) {
            param->KCQ2pe = real + I * imag;
            if (real != 0.0 || imag != 0.0) {
                *has_np_contributions = 1;
                printf("Set %-18s = %f + I*(%f)\n", "KCQ2pe", creal(param->KCQ2pe), cimag(param->KCQ2pe));
            }
        } else if (sscanf(line, "KCQ2pmu = %lf + I*(%lf)", &real, &imag) == 2) {
            param->KCQ2pmu = real + I * imag;
            if (real != 0.0 || imag != 0.0) {
                *has_np_contributions = 1;
                printf("Set %-18s = %f + I*(%f)\n", "KCQ2pmu", creal(param->KCQ2pmu), cimag(param->KCQ2pmu));
            }
        } else if (sscanf(line, "KCLe = %lf + I*(%lf)", &real, &imag) == 2) {
            param->KCLe = real + I * imag;
            if (real != 0.0 || imag != 0.0) {
                *has_np_contributions = 1;
                printf("Set %-18s = %f + I*(%f)\n", "KCLe", creal(param->KCLe), cimag(param->KCLe));
            }
        } else if (sscanf(line, "KCLmu = %lf + I*(%lf)", &real, &imag) == 2) {
            param->KCLmu = real + I * imag;
            if (real != 0.0 || imag != 0.0) {
                *has_np_contributions = 1;
                printf("Set %-18s = %f + I*(%f)\n", "KCLmu", creal(param->KCLmu), cimag(param->KCLmu));
            }
        } else if (sscanf(line, "KCLtau = %lf + I*(%lf)", &real, &imag) == 2) {
            param->KCLtau = real + I * imag;
            if (real != 0.0 || imag != 0.0) {
                *has_np_contributions = 1;
                printf("Set %-18s = %f + I*(%f)\n", "KCLtau", creal(param->KCLtau), cimag(param->KCLtau));
            }
        } else if (sscanf(line, "KCRe = %lf + I*(%lf)", &real, &imag) == 2) {
            param->KCRe = real + I * imag;
            if (real != 0.0 || imag != 0.0) {
                *has_np_contributions = 1;
                printf("Set %-18s = %f + I*(%f)\n", "KCRe", creal(param->KCRe), cimag(param->KCRe));
            }
        } else if (sscanf(line, "KCRmu = %lf + I*(%lf)", &real, &imag) == 2) {
            param->KCRmu = real + I * imag;
            if (real != 0.0 || imag != 0.0) {
                *has_np_contributions = 1;
                printf("Set %-18s = %f + I*(%f)\n", "KCRmu", creal(param->KCRmu), cimag(param->KCRmu));
            }
        } else if (sscanf(line, "KCRtau = %lf + I*(%lf)", &real, &imag) == 2) {
            param->KCRtau = real + I * imag;
            if (real != 0.0 || imag != 0.0) {
                *has_np_contributions = 1;
                printf("Set %-18s = %f + I*(%f)\n", "KCRtau", creal(param->KCRtau), cimag(param->KCRtau));
            }
        } else {
            // Only print "Invalid line" for non-empty, non-comment lines that don't match any pattern
            if (strlen(line) > 1 && line[0] != '/' && line[0] != '\n') {
                printf("Invalid line: %s", line);
            }
        }

        // Break if no more comments are expected in this line
        if (!in_comment) {
            break;
        }
    }
}


// Helper function to get the exponent of a number in scientific notation
int get_exponent(double number) {
    if (number == 0.0) return 0;
    return (int)floor(log10(fabs(number)));
}

// Helper function to format uncertainty with matched exponent and adjustable precision
void format_with_matched_exponent(double main_value, double uncertainty, char *formatted_str, int precision) {
    int main_exponent = get_exponent(main_value);
    double scaled_uncertainty = uncertainty / pow(10, main_exponent);

    // Create the format string with the specified precision
    char format[10];
    snprintf(format, sizeof(format), "%%.%df", precision);  // e.g., "%.3f" if precision is 3
    snprintf(formatted_str, 20, format, scaled_uncertainty);
}

// Main pretty-printing function with adjustable precision
void pretty_print(char *str, double a, double b, int precision, int padding_size) {
    char format[50];            // Buffer for dynamic format string
    char uncertainty_str[20];   // Buffer for formatted uncertainty

    // Dynamically construct the format string based on the padding size
    if (padding_size > 0) {
        snprintf(format, sizeof(format), "%%-%ds = ", padding_size);  // Left-aligned with padding
    } else {
        strcpy(format, "%s = ");  // No padding
    }

    // Print the formatted string
    printf(format, str);

    // Create format strings for the main value and uncertainty
    char main_format[10];
    snprintf(main_format, sizeof(main_format), "%% .%df", precision);  // e.g., "%.3f" if precision is 3

    char scientific_format[10];
    snprintf(scientific_format, sizeof(scientific_format), "%% .%de", precision);  // e.g., "%.3e" if precision is 3

    // Print the main value `a`
    if (fabs(a) >= 0.0001) {
        printf(main_format, a);  // Regular floating-point notation
    } else {
        printf(scientific_format, a);  // Scientific notation for small values
    }

    printf(" ± ");

    // Print the uncertainty `b`
    if (fabs(a) >= 0.0001) {
		printf("\b");//backspace for removing extra space not needed for uncertainties (always positive)
        printf(main_format, b);  // Regular notation if main value is in regular notation
    } else {
        // Match exponent for `b` if main value `a` is in scientific notation
        format_with_matched_exponent(a, b, uncertainty_str, precision);
        printf("%se%d", uncertainty_str, get_exponent(a));
    }

    printf("\n");
}

void save_predictions_to_csv(const char* filename, char names[][MAX_OBS_NAME_LEN],
                           double* predictions, double** covariance_th, int nobs) {
    // Ensure output directory exists
    struct stat st = {0};
    if (stat("dataOutput", &st) == -1) {
        mkdir("dataOutput", 0700); // Create directory with read/write/search permissions for owner
    }

    char csv_filename[512];
    snprintf(csv_filename, sizeof(csv_filename), "dataOutput/%s.csv", filename);

    // Check if file exists
    FILE* test = fopen(csv_filename, "r");
    if (test != NULL) {
        fclose(test);
        char response;
        printf("\nFile %s already exists. Overwrite? (y/n): ", csv_filename);
        scanf(" %c", &response);
        if (response != 'y' && response != 'Y') {
            printf("Predictions not saved.\n");
            return;
        }
    }

    FILE* file = fopen(csv_filename, "w");
    if (!file) {
        perror("Error opening file for writing");
        return;
    }


    // Write data
    for (int i = 0; i < nobs; i++) {
        fprintf(file, "%s,%.6e,%.6e\n", names[i], predictions[i], sqrt(covariance_th[i][i]));
    }

    fclose(file);
    printf("\nPredictions saved to %s\n", csv_filename);
}
