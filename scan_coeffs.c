#include "src/include.h"
#include <sys/time.h>

/*-----------------------------------------------------------------*/
/* Calculation of chi2 scan with choice of  NP Wilson coefficients */
/*-----------------------------------------------------------------*/

/* Change the number of point and scan range with NUM_POINTS, RANGE_MIN and RANGE_MAX below  */
#define NUM_POINTS 1000
#define RANGE_MIN -2.5
#define RANGE_MAX 2.5

/* Change numerical value for modifying precision in the float type output (WC and chi2) and exponent type output (obs prediction) */
#define FMT_FLOAT "%.6f"
#define FMT_EXPONENT "%.6e"

#define DELTA_CHI2_MAX 106.0  // 106.0 ~ 10 sigma for 2dim fit

// Set to 0 to use SM covariance, to 1 to recalculate for each point  (makes the running time significantly longer)
#define RECALCULATE_COV_MATRIX 0


// Global variables for time tracking
static struct timeval start_time;
static int is_first_progress = 1;

void print_progress_bar(int iteration, int total, int show_time);
int is_file_non_empty(const char* filename);
int read_coeff(double parameter_values[], struct parameters *param, char coeff_names[][50], int *ncoeff);
// void write_header(FILE *output_file, char coeff_names[][50], int ncoeffs);
// void write_output_header(FILE *output_file, char coeff_names[][50], int ncoeffs,
//                         char obs_names[][50], int nobs, int record_predictions);
void print_current_time_in_timezone();
void apply_coeff_value(struct parameters *param, const char* name, double value);
void get_datetime_string(char *buffer, size_t size);
void display_scan_instructions();

void display_scan_instructions() {
    printf("\n");
    printf("***************************************************************************************\n");
    printf("*************** Wilson Coefficient Scanner for Flavour Observables ********************\n");
    printf("***************************************************************************************\n");
    printf("* \033[1;36mProgram Overview:\033[0m                                                                   *\n");
    printf("* - Performs random scan over Wilson coefficients and calculates chi2 values          *\n");
    printf("* - Outputs results in CSV format for easy analysis                                   *\n");
    printf("* - Includes Standard Model prediction as first data point                            *\n");
    printf("*                                                                                     *\n");
    printf("* \033[1;36mConfiguration Options (edit in code):\033[0m                                               *\n");
    printf("* - #define NUM_POINTS %-5d    // Number of random points to generate                *\n", NUM_POINTS);
    printf("* - #define RANGE_MIN %-5.1f     // Minimum value for coefficients                     *\n", RANGE_MIN);
    printf("* - #define RANGE_MAX %-5.1f     // Maximum value for coefficients                     *\n", RANGE_MAX);
    printf("*                                                                                     *\n");
    printf("* \033[1;36mAdditional Options:\033[0m                                                                 *\n");
    printf("* + \033[1;32m-p\033[0m              : Record predictions for all observables (creates larger file)    *\n");
    printf("* + \033[1;32m-pe\033[0m             : Record predictions AND errors for all observables               *\n");
    printf("*                                                                                     *\n");
    printf("* \033[1;36mRequired Arguments:\033[0m                                                                 *\n");
    printf("* 1) \033[1;32mobs_file\033[0m       : File with list of observables to calculate                      *\n");
    printf("* 2) \033[1;32mcoeff_file\033[0m     : File with list of Wilson coefficients to scan                   *\n");
    printf("*                                                                                     *\n");
    printf("* \033[1;36mOptional Arguments:\033[0m                                                                 *\n");
    printf("* \033[1;36mAdditional Options:\033[0m                                                                 *\n");
    printf("* 3) \033[1;32moutput_name\033[0m    : Base name for output file (default: auto-generated)             *\n");
    printf("* 4) \033[1;32mexp_file\033[0m       : Experimental values file                                        *\n");
    printf("* 5) \033[1;32mexp_corr_file\033[0m  : Experimental correlations file                                  *\n");
    printf("* 6) \033[1;32mnuis_file\033[0m      : Nuisance parameters file                                        *\n");
    printf("* 7) \033[1;32mnuis_corr_file\033[0m : Nuisance correlations file                                      *\n");
    printf("*                                                                                     *\n");
    printf("* \033[1;36mExample Usage:\033[0m                                                                      *\n");
    printf("* \033[1;32m./scan_coeffs.x -p chi2_input/obs_BKstrmm_Pi_LHCb.in chi2_input/coeff_list.txt     \033[0m *\n");
    printf("* \033[1;32m./scan_coeffs.x chi2_input/obs_BKstrmm_Pi_LHCb.in chi2_input/coeff_list.txt my_scan\033[0m *\n");
    printf("***************************************************************************************\n\n");
}


int main(int argc, char** argv) {

    // Including or not the individual predictions and errors
    int record_predictions = 0;
    int record_errors = 0;

    // Parse for -p and -pe flags before your existing argument processing
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-p") == 0) {
            record_predictions = 1;
            // Remove -p from arguments
            for (int j = i; j < argc - 1; j++) {
                argv[j] = argv[j + 1];
            }
            argc--;
            i--;
            printf("Will record predictions for all observables\n");
        }
        else if (strcmp(argv[i], "-pe") == 0) {
            record_predictions = 1;
            record_errors = 1;
            // Remove -pe from arguments
            for (int j = i; j < argc - 1; j++) {
                argv[j] = argv[j + 1];
            }
            argc--;
            i--;
            printf("Will record predictions AND errors for all observables\n");
        }
    }

    char obsfile[500], coefffile[500], expfile[500], corrfile[500], nuisfile[500], nuiscorrfile[500];
    char output_base[100];
	printf("SuperIso " VERSION " - %d\n\n", YEAR);

    // Show instructions if no arguments provided
    if(argc < 2) {
        display_scan_instructions();
        fprintf(stderr, "\033[1;31mError: At least 2 arguments are required.\033[0m\n");
        fprintf(stderr, "Example:\n");
        fprintf(stderr, "  \033[1;32m./scan_coeffs.x chi2_input/obs_BKstrmm_Pi_LHCb.in chi2_input/coeff_list.txt\033[0m\n");
        exit(EXIT_FAILURE);
    }

    // Ensure the `scans` directory exists (unchanged)
    if (access("scans", F_OK) != 0) {
		// Directory doesn't exist, attempt to create it
        if (mkdir("scans", 0755) != 0) {
            perror("Error creating 'scans' directory");
            return 1;
        }
        printf("Directory 'scans' created successfully.\n");
    }


  	if(argc<3)
  	{
		// Less than 2 arguments: Print usage example and exit
		fprintf(stderr, "Error: At least 2 arguments are required.\n");
		fprintf(stderr, "Example:\n");
		fprintf(stderr, "  ./scan_coeffs_general chi2_input/obs_BKstrmm_Pi_LHCb.in chi2_input/coeff_list.txt \n");
		fprintf(stderr, "  where in 'chi2_input/coeff_list.txt' relevant coefficients should be uncommented.\n");
		fprintf(stderr, "Optional arguments: output_file, exp_file, exp_cor_file, nuis_file, nuis_corr_file\n");
		exit(EXIT_FAILURE);
  	} else {
		// First argument: observable input file
		strncpy(obsfile, argv[1], sizeof(obsfile) - 1);
		obsfile[sizeof(obsfile) - 1] = '\0';
		printf("Reading observable list from %s\n", obsfile);

		// Second argument: coefficient input file
		strncpy(coefffile, argv[2], sizeof(coefffile) - 1);
		coefffile[sizeof(coefffile) - 1] = '\0';
		printf("Reading coefficients to scan from %s\n", coefffile);

	}

	// Third argument: output filename (optional)
    // Handle output filename
    if (argc >= 4) {
        // User-specified filename
        strncpy(output_base, argv[3], sizeof(output_base) - 1);
        output_base[sizeof(output_base) - 1] = '\0';
    } else {
        // Auto-generate chi2_<datetime>.csv filename
        char datetime[16]; // YYYYMMDD_HHMMSS + null terminator
        get_datetime_string(datetime, sizeof(datetime));
        snprintf(output_base, sizeof(output_base), "%s", datetime);
    }

    // Ensure .csv extension
    char chi2_output_filename[256];
    // Strip existing extension if present
    char *dot = strrchr(output_base, '.');
    if (dot) *dot = '\0';  // Truncate at last dot

    if (record_predictions) {
        if (record_errors) {
            snprintf(chi2_output_filename, sizeof(chi2_output_filename),
                    "scans/%s_with_predictions_and_errors.csv", output_base);
        } else {
            snprintf(chi2_output_filename, sizeof(chi2_output_filename),
                    "scans/%s_with_predictions.csv", output_base);
        }
    } else {
        snprintf(chi2_output_filename, sizeof(chi2_output_filename),
                "scans/%s.csv", output_base);
    }


    if(argc>=5) {
        sscanf(argv[4],"%s",expfile);
		printf("Reading user-provided experimental values from %s\n",expfile);
    } else {
		sprintf(expfile,"chi2_input/exp_values.in");
    }

    if(argc>=6) {
        sscanf(argv[5],"%s",corrfile);
		printf("Reading user-provided experimental correlations from %s\n",corrfile);
    } else {
		sprintf(corrfile,"chi2_input/exp_corr.in");
    }

    if(argc>=7) {
        sscanf(argv[6],"%s",nuisfile);
		printf("Reading user-provided nuisance parameters from %s\n",nuisfile);
    } else {
		sprintf(nuisfile,"chi2_input/nuisance.in");
    }

    if(argc>=8) {
        sscanf(argv[7],"%s",nuiscorrfile);
		printf("Reading user-provided nuisance parameter correlations from %s\n",nuiscorrfile);
    } else {
		sprintf(nuiscorrfile,"chi2_input/nuisance_corr.in");
    }


    // Load parameters and setup
    struct parameters param;
    char names[NBOBSMAX][50], coeff_names[50][50];
    int nobs = read_obs_list(obsfile, names);
    int ncoeff = read_obs_list(coefffile, coeff_names);

    // Check if any coefficients were read
    if (ncoeff == 0) {
        fprintf(stderr, "\033[1;31mError: No Wilson coefficients found in %s\033[0m\n", coefffile);
        fprintf(stderr, "Please uncomment at least one Wilson coefficient in the file to scan.\n");
        fprintf(stderr, "Example coefficient lines to uncomment for scan over real C9 and C10:\n");
        fprintf(stderr, "# dC9_real\n");
        fprintf(stderr, "# dC19_real\n");
        return 1;
    }


    Init_param(&param);
    slha_adjust(&param);
    param.SM = 1;

	snprintf(param.exp_values, sizeof(param.exp_values), "%s", expfile);
	snprintf(param.exp_corr, sizeof(param.exp_corr), "%s", corrfile);
	snprintf(param.nuisance_values, sizeof(param.nuisance_values), "%s", nuisfile);
	snprintf(param.nuisance_corr, sizeof(param.nuisance_corr), "%s", nuiscorrfile);

    double* predictions = NULL;
    get_predictions(names, &nobs, &predictions, &param);

	// // Allocate and calculate covariance matrices
	// printf("Covariance matrix being calculated for the first time (can take a while)...\n");
	// double *central_exp, **covariance_exp, **covariance_th, **covariance_tot, **inv_cov_tot;
	// get_exp_covariance(&covariance_exp, &central_exp, names, &nobs, &param);
	// get_th_covariance(&covariance_th, names, &nobs, &param);
	// get_covtot(&covariance_th, &covariance_exp, &covariance_tot, nobs);
    // Allocate and calculate covariance matrices
    printf("Covariance matrix being calculated for first time (can take a while)...\n");
    double *central_exp, **covariance_exp, **covariance_th, **covariance_tot, **inv_cov_tot;
    get_exp_covariance(&covariance_exp, &central_exp, names, &nobs, &param);
    get_th_covariance(&covariance_th, names, &nobs, &param);
    get_covtot(&covariance_th, &covariance_exp, &covariance_tot, nobs);

    // Invert covariance matrix
    if (get_invcovtot(&covariance_tot, &inv_cov_tot, nobs) < 0) {
        fprintf(stderr, "Error: Failed to invert covariance matrix.\n");
        return -3;
    }

    // Store SM covariance for later comparison if needed
    double **covariance_th_SM = NULL;
    if (RECALCULATE_COV_MATRIX) {
        // Allocate space to store SM covariance
        covariance_th_SM = (double **)malloc(nobs * sizeof(double *));
        for (int i = 0; i < nobs; i++) {
            covariance_th_SM[i] = (double *)malloc(nobs * sizeof(double));
            for (int j = 0; j < nobs; j++) {
                covariance_th_SM[i][j] = covariance_th[i][j];
            }
        }
        printf("Will recalculate theoretical covariance matrix for each point\n");
    } else {
        printf("Using fixed SM theoretical covariance matrix for all points\n");
    }

	// Invert covariance matrix
	if (get_invcovtot(&covariance_tot, &inv_cov_tot, nobs) < 0) {
		fprintf(stderr, "Error: Failed to invert covariance matrix.\n");
		return -3;
	}

	// Prepare chi2 output file
	// char chi2_output_filename[150];
	snprintf(chi2_output_filename, sizeof(chi2_output_filename), "scans/chi2_%s.csv", output_base);

	// Check if output file exists
	int file_exists = is_file_non_empty(chi2_output_filename);
	char mode[3]; // To hold the mode string for fopen
	int write_header = 1; // Flag to determine whether to write the header


	if (file_exists) {
        char choice;
        printf("Warning: File '%s' exists and is non-empty.\n", chi2_output_filename);
        printf("Do you want to (o)verwrite, (a)ppend, or (c)ancel? (o/a/c): ");
        fflush(stdout);
        scanf(" %c", &choice);

        if (choice == 'o' || choice == 'O') {
            strcpy(mode, "w+"); // Overwrite mode
            write_header = 1;
        } else if (choice == 'a' || choice == 'A') {
            strcpy(mode, "a"); // Append mode
            write_header = 0;
        } else {
            printf("Exiting program without modifications.\n");
            return 0;
        }
    } else {
        // File does not exist or is empty
        strcpy(mode, "w+");
        write_header = 1;
    }


	// Open output file
    FILE* output = fopen(chi2_output_filename, mode);
    if (!output) {
        fprintf(stderr, "Error: Could not open file %s for writing.\n", chi2_output_filename);
        return -1;
    }

    // Write header if needed
    if (write_header) {
        // Write coefficient headers
        for (int i = 0; i < ncoeff; i++) {
            fprintf(output, "%s,", coeff_names[i]);
        }

        // Write chi2 header
        fprintf(output, "chi2");

        // Write observable headers if predictions are recorded
        if (record_predictions) {
            for (int i = 0; i < nobs; i++) {
                if (record_errors) {
                    fprintf(output, ",%s,%s_err", names[i], names[i]);
                } else {
                    fprintf(output, ",%s", names[i]);
                }
            }
        }
        fprintf(output, "\n");
    }

    // Calculate and write SM chi2
    double parameter_values[50] = {0};  // All zeros for SM
    read_coeff(parameter_values, &param, coeff_names, &ncoeff);

    get_predictions(names, &nobs, &predictions, &param);
    double chi2_SM = get_chi2(inv_cov_tot, predictions, central_exp, nobs);

    // Write SM line (use parameter_values directly)
    for (int j = 0; j < ncoeff; j++) {
        fprintf(output, FMT_FLOAT ",", parameter_values[j]);
    }
    fprintf(output, FMT_FLOAT, chi2_SM);  // Note: no \n yet

    // Write predictions and errors if enabled
    if (record_predictions) {
        for (int i = 0; i < nobs; i++) {
            fprintf(output, "," FMT_EXPONENT, predictions[i]);
            if (record_errors) {
                double uncertainty = sqrt(covariance_th[i][i]);
                fprintf(output, "," FMT_EXPONENT, uncertainty);
            }
        }
    }
    fprintf(output, "\n");
    free(predictions);


	// Variables for progress bar
    int total_iterations = NUM_POINTS;
    int current_iteration = 1;
    int nSaved = 0;  // Counter for saved points



    // Main scan loop (unchanged except debug prints)
    // printf("Scan for %d chi2 points for %d observables,\n", NUM_POINTS, nobs);
    // printf("varying %d Wilson coefficients:\n", ncoeff);
    // for (int i = 0; i < ncoeff; i++) {
    //     printf("%s ", coeff_names[i]);
    // }
    // printf("\n");
    // printf("each in the [%.1f, %.1f] range\n", RANGE_MIN, RANGE_MAX);
    // printf("\nProducing %d chi2 points with scan over %d observables\n", NUM_POINTS, nobs);
    // printf("Scan range: [%.1f, %.1f]\n", RANGE_MIN, RANGE_MAX);
    // printf("Producing %d chi2 points with scan over %d observables\n", NUM_POINTS, nobs);
    // printf("Scan range: [%.1f, %.1f]\n", RANGE_MIN, RANGE_MAX);
    // if (RECALCULATE_COV_MATRIX) {
    //     printf("Theoretical covariance matrix will be recalculated for each point\n");
    // } else {
    //     printf("Using fixed SM theoretical covariance matrix\n");
    // }

	// Call the function to print the current time in CET
	print_current_time_in_timezone();


	// double parameter_values[50]; // Array to store parameter values during the scan
	// Start measuring time
	clock_t start_time = clock();
	srand(time(NULL)); // Non-deterministic sequence
	// const unsigned int FIXED_SEED = 42;  // Any integer works (42 is classic)
	// srand(FIXED_SEED);                   // Deterministic sequence
	/* ==== Main scanning loop (OLD VERSION) ==== */
	// for (int i = 0; i < NUM_POINTS; i++) {
	// 	// Generate random parameter values
	// 	for (int j = 0; j < ncoeff; j++) {
	// 		parameter_values[j] = RANGE_MIN + (RANGE_MAX - RANGE_MIN) * rand() / (RAND_MAX + 1.0);
	// 	}
 //
	// 	// Initialize parameters based on scanned values
	// 	read_coeff(parameter_values, &param, coeff_names, &ncoeff);
	// 	/* Testing produced Wilson coefficients */
	// 	// for (int k = 0; k < ncoeff; k++) {
	// 	// 	printf("parameter_values[%d]=%.3f ", k, parameter_values[k]);
	// 	// }
	// 	// printf("\n");
	// 	// printf("C[7]=%.3f+I*(%.3f),\tC[8]=%.3f+I*(%.3f)\n", creal(param.deltaC[7]), cimag(param.deltaC[7]), creal(param.deltaC[8]), cimag(param.deltaC[8]));
	// 	// printf("Cp[7]=%.3f+I*(%.3f),\tCp[8]=%.3f+I*(%.3f)\n", creal(param.deltaCp[7]), cimag(param.deltaCp[7]), creal(param.deltaCp[8]), cimag(param.deltaCp[8]));
	// 	// printf("C[9]=%.3f+I*(%.3f),\tC[19]=%.3f+I*(%.3f),\tC[29]=%.3f+I*(%.3f)\n", creal(param.deltaC[9]), cimag(param.deltaC[9]), creal(param.deltaC[19]), cimag(param.deltaC[19]), creal(param.deltaC[29]),cimag(param.deltaC[29]));
	// 	// printf("C[10]=%.3f+I*(%.3f),\tC[20]=%.3f+I*(%.3f),\tC[30]=%.3f+I*(%.3f)\n", creal(param.deltaC[10]), cimag(param.deltaC[10]), creal(param.deltaC[20]), cimag(param.deltaC[20]), creal(param.deltaC[30]),cimag(param.deltaC[30]));
	// 	// printf("Cp[9]=%.3f+I*(%.3f),\tCp[19]=%.3f+I*(%.3f),\tCp[29]=%.3f+I*(%.3f)\n", creal(param.deltaCp[9]), cimag(param.deltaCp[9]), creal(param.deltaCp[19]), cimag(param.deltaCp[19]), creal(param.deltaCp[29]),cimag(param.deltaCp[29]));
	// 	// printf("Cp[10]=%.3f+I*(%.3f),\tCp[20]=%.3f+I*(%.3f),\tCp[30]=%.3f+I*(%.3f)\n", creal(param.deltaCp[10]), cimag(param.deltaCp[10]), creal(param.deltaCp[20]), cimag(param.deltaCp[20]), creal(param.deltaCp[30]),cimag(param.deltaCp[30]));
	// 	// printf("CQ[1]=%.3f+I*(%.3f),\tCQ[3]=%.3f+I*(%.3f),\tCQ[5]=%.3f+I*(%.3f)\n", creal(param.deltaCQ[1]), cimag(param.deltaCQ[1]), creal(param.deltaCQ[3]),cimag(param.deltaCQ[3]), creal(param.deltaCQ[5]), cimag(param.deltaCQ[5]));
	// 	// printf("CQ[2]=%.3f+I*(%.3f),\tCQ[4]=%.3f+I*(%.3f),\tCQ[6]=%.3f+I*(%.3f)\n", creal(param.deltaCQ[2]), cimag(param.deltaCQ[2]),creal(param.deltaCQ[4]), cimag(param.deltaCQ[4]), creal(param.deltaCQ[6]),cimag(param.deltaCQ[6]));
	// 	// printf("CQp[1]=%.3f+I*(%.3f),\tCQp[3]=%.3f+I*(%.3f),\tCQp[5]=%.3f+I*(%.3f)\n", creal(param.deltaCQp[1]), cimag(param.deltaCQp[1]), creal(param.deltaCQp[3]),cimag(param.deltaCQp[3]), creal(param.deltaCQp[5]), cimag(param.deltaCQp[5]));
	// 	// printf("CQp[2]=%.3f+I*(%.3f),\tCQp[4]=%.3f+I*(%.3f),\tCQp[6]=%.3f+I*(%.3f)\n", creal(param.deltaCQp[2]), cimag(param.deltaCQp[2]),creal(param.deltaCQp[4]), cimag(param.deltaCQp[4]), creal(param.deltaCQp[6]),cimag(param.deltaCQp[6]));
 //
 //
 //
 //
 //
	// 	// Get predictions and chi2
	// 	get_predictions(names, &nobs, &predictions, &param);
	// 	double chi2 = get_chi2(inv_cov_tot, predictions, central_exp, nobs);
 //
	// 	// Write results to output file
 //        for (int j = 0; j < ncoeff; j++) {
 //            fprintf(output, FMT_FLOAT ",", parameter_values[j]);
 //        }
 //        fprintf(output, FMT_FLOAT, chi2);  // No \n yet
 //
 //        // Add predictions if enabled
 //        if (record_predictions) {
 //            for (int i = 0; i < nobs; i++) {
 //                fprintf(output, "," FMT_EXPONENT, predictions[i]);
 //            }
 //        }
 //        fprintf(output, "\n");
	// 	fflush(output);
 //
	// 	free(predictions);
	// 	// print_progress_bar(i + 1, NUM_POINTS);
 //        print_progress_bar(i + 1, NUM_POINTS, 1); // The '1' enables time display
	// }


    // Main scanning loop
    if (ncoeff == 1) {
        // Automatic ordered scan for 1D case
        printf("Single coefficient detected - performing ORDERED scan for %s\n", coeff_names[0]);
        printf("Range: [%.1f, %.1f] with %d points\n", RANGE_MIN, RANGE_MAX, NUM_POINTS);

        for (int i = 0; i < NUM_POINTS; i++) {
            // Linear spacing from RANGE_MIN to RANGE_MAX
            parameter_values[0] = RANGE_MIN + (RANGE_MAX - RANGE_MIN) * i / (NUM_POINTS - 1.0);

            // Initialize parameters based on scanned values
            read_coeff(parameter_values, &param, coeff_names, &ncoeff);

            // // Get predictions and chi2
            // get_predictions(names, &nobs, &predictions, &param);
            // double chi2 = get_chi2(inv_cov_tot, predictions, central_exp, nobs);

            // Get predictions
            get_predictions(names, &nobs, &predictions, &param);

            // Recalculate covariance matrix if enabled
            if (RECALCULATE_COV_MATRIX) {
                // Free old theoretical covariance
                for (int i = 0; i < nobs; i++) free(covariance_th[i]);
                free(covariance_th);

                // Free old total covariance and inverse
                for (int i = 0; i < nobs; i++) free(covariance_tot[i]);
                free(covariance_tot);
                for (int i = 0; i < nobs; i++) free(inv_cov_tot[i]);
                free(inv_cov_tot);

                // Recalculate theoretical covariance
                get_th_covariance(&covariance_th, names, &nobs, &param);
                get_covtot(&covariance_th, &covariance_exp, &covariance_tot, nobs);

                // Re-invert covariance matrix
                if (get_invcovtot(&covariance_tot, &inv_cov_tot, nobs) < 0) {
                    fprintf(stderr, "Error: Failed to invert covariance matrix at point %d.\n", i);
                    continue; // Skip this point
                }
            }

            double chi2 = get_chi2(inv_cov_tot, predictions, central_exp, nobs);



            // Write results (with predictions if enabled)
            // Only write if chi2 is within threshold
            if (chi2 <= chi2_SM + DELTA_CHI2_MAX) {
                // Write results (with predictions if enabled)
                fprintf(output, FMT_FLOAT ",", parameter_values[0]);
                fprintf(output, FMT_FLOAT, chi2);

                if (record_predictions) {
                    for (int j = 0; j < nobs; j++) {
                        fprintf(output, "," FMT_EXPONENT, predictions[j]);
                        if (record_errors) {
                            double uncertainty;
                            if (RECALCULATE_COV_MATRIX) {
                                uncertainty = sqrt(covariance_th[j][j]);  // Use current covariance
                            } else {
                                uncertainty = sqrt(covariance_th_SM[j][j]);  // Use SM covariance
                            }
                            fprintf(output, "," FMT_EXPONENT, uncertainty);
                        }
                    }
                }
                fprintf(output, "\n");
                fflush(output);
                nSaved++;  // Increment counter
            }

            free(predictions);
            print_progress_bar(i + 1, NUM_POINTS, 1);
        }
    } else {
        // Random scan for multi-dimensional case
        printf("Multiple coefficients detected - performing RANDOM scan over %d coefficients:\n", ncoeff);
        for (int i = 0; i < ncoeff; i++) {
            printf("%s ", coeff_names[i]);
        }
        printf("\n");
        printf("Range: [%.1f, %.1f] with %d points\n", RANGE_MIN, RANGE_MAX, NUM_POINTS);

        for (int i = 0; i < NUM_POINTS; i++) {
            // Generate random parameter values
            for (int j = 0; j < ncoeff; j++) {
                parameter_values[j] = RANGE_MIN + (RANGE_MAX - RANGE_MIN) * rand() / (RAND_MAX + 1.0);
            }

            // Initialize parameters based on scanned values
            read_coeff(parameter_values, &param, coeff_names, &ncoeff);

            // // Get predictions and chi2
            // get_predictions(names, &nobs, &predictions, &param);
            // double chi2 = get_chi2(inv_cov_tot, predictions, central_exp, nobs);

            // Get predictions
            get_predictions(names, &nobs, &predictions, &param);

            // Recalculate covariance matrix if enabled
            if (RECALCULATE_COV_MATRIX) {
                // Free old theoretical covariance
                for (int i = 0; i < nobs; i++) free(covariance_th[i]);
                free(covariance_th);

                // Free old total covariance and inverse
                for (int i = 0; i < nobs; i++) free(covariance_tot[i]);
                free(covariance_tot);
                for (int i = 0; i < nobs; i++) free(inv_cov_tot[i]);
                free(inv_cov_tot);

                // Recalculate theoretical covariance
                get_th_covariance(&covariance_th, names, &nobs, &param);
                get_covtot(&covariance_th, &covariance_exp, &covariance_tot, nobs);

                // Re-invert covariance matrix
                if (get_invcovtot(&covariance_tot, &inv_cov_tot, nobs) < 0) {
                    fprintf(stderr, "Error: Failed to invert covariance matrix at point %d.\n", i);
                    continue; // Skip this point
                }
            }

            double chi2 = get_chi2(inv_cov_tot, predictions, central_exp, nobs);


            // Write results (with predictions if enabled)
            // Only write if chi2 is within threshold
            if (chi2 <= chi2_SM + DELTA_CHI2_MAX) {
                // Write results (with predictions if enabled)
                for (int j = 0; j < ncoeff; j++) {
                    fprintf(output, FMT_FLOAT ",", parameter_values[j]);
                }
                fprintf(output, FMT_FLOAT, chi2);

                if (record_predictions) {
                    for (int j = 0; j < nobs; j++) {
                        fprintf(output, "," FMT_EXPONENT, predictions[j]);
                        if (record_errors) {
                            double uncertainty;
                            if (RECALCULATE_COV_MATRIX) {
                                uncertainty = sqrt(covariance_th[j][j]);  // Use current covariance
                            } else {
                                uncertainty = sqrt(covariance_th_SM[j][j]);  // Use SM covariance
                            }
                            fprintf(output, "," FMT_EXPONENT, uncertainty);
                        }
                    }
                }
                fprintf(output, "\n");
                fflush(output);
                nSaved++;  // Increment counter
            }

            free(predictions);
            print_progress_bar(i + 1, NUM_POINTS, 1);
        }
    }
	fclose(output);

	// // Free allocated memory
	// for (int ie = 0; ie < nobs; ie++) free(covariance_th[ie]);
	// free(covariance_th);
 //
	// for (int ie = 0; ie < nobs; ie++) free(covariance_exp[ie]);
	// free(covariance_exp);
 //
	// free(central_exp);
 //
	// for (int ie = 0; ie < nobs; ie++) free(covariance_tot[ie]);
	// free(covariance_tot);
 //
	// for (int ie = 0; ie < nobs; ie++) free(inv_cov_tot[ie]);
	// free(inv_cov_tot);

    // Free allocated memory
    for (int ie = 0; ie < nobs; ie++) free(covariance_th[ie]);
    free(covariance_th);

    // Free SM covariance if allocated
    if (RECALCULATE_COV_MATRIX && covariance_th_SM) {
        for (int ie = 0; ie < nobs; ie++) free(covariance_th_SM[ie]);
        free(covariance_th_SM);
    }

    for (int ie = 0; ie < nobs; ie++) free(covariance_exp[ie]);
    free(covariance_exp);
    free(central_exp);

    for (int ie = 0; ie < nobs; ie++) free(covariance_tot[ie]);
    free(covariance_tot);

    for (int ie = 0; ie < nobs; ie++) free(inv_cov_tot[ie]);
    free(inv_cov_tot);

	// End measuring time
	clock_t end_time = clock();

	// Calculate the time difference in seconds
	double time_spent = (double)(end_time - start_time) / CLOCKS_PER_SEC;

	// Print the time taken
    printf("\nScan completed. Results written to %s\n", chi2_output_filename);
    printf("Scan took %.1f seconds for %d data points (%d saved with chi2 < ' %.1f + %.1f')\n\n",
        time_spent, NUM_POINTS, nSaved, chi2_SM, DELTA_CHI2_MAX);

	return 0;
}

void print_progress_bar(int iteration, int total, int show_time) {
    int barWidth = 50;
    float progress = (float)iteration / total;

    // Initialize start time on first call
    if (is_first_progress) {
        gettimeofday(&start_time, NULL);
        is_first_progress = 0;
    }

    printf("[");
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) printf("=");
        else if (i == pos) printf(">");
        else printf(" ");
    }

    // Calculate time information if requested
    if (show_time) {
        struct timeval current_time;
        gettimeofday(&current_time, NULL);
        double elapsed = (current_time.tv_sec - start_time.tv_sec) +
                        (current_time.tv_usec - start_time.tv_usec) / 1000000.0;
        double remaining = (total - iteration) * (elapsed / iteration);

        printf("] %d%% (%.1fs/%.1fs)",
              (int)(progress * 100),
              elapsed,
              remaining);
    } else {
        printf("] %d%%", (int)(progress * 100));
    }

    printf("\r");
    fflush(stdout);

    // Print newline when complete
    if (iteration == total) {
        printf("\n");
        is_first_progress = 1; // Reset for next use
    }
}

// Check if a file is non-empty
int is_file_non_empty(const char* filename) {
	struct stat file_stat;
	if (stat(filename, &file_stat) != 0) {
		// File doesn't exist or can't be accessed
		return 0;
	}
	return file_stat.st_size > 0;
}


void apply_coeff_value(struct parameters *param, const char* name, double value) {
    // Check if coefficient is flavor-specific (contains "_mu_", "_e_", or "_tau_")
    int is_flavor_specific = (strstr(name, "_mu_") || strstr(name, "_e_") || strstr(name, "_tau_"));

    /* ***** Handle B-decays (b->s) coefficients ***** */
    // Handle real and non-primed
    if (strstr(name, "_real") && !strstr(name, "_prime")) {
        if (strstr(name, "dC7_")) {
            param->deltaC[7] = value + cimag(param->deltaC[7]) * I;
        }
        else if (strstr(name, "dC8_")) {
            param->deltaC[8] = value + cimag(param->deltaC[8]) * I;
        }
        else if (strstr(name, "dC9_")) {
            if (!is_flavor_specific) {
                // Universal: set for all flavors
                param->deltaC[9]  = value + cimag(param->deltaC[9]) * I;   // mu
                param->deltaC[19] = value + cimag(param->deltaC[19]) * I;  // e
                param->deltaC[29] = value + cimag(param->deltaC[29]) * I;  // tau
            }
            else if (strstr(name, "_mu_")) {
                param->deltaC[9] = value + cimag(param->deltaC[9]) * I;
            }
            else if (strstr(name, "_e_")) {
                param->deltaC[19] = value + cimag(param->deltaC[19]) * I;
            }
            else if (strstr(name, "_tau_")) {
                param->deltaC[29] = value + cimag(param->deltaC[29]) * I;
            }
        }
        else if (strstr(name, "dC10_")) {
            if (!is_flavor_specific) {
                // Universal: set for all flavors
                param->deltaC[10]  = value + cimag(param->deltaC[10]) * I;   // mu
                param->deltaC[20] = value + cimag(param->deltaC[20]) * I;  // e
                param->deltaC[30] = value + cimag(param->deltaC[30]) * I;  // tau
            }
            else if (strstr(name, "_mu_")) {
                param->deltaC[10] = value + cimag(param->deltaC[10]) * I;
            }
            else if (strstr(name, "_e_")) {
                param->deltaC[20] = value + cimag(param->deltaC[20]) * I;
            }
            else if (strstr(name, "_tau_")) {
                param->deltaC[30] = value + cimag(param->deltaC[30]) * I;
            }
        }
        else if (strstr(name, "dCQ1_")) {
            if (!is_flavor_specific) {
                // Universal: set for all flavors
                param->deltaCQ[1] = value + cimag(param->deltaCQ[1]) * I;   // mu
                param->deltaCQ[3] = value + cimag(param->deltaC[3]) * I;  // e
                param->deltaCQ[5] = value + cimag(param->deltaC[5]) * I;  // tau
            }
            else if (strstr(name, "_mu_")) {
                param->deltaCQ[1] = value + cimag(param->deltaCQ[1]) * I;
            }
            else if (strstr(name, "_e_")) {
                param->deltaCQ[3] = value + cimag(param->deltaCQ[3]) * I;
            }
            else if (strstr(name, "_tau_")) {
                param->deltaCQ[5] = value + cimag(param->deltaCQ[5]) * I;
            }
        }
        else if (strstr(name, "dCQ2_")) {
            if (!is_flavor_specific) {
                // Universal: set for all flavors
                param->deltaCQ[2] = value + cimag(param->deltaCQ[2]) * I;   // mu
                param->deltaCQ[4] = value + cimag(param->deltaCQ[4]) * I;  // e
                param->deltaCQ[6] = value + cimag(param->deltaCQ[6]) * I;  // tau
            }
            else if (strstr(name, "_mu_")) {
                param->deltaCQ[2] = value + cimag(param->deltaCQ[2]) * I;
            }
            else if (strstr(name, "_e_")) {
                param->deltaCQ[4] = value + cimag(param->deltaCQ[4]) * I;
            }
            else if (strstr(name, "_tau_")) {
                param->deltaCQ[6] = value + cimag(param->deltaCQ[6]) * I;
            }
        }
    }
    // Handle real and primed
    if (strstr(name, "_real") && strstr(name, "_prime")) {
        if (strstr(name, "dC7_")) {
            param->deltaCp[7] = value + cimag(param->deltaCp[7]) * I;
        }
        else if (strstr(name, "dC8_")) {
            param->deltaCp[8] = value + cimag(param->deltaCp[8]) * I;
        }
        else if (strstr(name, "dC9_")) {
            if (!is_flavor_specific) {
                // Universal: set for all flavors
                param->deltaCp[9]  = value + cimag(param->deltaCp[9]) * I;   // mu
                param->deltaCp[19] = value + cimag(param->deltaCp[19]) * I;  // e
                param->deltaCp[29] = value + cimag(param->deltaCp[29]) * I;  // tau
            }
            else if (strstr(name, "_mu_")) {
                param->deltaCp[9] = value + cimag(param->deltaCp[9]) * I;
            }
            else if (strstr(name, "_e_")) {
                param->deltaCp[19] = value + cimag(param->deltaCp[19]) * I;
            }
            else if (strstr(name, "_tau_")) {
                param->deltaCp[29] = value + cimag(param->deltaCp[29]) * I;
            }
        }
        else if (strstr(name, "dC10_")) {
            if (!is_flavor_specific) {
                // Universal: set for all flavors
                param->deltaCp[10]  = value + cimag(param->deltaCp[10]) * I;   // mu
                param->deltaCp[20] = value + cimag(param->deltaCp[20]) * I;  // e
                param->deltaCp[30] = value + cimag(param->deltaCp[30]) * I;  // tau
            }
            else if (strstr(name, "_mu_")) {
                param->deltaCp[10] = value + cimag(param->deltaCp[10]) * I;
            }
            else if (strstr(name, "_e_")) {
                param->deltaCp[20] = value + cimag(param->deltaCp[20]) * I;
            }
            else if (strstr(name, "_tau_")) {
                param->deltaCp[30] = value + cimag(param->deltaCp[30]) * I;
            }
        }
        else if (strstr(name, "dCQ1_")) {
            if (!is_flavor_specific) {
                // Universal: set for all flavors
                param->deltaCQp[1] = value + cimag(param->deltaCQp[1]) * I;   // mu
                param->deltaCQp[3] = value + cimag(param->deltaCQp[3]) * I;  // e
                param->deltaCQp[5] = value + cimag(param->deltaCQp[5]) * I;  // tau
            }
            else if (strstr(name, "_mu_")) {
                param->deltaCQp[1] = value + cimag(param->deltaCQp[1]) * I;
            }
            else if (strstr(name, "_e_")) {
                param->deltaCQp[3] = value + cimag(param->deltaCQp[3]) * I;
            }
            else if (strstr(name, "_tau_")) {
                param->deltaCQp[5] = value + cimag(param->deltaCQp[5]) * I;
            }
        }
        else if (strstr(name, "dCQ2_")) {
            if (!is_flavor_specific) {
                // Universal: set for all flavors
                param->deltaCQp[2] = value + cimag(param->deltaCQp[2]) * I;   // mu
                param->deltaCQp[4] = value + cimag(param->deltaCQp[4]) * I;  // e
                param->deltaCQp[6] = value + cimag(param->deltaCQp[6]) * I;  // tau
            }
            else if (strstr(name, "_mu_")) {
                param->deltaCQp[2] = value + cimag(param->deltaCQp[2]) * I;
            }
            else if (strstr(name, "_e_")) {
                param->deltaCQp[4] = value + cimag(param->deltaCQp[4]) * I;
            }
            else if (strstr(name, "_tau_")) {
                param->deltaCQp[6] = value + cimag(param->deltaCQp[6]) * I;
            }
        }
    }
    // Handle imaginary parts (same logic as above but with _imag)
    else if (strstr(name, "_imag") && !strstr(name, "_prime")) {
        if (strstr(name, "dC7_")) {
            param->deltaC[7] = creal(param->deltaC[7]) + value * I;
        }
        else if (strstr(name, "dC8_")) {
            param->deltaC[8] = creal(param->deltaC[8]) + value * I;
        }
        else if (strstr(name, "dC9_")) {
            if (!is_flavor_specific) {
                // Universal: set for all flavors
                param->deltaC[9]  = creal(param->deltaC[9]) + value * I;   // mu
                param->deltaC[19] = creal(param->deltaC[19]) + value * I;  // e
                param->deltaC[29] = creal(param->deltaC[29]) + value * I;  // tau
            }
            else if (strstr(name, "_mu_")) {
                param->deltaC[9] = creal(param->deltaC[9]) + value * I;
            }
            else if (strstr(name, "_e_")) {
                param->deltaC[19] = creal(param->deltaC[19]) + value * I;
            }
            else if (strstr(name, "_tau_")) {
                param->deltaC[29] = creal(param->deltaC[29]) + value * I;
            }
        }
        else if (strstr(name, "dC10_")) {
            if (!is_flavor_specific) {
                // Universal: set for all flavors
                param->deltaC[10]  = creal(param->deltaC[10]) + value * I;   // mu
                param->deltaC[20] = creal(param->deltaC[20]) + value * I;  // e
                param->deltaC[30] = creal(param->deltaC[30]) + value * I;  // tau
            }
            else if (strstr(name, "_mu_")) {
                param->deltaC[10] = creal(param->deltaC[10]) + value * I;
            }
            else if (strstr(name, "_e_")) {
                param->deltaC[20] = creal(param->deltaC[20]) + value * I;
            }
            else if (strstr(name, "_tau_")) {
                param->deltaC[30] = creal(param->deltaC[30]) + value * I;
            }
        }
        else if (strstr(name, "dCQ1_")) {
            if (!is_flavor_specific) {
                // Universal: set for all flavors
                param->deltaCQ[1] = creal(param->deltaCQ[1]) + value * I;   // mu
                param->deltaCQ[3] = creal(param->deltaC[3]) + value * I;  // e
                param->deltaCQ[5] = creal(param->deltaC[5]) + value * I;  // tau
            }
            else if (strstr(name, "_mu_")) {
                param->deltaCQ[1] = creal(param->deltaCQ[1]) + value * I;
            }
            else if (strstr(name, "_e_")) {
                param->deltaCQ[3] = creal(param->deltaCQ[3]) + value * I;
            }
            else if (strstr(name, "_tau_")) {
                param->deltaCQ[5] = creal(param->deltaCQ[5]) + value * I;
            }
        }
        else if (strstr(name, "dCQ2_")) {
            if (!is_flavor_specific) {
                // Universal: set for all flavors
                param->deltaCQ[2] = creal(param->deltaCQ[2]) + value * I;   // mu
                param->deltaCQ[4] = creal(param->deltaCQ[4]) + value * I;  // e
                param->deltaCQ[6] = creal(param->deltaCQ[6]) + value * I;  // tau
            }
            else if (strstr(name, "_mu_")) {
                param->deltaCQ[2] = creal(param->deltaCQ[2]) + value * I;
            }
            else if (strstr(name, "_e_")) {
                param->deltaCQ[4] = creal(param->deltaCQ[4]) + value * I;
            }
            else if (strstr(name, "_tau_")) {
                param->deltaCQ[6] = creal(param->deltaCQ[6]) + value * I;
            }
        }
    }
    else if (strstr(name, "_imag") && strstr(name, "_prime")) {
        // Handle primed coefficients (similar logic)
        if (strstr(name, "dC7_prime_")) {
            param->deltaCp[7] = creal(param->deltaCp[7]) + value * I;
        }
        else if (strstr(name, "dC8_")) {
            param->deltaCp[8] = creal(param->deltaCp[8]) + value * I;
        }
        else if (strstr(name, "dC9_")) {
            if (!is_flavor_specific) {
                // Universal: set for all flavors
                param->deltaCp[9]  = creal(param->deltaCp[9]) + value * I;   // mu
                param->deltaCp[19] = creal(param->deltaCp[19]) + value * I;  // e
                param->deltaCp[29] = creal(param->deltaCp[29]) + value * I;  // tau
            }
            else if (strstr(name, "_mu_")) {
                param->deltaCp[9] = creal(param->deltaCp[9]) + value * I;
            }
            else if (strstr(name, "_e_")) {
                param->deltaCp[19] = creal(param->deltaCp[19]) + value * I;
            }
            else if (strstr(name, "_tau_")) {
                param->deltaCp[29] = creal(param->deltaCp[29]) + value * I;
            }
        }
        else if (strstr(name, "dC10_")) {
            if (!is_flavor_specific) {
                // Universal: set for all flavors
                param->deltaCp[10]  = creal(param->deltaCp[10]) + value * I;   // mu
                param->deltaCp[20] = creal(param->deltaCp[20]) + value * I;  // e
                param->deltaCp[30] = creal(param->deltaCp[30]) + value * I;  // tau
            }
            else if (strstr(name, "_mu_")) {
                param->deltaCp[10] = creal(param->deltaCp[10]) + value * I;
            }
            else if (strstr(name, "_e_")) {
                param->deltaCp[20] = creal(param->deltaCp[20]) + value * I;
            }
            else if (strstr(name, "_tau_")) {
                param->deltaCp[30] = creal(param->deltaCp[30]) + value * I;
            }
        }
        else if (strstr(name, "dCQ1_")) {
            if (!is_flavor_specific) {
                // Universal: set for all flavors
                param->deltaCQp[1] = creal(param->deltaCQp[1]) + value * I;   // mu
                param->deltaCQp[3] = creal(param->deltaCQp[3]) + value * I;  // e
                param->deltaCQp[5] = creal(param->deltaCQp[5]) + value * I;  // tau
            }
            else if (strstr(name, "_mu_")) {
                param->deltaCQp[1] = creal(param->deltaCQp[1]) + value * I;
            }
            else if (strstr(name, "_e_")) {
                param->deltaCQp[3] = creal(param->deltaCQp[3]) + value * I;
            }
            else if (strstr(name, "_tau_")) {
                param->deltaCQp[5] = creal(param->deltaCQp[5]) + value * I;
            }
        }
        else if (strstr(name, "dCQ2_")) {
            if (!is_flavor_specific) {
                // Universal: set for all flavors
                param->deltaCQp[2] = creal(param->deltaCQp[2]) + value * I;   // mu
                param->deltaCQp[4] = creal(param->deltaCQp[4]) + value * I;  // e
                param->deltaCQp[6] = creal(param->deltaCQp[6]) + value * I;  // tau
            }
            else if (strstr(name, "_mu_")) {
                param->deltaCQp[2] = creal(param->deltaCQp[2]) + value * I;
            }
            else if (strstr(name, "_e_")) {
                param->deltaCQp[4] = creal(param->deltaCQp[4]) + value * I;
            }
            else if (strstr(name, "_tau_")) {
                param->deltaCQp[6] = creal(param->deltaCQp[6]) + value * I;
            }
        }
    }


    /* **** Handle kaon (s->d) coefficients **** */
    else if (strstr(name, "dKC")) {
        // Check if coefficient is flavor-specific (contains "_e", "_mu", or "_tau")
        int is_flavor_specific = (strstr(name, "_e_") || strstr(name, "_mu_") || strstr(name, "_tau_"));

        // Handle real and non-primed
        if (strstr(name, "_real") && !strstr(name, "_prime")) {
            if (strstr(name, "dKC9_")) {
                if (!is_flavor_specific) {
                    // Universal: set for both flavors
                    param->KC9e = value + cimag(param->KC9e) * I;
                    param->KC9mu = value + cimag(param->KC9mu) * I;
                }
                else if (strstr(name, "_e_")) {
                    param->KC9e = value + cimag(param->KC9e) * I;
                }
                else if (strstr(name, "_mu_")) {
                    param->KC9mu = value + cimag(param->KC9mu) * I;
                }
            }
            else if (strstr(name, "dKC10_")) {
                if (!is_flavor_specific) {
                    // Universal: set for both flavors
                    param->KC10e = value + cimag(param->KC10e) * I;
                    param->KC10mu = value + cimag(param->KC10mu) * I;
                }
                else if (strstr(name, "_e_")) {
                    param->KC10e = value + cimag(param->KC10e) * I;
                }
                else if (strstr(name, "_mu_")) {
                    param->KC10mu = value + cimag(param->KC10mu) * I;
                }
            }
            else if (strstr(name, "dKCQ1_")) {
                if (!is_flavor_specific) {
                    // Universal: set for both flavors
                    param->KCQ1e = value + cimag(param->KCQ1e) * I;
                    param->KCQ1mu = value + cimag(param->KCQ1mu) * I;
                }
                else if (strstr(name, "_e_")) {
                    param->KCQ1e = value + cimag(param->KCQ1e) * I;
                }
                else if (strstr(name, "_mu_")) {
                    param->KCQ1mu = value + cimag(param->KCQ1mu) * I;
                }
            }
            else if (strstr(name, "dKCQ2_")) {
                if (!is_flavor_specific) {
                    // Universal: set for both flavors
                    param->KCQ2e = value + cimag(param->KCQ2e) * I;
                    param->KCQ2mu = value + cimag(param->KCQ2mu) * I;
                }
                else if (strstr(name, "_e_")) {
                    param->KCQ2e = value + cimag(param->KCQ2e) * I;
                }
                else if (strstr(name, "_mu_")) {
                    param->KCQ2mu = value + cimag(param->KCQ2mu) * I;
                }
            }
            else if (strstr(name, "dKCL_")) {
                if (!is_flavor_specific) {
                    // Universal: set for all flavors
                    param->KCLe = value + cimag(param->KCLe) * I;
                    param->KCLmu = value + cimag(param->KCLmu) * I;
                    param->KCLtau = value + cimag(param->KCLtau) * I;
                }
                else if (strstr(name, "_e_")) {
                    param->KCLe = value + cimag(param->KCLe) * I;
                }
                else if (strstr(name, "_mu_")) {
                    param->KCLmu = value + cimag(param->KCLmu) * I;
                }
                else if (strstr(name, "_tau_")) {
                    param->KCLtau = value + cimag(param->KCLtau) * I;
                }
            }
            else if (strstr(name, "dKCR_")) {
                if (!is_flavor_specific) {
                    // Universal: set for all flavors
                    param->KCRe = value + cimag(param->KCRe) * I;
                    param->KCRmu = value + cimag(param->KCRmu) * I;
                    param->KCRtau = value + cimag(param->KCRtau) * I;
                }
                else if (strstr(name, "_e_")) {
                    param->KCRe = value + cimag(param->KCRe) * I;
                }
                else if (strstr(name, "_mu_")) {
                    param->KCRmu = value + cimag(param->KCRmu) * I;
                }
                else if (strstr(name, "_tau_")) {
                    param->KCRtau = value + cimag(param->KCRtau) * I;
                }
            }
        }
        // Handle real and primed
        else if (strstr(name, "_real") && strstr(name, "_prime")) {
            if (strstr(name, "dKC9_")) {
                if (!is_flavor_specific) {
                    // Universal: set for both flavors
                    param->KC9pe = value + cimag(param->KC9pe) * I;
                    param->KC9pmu = value + cimag(param->KC9pmu) * I;
                }
                else if (strstr(name, "_e_")) {
                    param->KC9pe = value + cimag(param->KC9pe) * I;
                }
                else if (strstr(name, "_mu_")) {
                    param->KC9pmu = value + cimag(param->KC9pmu) * I;
                }
            }
            else if (strstr(name, "dKC10_")) {
                if (!is_flavor_specific) {
                    // Universal: set for both flavors
                    param->KC10pe = value + cimag(param->KC10pe) * I;
                    param->KC10pmu = value + cimag(param->KC10pmu) * I;
                }
                else if (strstr(name, "_e_")) {
                    param->KC10pe = value + cimag(param->KC10pe) * I;
                }
                else if (strstr(name, "_mu_")) {
                    param->KC10pmu = value + cimag(param->KC10pmu) * I;
                }
            }
            else if (strstr(name, "dKCQ1_")) {
                if (!is_flavor_specific) {
                    // Universal: set for both flavors
                    param->KCQ1pe = value + cimag(param->KCQ1pe) * I;
                    param->KCQ1pmu = value + cimag(param->KCQ1pmu) * I;
                }
                else if (strstr(name, "_e_")) {
                    param->KCQ1pe = value + cimag(param->KCQ1pe) * I;
                }
                else if (strstr(name, "_mu_")) {
                    param->KCQ1pmu = value + cimag(param->KCQ1pmu) * I;
                }
            }
            else if (strstr(name, "dKCQ2_")) {
                if (!is_flavor_specific) {
                    // Universal: set for both flavors
                    param->KCQ2pe = value + cimag(param->KCQ2pe) * I;
                    param->KCQ2pmu = value + cimag(param->KCQ2pmu) * I;
                }
                else if (strstr(name, "_e_")) {
                    param->KCQ2pe = value + cimag(param->KCQ2pe) * I;
                }
                else if (strstr(name, "_mu_")) {
                    param->KCQ2pmu = value + cimag(param->KCQ2pmu) * I;
                }
            }
        }
        // Handle imaginary parts (same logic as above but with _imag)
        else if (strstr(name, "_imag") && !strstr(name, "_prime")) {
            if (strstr(name, "dKC9_")) {
                if (!is_flavor_specific) {
                    // Universal: set for both flavors
                    param->KC9e = creal(param->KC9e) + value * I;
                    param->KC9mu = creal(param->KC9mu) + value * I;
                }
                else if (strstr(name, "_e_")) {
                    param->KC9e = creal(param->KC9e) + value * I;
                }
                else if (strstr(name, "_mu_")) {
                    param->KC9mu = creal(param->KC9mu) + value * I;
                }
            }
            else if (strstr(name, "dKC10_")) {
                if (!is_flavor_specific) {
                    // Universal: set for both flavors
                    param->KC10e = creal(param->KC10e) + value * I;
                    param->KC10mu = creal(param->KC10mu) + value * I;
                }
                else if (strstr(name, "_e_")) {
                    param->KC10e = creal(param->KC10e) + value * I;
                }
                else if (strstr(name, "_mu_")) {
                    param->KC10mu = creal(param->KC10mu) + value * I;
                }
            }
            else if (strstr(name, "dKCQ1_")) {
                if (!is_flavor_specific) {
                    // Universal: set for both flavors
                    param->KCQ1e = creal(param->KCQ1e) + value * I;
                    param->KCQ1mu = creal(param->KCQ1mu) + value * I;
                }
                else if (strstr(name, "_e_")) {
                    param->KCQ1e = creal(param->KCQ1e) + value * I;
                }
                else if (strstr(name, "_mu_")) {
                    param->KCQ1mu = creal(param->KCQ1mu) + value * I;
                }
            }
            else if (strstr(name, "dKCQ2_")) {
                if (!is_flavor_specific) {
                    // Universal: set for both flavors
                    param->KCQ2e = creal(param->KCQ2e) + value * I;
                    param->KCQ2mu = creal(param->KCQ2mu) + value * I;
                }
                else if (strstr(name, "_e_")) {
                    param->KCQ2e = creal(param->KCQ2e) + value * I;
                }
                else if (strstr(name, "_mu_")) {
                    param->KCQ2mu = creal(param->KCQ2mu) + value * I;
                }
            }
            else if (strstr(name, "dKCL_")) {
                if (!is_flavor_specific) {
                    // Universal: set for all flavors
                    param->KCLe = creal(param->KCLe) + value * I;
                    param->KCLmu = creal(param->KCLmu) + value * I;
                    param->KCLtau = creal(param->KCLtau) + value * I;
                }
                else if (strstr(name, "_e_")) {
                    param->KCLe = creal(param->KCLe) + value * I;
                }
                else if (strstr(name, "_mu_")) {
                    param->KCLmu = creal(param->KCLmu) + value * I;
                }
                else if (strstr(name, "_tau_")) {
                    param->KCLtau = creal(param->KCLtau) + value * I;
                }
            }
            else if (strstr(name, "dKCR_")) {
                if (!is_flavor_specific) {
                    // Universal: set for all flavors
                    param->KCRe = creal(param->KCRe) + value * I;
                    param->KCRmu = creal(param->KCRmu) + value * I;
                    param->KCRtau = creal(param->KCRtau) + value * I;
                }
                else if (strstr(name, "_e_")) {
                    param->KCRe = creal(param->KCRe) + value * I;
                }
                else if (strstr(name, "_mu_")) {
                    param->KCRmu = creal(param->KCRmu) + value * I;
                }
                else if (strstr(name, "_tau_")) {
                    param->KCRtau = creal(param->KCRtau) + value * I;
                }
            }
        }
        else if (strstr(name, "_imag") && strstr(name, "_prime")) {
            // Handle primed coefficients imaginary parts
            if (strstr(name, "dKC9_")) {
                if (!is_flavor_specific) {
                    // Universal: set for both flavors
                    param->KC9pe = creal(param->KC9pe) + value * I;
                    param->KC9pmu = creal(param->KC9pmu) + value * I;
                }
                else if (strstr(name, "_e_")) {
                    param->KC9pe = creal(param->KC9pe) + value * I;
                }
                else if (strstr(name, "_mu_")) {
                    param->KC9pmu = creal(param->KC9pmu) + value * I;
                }
            }
            else if (strstr(name, "dKC10_")) {
                if (!is_flavor_specific) {
                    // Universal: set for both flavors
                    param->KC10pe = creal(param->KC10pe) + value * I;
                    param->KC10pmu = creal(param->KC10pmu) + value * I;
                }
                else if (strstr(name, "_e_")) {
                    param->KC10pe = creal(param->KC10pe) + value * I;
                }
                else if (strstr(name, "_mu_")) {
                    param->KC10pmu = creal(param->KC10pmu) + value * I;
                }
            }
            else if (strstr(name, "dKCQ1_")) {
                if (!is_flavor_specific) {
                    // Universal: set for both flavors
                    param->KCQ1pe = creal(param->KCQ1pe) + value * I;
                    param->KCQ1pmu = creal(param->KCQ1pmu) + value * I;
                }
                else if (strstr(name, "_e_")) {
                    param->KCQ1pe = creal(param->KCQ1pe) + value * I;
                }
                else if (strstr(name, "_mu_")) {
                    param->KCQ1pmu = creal(param->KCQ1pmu) + value * I;
                }
            }
            else if (strstr(name, "dKCQ2_")) {
                if (!is_flavor_specific) {
                    // Universal: set for both flavors
                    param->KCQ2pe = creal(param->KCQ2pe) + value * I;
                    param->KCQ2pmu = creal(param->KCQ2pmu) + value * I;
                }
                else if (strstr(name, "_e_")) {
                    param->KCQ2pe = creal(param->KCQ2pe) + value * I;
                }
                else if (strstr(name, "_mu_")) {
                    param->KCQ2pmu = creal(param->KCQ2pmu) + value * I;
                }
            }
        }
    }

}


// Initialize parameters
int read_coeff(double parameter_values[], struct parameters *param, char coeff_names[][50], int *ncoeff_ptr) {
    for (int i = 0; i < *ncoeff_ptr; i++) {
        apply_coeff_value(param, coeff_names[i], parameter_values[i]);
    }
    return 0;
}

    // void write_header(FILE *output_file, char coeff_names[][50], int ncoeffs) {
    //     for (int i = 0; i < ncoeffs; i++) {
    //         fprintf(output_file, "%s,", coeff_names[i]);
    //     }
    //     fprintf(output_file, "chi2\n");
    // }

void print_current_time_in_timezone() {
    time_t rawtime;
    struct tm *timeinfo;

    // Get the current system time
    time(&rawtime);

    // Convert the time to local time
    timeinfo = localtime(&rawtime);

    // Print the current time in a human-readable format
    printf("\nStarted at: %s", asctime(timeinfo));
}
// Generate datetime string in compact YYYYMMDD_HHMMSS format
void get_datetime_string(char *buffer, size_t size) {
    time_t now = time(NULL);
    struct tm *tm_info = localtime(&now);
    strftime(buffer, size, "%Y%m%d_%H%M%S", tm_info);  // e.g. "20240425_142305"
}

// // Rename the function
// void write_output_header(FILE *output_file, char coeff_names[][50], int ncoeffs,
//                         char obs_names[][50], int nobs, int record_predictions) {
//     // Write coefficient headers
//     for (int i = 0; i < ncoeffs; i++) {
//         fprintf(output_file, "%s,", coeff_names[i]);
//     }
//
//     // Write chi2 header
//     fprintf(output_file, "chi2");
//
//     // Write observable headers if predictions are recorded
//     if (record_predictions) {
//         for (int i = 0; i < nobs; i++) {
//             fprintf(output_file, ",%s", obs_names[i]);
//         }
//     }
//     fprintf(output_file, "\n");
// }
