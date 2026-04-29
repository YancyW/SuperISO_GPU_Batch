#include "src/include.h"

#define WORKER_LINE_MAX 4096

static int evaluate_slha_file(FILE *out, const char *name)
{
	int test = test_slha((char *)name);

	fprintf(out, "\n");
	fprintf(out, "SuperIso " VERSION " - %d\n\n", YEAR);
	fprintf(out, "SLHA input file\n\n");

	if(test > 0)
	{
		if(test == 2) fprintf(out, "WARNING: only tested in the MFV scenario!\n\n");
	}
	else
	{
		if(test == -1) fprintf(out, "Invalid point\n\n");
		else if(test == -2) fprintf(out, "Model not yet implemented\n\n");
		else if(test == -3) fprintf(out, "Invalid SLHA file\n\n");
		else if(test == -4) fprintf(out, "SLHA file absent\n\n");
		return 1;
	}

	char obsfile[500], expfile[500], corrfile[500], nuisfile[500], nuiscorrfile[500];
	sprintf(obsfile, "chi2_input/myobs.in");
	sprintf(expfile, "chi2_input/exp_values.in");
	sprintf(corrfile, "chi2_input/exp_corr.in");
	sprintf(nuisfile, "chi2_input/nuisance.in");
	sprintf(nuiscorrfile, "chi2_input/nuisance_corr.in");

	struct parameters param;
	int ie = 0;
	char names[NBOBSMAX][50];

	int nobs = read_obs_list(obsfile, names);
	fprintf(out, "Reading observable list from %s\n", obsfile);

	Init_param(&param);
	Les_Houches_Reader((char *)name, &param);

	sprintf(param.exp_values, "%s", expfile);
	sprintf(param.exp_corr, "%s", corrfile);
	sprintf(param.nuisance_values, "%s", nuisfile);
	sprintf(param.nuisance_corr, "%s", nuiscorrfile);

	fprintf(out, "\n");

	double *predictions;
	get_predictions(names, &nobs, &predictions, &param);

	double **covariance_th;
	get_th_covariance(&covariance_th, names, &nobs, &param);

	double *central_exp, **covariance_exp;
	get_exp_covariance(&covariance_exp, &central_exp, names, &nobs, &param);

	for(ie = 0; ie < nobs; ie++)
		fprintf(out, "\t%-25s =\t%.5e +/- %.5e\t(exp: %.5e +/- %.5e)\n",
			names[ie], predictions[ie], sqrt(covariance_th[ie][ie]), central_exp[ie], sqrt(covariance_exp[ie][ie]));

	double **covariance_tot;
	get_covtot(&covariance_th, &covariance_exp, &covariance_tot, nobs);

	double **inv_cov_tot;
	if(get_invcovtot(&covariance_tot, &inv_cov_tot, nobs) < 0)
	{
		fprintf(out, "Invalid total covariance matrix\n");
		return 1;
	}

	double chi2 = get_chi2(inv_cov_tot, predictions, central_exp, nobs);
	fprintf(out, "\nchi2 = %.2f\n", chi2);
	fprintf(out, "n_obs=%d\n\n", nobs);

	write_covariance("slha_covariance.out", covariance_tot, names, nobs);
	fprintf(out, "Total covariance matrix written in slha_covariance.out\n");

	for(ie = 0; ie < nobs; ie++) free(covariance_th[ie]);
	free(covariance_th);

	for(ie = 0; ie < nobs; ie++) free(covariance_exp[ie]);
	free(covariance_exp);

	free(central_exp);

	for(ie = 0; ie < nobs; ie++) free(covariance_tot[ie]);
	free(covariance_tot);

	for(ie = 0; ie < nobs; ie++) free(inv_cov_tot[ie]);
	free(inv_cov_tot);

	return 0;
}

static void trim_newline(char *text)
{
	size_t len = strlen(text);
	while(len > 0 && (text[len - 1] == '\n' || text[len - 1] == '\r'))
	{
		text[len - 1] = '\0';
		len--;
	}
}

static int process_worker_request(const char *request)
{
	char buffer[WORKER_LINE_MAX];
	char *tab = NULL;
	char *request_id = NULL;
	char *file_name = NULL;
	int rc = 1;

	strncpy(buffer, request, sizeof(buffer) - 1);
	buffer[sizeof(buffer) - 1] = '\0';
	trim_newline(buffer);

	tab = strchr(buffer, '\t');
	if(tab)
	{
		*tab = '\0';
		request_id = buffer;
		file_name = tab + 1;
	}
	else
	{
		request_id = "request";
		file_name = buffer;
	}

	printf("===BEGIN\t%s===\n", request_id);
	fflush(stdout);
	rc = evaluate_slha_file(stdout, file_name);
	printf("===END\t%s\t%d===\n", request_id, rc);
	fflush(stdout);
	return rc;
}

int main(int argc, char **argv)
{
	if(argc == 2)
	{
		return evaluate_slha_file(stdout, argv[1]);
	}

	if(argc != 1)
	{
		fprintf(stderr, "Usage: %s [slha_file]\n", argv[0]);
		fprintf(stderr, "Or run without arguments and stream '<id>\\t<slha_file>' requests on stdin.\n");
		return 1;
	}

	char line[WORKER_LINE_MAX];
	while(fgets(line, sizeof(line), stdin) != NULL)
	{
		if(line[0] == '\n' || line[0] == '\r' || line[0] == '\0') continue;
		process_worker_request(line);
	}

	return 0;
}