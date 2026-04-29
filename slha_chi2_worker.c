#include "src/include.h"

#define WORKER_LINE_MAX 4096

typedef struct worker_context
{
	char obsfile[500];
	char expfile[500];
	char corrfile[500];
	char nuisfile[500];
	char nuiscorrfile[500];
	char names[NBOBSMAX][50];
	char namenuisance[NNUIS + 1][50];
	int nobs;
	struct nuisance nuisparam;
	double *central_exp;
	double **covariance_exp;
	double **nuiscorr;
} worker_context;

static void free_matrix(double **matrix, int rows)
{
	int ie;
	if(!matrix) return;
	for(ie = 0; ie < rows; ie++) free(matrix[ie]);
	free(matrix);
}

static void free_worker_context(worker_context *ctx)
{
	if(!ctx) return;
	free(ctx->central_exp);
	ctx->central_exp = NULL;
	free_matrix(ctx->covariance_exp, ctx->nobs);
	ctx->covariance_exp = NULL;
	free_matrix(ctx->nuiscorr, NNUIS + 1);
	ctx->nuiscorr = NULL;
}

static int init_worker_context(worker_context *ctx)
{
	int ie, je;
	struct parameters param;

	memset(ctx, 0, sizeof(*ctx));

	sprintf(ctx->obsfile, "chi2_input/myobs.in");
	sprintf(ctx->expfile, "chi2_input/exp_values.in");
	sprintf(ctx->corrfile, "chi2_input/exp_corr.in");
	sprintf(ctx->nuisfile, "chi2_input/nuisance.in");
	sprintf(ctx->nuiscorrfile, "chi2_input/nuisance_corr.in");

	ctx->nobs = read_obs_list(ctx->obsfile, ctx->names);

	Init_param(&param);
	sprintf(param.exp_values, "%s", ctx->expfile);
	sprintf(param.exp_corr, "%s", ctx->corrfile);
	sprintf(param.nuisance_values, "%s", ctx->nuisfile);
	sprintf(param.nuisance_corr, "%s", ctx->nuiscorrfile);

	read_nuisance(param.nuisance_values, &ctx->nuisparam);
	observables(0, NULL, 0, NULL, NULL, &ctx->nuisparam, ctx->namenuisance, &param);

	ctx->nuiscorr = (double **) malloc((NNUIS + 1) * sizeof(double *));
	for(ie = 0; ie <= NNUIS; ie++) ctx->nuiscorr[ie] = (double *) malloc((NNUIS + 1) * sizeof(double));
	for(ie = 0; ie <= NNUIS; ie++) for(je = 0; je <= NNUIS; je++) ctx->nuiscorr[ie][je] = (ie == je);
	read_correlation(param.nuisance_corr, ctx->nuiscorr, ctx->namenuisance, NNUIS);

	get_exp_covariance(&ctx->covariance_exp, &ctx->central_exp, ctx->names, &ctx->nobs, &param);
	return 1;
}

static int evaluate_slha_file(FILE *out, const char *name, worker_context *ctx)
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

	struct parameters param;
	int ie = 0;
	int nobs = ctx->nobs;
	fprintf(out, "Reading observable list from %s\n", ctx->obsfile);

	Init_param(&param);
	Les_Houches_Reader((char *)name, &param);

	sprintf(param.exp_values, "%s", ctx->expfile);
	sprintf(param.exp_corr, "%s", ctx->corrfile);
	sprintf(param.nuisance_values, "%s", ctx->nuisfile);
	sprintf(param.nuisance_corr, "%s", ctx->nuiscorrfile);

	fprintf(out, "\n");

	double *predictions;
	get_predictions_nuisance(ctx->names, &nobs, &predictions, &param, &ctx->nuisparam);

	double **covariance_th;
	get_th_covariance_nuisance(&covariance_th, ctx->names, &nobs, &param, &ctx->nuisparam, ctx->nuiscorr);

	for(ie = 0; ie < nobs; ie++)
		fprintf(out, "\t%-25s =\t%.5e +/- %.5e\t(exp: %.5e +/- %.5e)\n",
			ctx->names[ie], predictions[ie], sqrt(covariance_th[ie][ie]), ctx->central_exp[ie], sqrt(ctx->covariance_exp[ie][ie]));

	double **covariance_tot;
	get_covtot(&covariance_th, &ctx->covariance_exp, &covariance_tot, nobs);

	double **inv_cov_tot;
	if(get_invcovtot(&covariance_tot, &inv_cov_tot, nobs) < 0)
	{
		fprintf(out, "Invalid total covariance matrix\n");
		return 1;
	}

	double chi2 = get_chi2(inv_cov_tot, predictions, ctx->central_exp, nobs);
	fprintf(out, "\nchi2 = %.2f\n", chi2);
	fprintf(out, "n_obs=%d\n\n", nobs);

	write_covariance("slha_covariance.out", covariance_tot, ctx->names, nobs);
	fprintf(out, "Total covariance matrix written in slha_covariance.out\n");

	free(predictions);

	for(ie = 0; ie < nobs; ie++) free(covariance_th[ie]);
	free(covariance_th);

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

static int process_worker_request(const char *request, worker_context *ctx)
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
	rc = evaluate_slha_file(stdout, file_name, ctx);
	printf("===END\t%s\t%d===\n", request_id, rc);
	fflush(stdout);
	return rc;
}

int main(int argc, char **argv)
{
	worker_context ctx;
	if(!init_worker_context(&ctx))
	{
		fprintf(stderr, "Failed to initialize SuperIso worker context\n");
		return 1;
	}

	if(argc == 2)
	{
		int rc = evaluate_slha_file(stdout, argv[1], &ctx);
		free_worker_context(&ctx);
		return rc;
	}

	if(argc != 1)
	{
		free_worker_context(&ctx);
		fprintf(stderr, "Usage: %s [slha_file]\n", argv[0]);
		fprintf(stderr, "Or run without arguments and stream '<id>\\t<slha_file>' requests on stdin.\n");
		return 1;
	}

	char line[WORKER_LINE_MAX];
	while(fgets(line, sizeof(line), stdin) != NULL)
	{
		if(line[0] == '\n' || line[0] == '\r' || line[0] == '\0') continue;
		process_worker_request(line, &ctx);
	}

	free_worker_context(&ctx);
	return 0;
}