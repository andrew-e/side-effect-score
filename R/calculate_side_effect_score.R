#' calculate_pki_nmues: calculating pKi values from a set of Ki(nM) per drug/receptor.  Using the following formula
#'   pKi = -log10(mean(KiM)) - 4, where KiM = Ki(nM) * 1e-8
#' @param binding_affinities: binding affinity data with columns drug, receptor, ki_nm
#' @returns: data frame with pki_mean and pki_se
convert_ki_to_mean_pki <- function(binding_affinities) {
  DEFAULT_KI_MEAN <- mean(binding_affinities$ki_nm)
  DEFAULT_KI_SE <- sd(binding_affinities$ki_nm) / sqrt(nrow(binding_affinities))

  pkis <- dplyr::select(binding_affinities, receptor, drug, ki_nm) |>
    dplyr::mutate(ki_nm = -log10(ki_nm*1e-8) - 4) |>
    dplyr::group_by(receptor, drug) |>
    dplyr::summarise_each(list(
      pki_mean = function(x) { return(mean(x, trim = 0.2)) },
      pki_se = function(x) { sd(x)/sqrt(length(x)) })
    )

  #set default values here...
  pkis$pki_mean[is.na(pkis$pki_mean)] <- DEFAULT_KI_MEAN
  pkis$pki_se[is.na(pkis$pki_se)] <- DEFAULT_KI_SE
  return(pkis)
}

#' mr_wald_ratio: formula for calculating the wald ratio, using an adapted formula as found in the TwoSampleMR package
#'  BIV = BZY / BZX, where BZY is the PheWAS result, and BZX is the QTL result
#' @param phewas: data frame of phewas results
#' @returns data frame with addition wald ratio results
calculate_wald_ratio <- function(phewas) {
  phewas$wald_ratio_beta <- phewas$phewas_beta / phewas$qtl_beta
  phewas$wald_ratio_se <- phewas$phewas_se / abs(phewas$qtl_beta)
  phewas$wald_ratio_pval <- stats::pnorm(abs(phewas$wald_ratio_beta) / phewas$wald_ratio_se, lower.tail=FALSE) * 2

  return(phewas)
}

#' calulcate_bootstrap_values: calculating a bootstrapped result based on two paris of mean/standard error
#' @param pki_values
#' @param mr_results
#' @returns data frame of bootstrap_mean_abs and bootstrap_se attached to the existing data frame
#'
calculate_bootstrap_values <- function(pki_values, mr_results) {
  #Only choose a single wald ratio result per receptor/side effect pair.  Use the largest effect size
  results <- merge(mr_results, pki_values, by="receptor") |>
    dplyr::arrange(side_effect, receptor, dplyr::desc(abs(wald_ratio_beta)))
  results <- results[!duplicated(results[,c("receptor", "side_effect")])]

  #Perform bootstrapping for every receptor binding affinity / side effect mr result pair
  bootstrap_data <- results[, c("wald_ratio_beta", "wald_ratio_se", "pki_mean", "pki_se")]
  bootstrap_results <- apply(bootstrap_data, 1, function(result) {
    first_bootstrap <- rnorm(1000, mean=result[['wald_ratio_beta']], sd=result[['wald_ratio_se']])
    second_bootstrap <- rnorm(1000, mean=result[['pki_mean']], sd=result[['pki_se']])
    joined_bootstrap <- first_bootstrap * second_bootstrap
    return(list(bootstrap_mean_abs=abs(mean(joined_bootstrap)), bootstrap_se=sd(joined_bootstrap) / sqrt(2)))
  }) |> dplyr::bind_rows()

  results <- dplyr::bind_cols(results, bootstrap_results)

  return(results)
}


binding_affinities <- data.table::fread("data/ki_summary_data.tsv") |>
  dplyr::filter(ki_nm < 10000)

#phewas_beta, phewas_se, qtl_beta
phewas_results <- data.table::fread("data/phewas_and_qtl_data.tsv")

pki_values <- convert_ki_to_mean_pki(binding_affinities)
mr_results <- calculate_wald_ratio(phewas_results)
results <- calculate_bootstrap_values(pki_values, mr_results)

#Now you can aggregate the results by receptor, or by drug
score_by_receptor <- dplyr::group_by(results, receptor) |>
  dplyr::summarise(score=sum(bootstrap_mean_abs))

score_by_side_effect <- dplyr::group_by(results, side_effect) |>
  dplyr::summarise(score=sum(bootstrap_mean_abs))