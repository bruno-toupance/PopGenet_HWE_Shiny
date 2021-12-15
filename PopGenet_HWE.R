#==============================================================================
#    PopGenet_HWE.R : Population Genetics HWE
#    Copyright (C) 2021  Bruno Toupance <bruno.toupance@mnhn.fr>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by	
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#==============================================================================



library(adegenet)
library(pegas)
library(HardyWeinberg)


#==============================================================================
# check_integer
#==============================================================================
check_integer <- function(
	nb = 1, 
	label = "nb", 
	check_list = list(check_flag = TRUE, check_msg = "")) 
{
	if (is.numeric(nb)) {
		if (nb == round(nb)) {
			if (nb < 0) {
				check_list$check_flag <- FALSE
				check_list$check_msg <- sprintf("%s\nFAIL: [%s] out of bounds", check_list$check_msg, label)
			}
		} else {
			check_list$check_flag <- FALSE
			check_list$check_msg <- sprintf("%s\nFAIL: [%s] not integer", check_list$check_msg, label)
		}
	} else {
		check_list$check_flag <- FALSE
		check_list$check_msg <- sprintf("%s\nFAIL: [%s] not numeric", check_list$check_msg, label)
	}
#------------------------------------------------------------------------------
	return(check_list)
}




#==============================================================================
# check_parameters
#==============================================================================
check_parameters <- function(
	nb_AA = 10, nb_Aa = 10, nb_aa = 10,
	pegas_B = 10) 
{
	check_list <- list(check_flag = TRUE, check_msg = "")
#------------------------------------------------------------------------------
	check_list <- check_integer(nb_AA, label = "nb_AA", check_list)
	check_list <- check_integer(nb_Aa, label = "nb_Aa", check_list)
	check_list <- check_integer(nb_aa, label = "nb_aa", check_list)
	check_list <- check_integer(pegas_B, label = "pegas_B", check_list)
#------------------------------------------------------------------------------
	if (check_list$check_flag) {
		nb_ind <- nb_AA + nb_Aa + nb_aa
		if (nb_ind == 0) {
			check_list$check_flag <- FALSE
			check_list$check_msg <- sprintf("%s\nFAIL: no individual", check_list$check_msg)
		}
	}
#------------------------------------------------------------------------------
	return(
		list(
			msg = check_list$check_msg, 
			flag = check_list$check_flag, 
			nb_AA = nb_AA, 
			nb_Aa = nb_Aa, 
			nb_aa = nb_aa,
			pegas_B = pegas_B
		)
	)
}





#==============================================================================
# compute_pop
#  Compute population statistics from genotypic structure
#==============================================================================
compute_pop <- function(
	nb_AA = 363, 
	nb_Aa = 634, 
	nb_aa = 282,
	pegas_perm_flag = FALSE,
	pegas_B = 1000,
	HardyWeinberg_perm_flag = FALSE,
	HardyWeinberg_B = 1000)
{
	param_checking <- check_parameters(nb_AA, nb_Aa, nb_aa, pegas_B)
	# print(param_checking)

#------------------------------------------------------------------------------
	if (param_checking$flag) {
	
# Build the genotype and allele names
		genotype_name <- c("AA", "Aa", "aa")
		allele_name <- substr(genotype_name[c(1, 3)], start = 1, stop = 1)
		
# Compute the genotype count anf frequency
		genotype_count <- c(nb_AA, nb_Aa, nb_aa)
		names(genotype_count) <- genotype_name
		nb_ind <- sum(genotype_count)
		genotype_freq <- genotype_count / nb_ind
		
# Compute the allele count and frequency
		allele_count <- as.integer(2 * genotype_count[c(1, 3)] + genotype_count[2])
		names(allele_count) <- allele_name
		nb_copy <- sum(allele_count)
		allele_freq <- allele_count / nb_copy
		
# Classical chi-squared computation
		p <- allele_freq[1]
		hwe_exp_freq <- c(p^2, 2 * p * (1 - p), (1 - p)^2)
		hwe_exp_count <- hwe_exp_freq * nb_ind
		hwe_obs_count <- genotype_count
		hwe_chisq_resid <- (hwe_obs_count - hwe_exp_count) / sqrt(hwe_exp_count)
		hwe_chisq_stat <- sum(hwe_chisq_resid^2)
		hwe_chisq_df <- 1
		hwe_chisq_pvalue <- 1 - pchisq(hwe_chisq_stat, df = hwe_chisq_df)
		

# Computation with "pegas" package
		ind_vec <- rep(genotype_name, times = genotype_count)
		ind_lab <- paste("ind_", 1:nb_ind, sep = "")
		ind_df <- data.frame(loc_1 = ind_vec)
		ind_gi <- df2genind(ind_df, sep = "", ploidy = 2)

		if (pegas_perm_flag) {
			hwe_exact_test = hw.test(ind_gi, B = pegas_B)
		} else {
			hwe_exact_test = NULL
		}

		pegas_data = list(
			pegas_B = pegas_B,
			ind_vec = ind_vec,
			ind_lab = ind_lab,
			ind_gi = ind_gi,
			hwe_exact_test = hwe_exact_test
		)
		
		
# Computation with "HardyWeinberg" package
#   see Graffelman J 2021 Exploring Diallelic Genetic Markers: The HardyWeinberg Package
		HW_test_1 <- HWChisq(genotype_count, cc = 0.0, verbose = FALSE) 
		HW_test_2 <- HWChisq(genotype_count, cc = 0.5, verbose = FALSE) 
		HW_test_3 <- HWLratio(genotype_count, verbose = FALSE) 
		HW_test_4 <- HWExact(genotype_count, pvaluetype = "selome", verbose = FALSE) 
		HW_test_5 <- HWExact(genotype_count, pvaluetype = "dost", verbose = FALSE) 
		HW_test_6 <- HWExact(genotype_count, pvaluetype = "midp", verbose = FALSE) 
		if (HardyWeinberg_perm_flag) {
			HW_test_7 <- HWPerm(genotype_count, nperm = HardyWeinberg_B, verbose = FALSE)
		} else {
			HW_test_7 <- data.frame(stat = NA, pval = NA)
		}

		HardyWeinberg_data <- data.frame(
			Test = c(
				"Chi-square test:", 
				"Chi-square test with continuity correction:",
				"Likelihood-ratio test:", 
				"Exact test with selome p-value:",
				"Exact test with dost p-value:",
				"Exact test with mid p-value:", 
				"Permutation test:"
			),
			Statistic = c(
				HW_test_1$chisq, 
				HW_test_2$chisq, 
				HW_test_3$G2, 
				NA, 
				NA, 
				NA, 
				HW_test_7$stat
			), 
			p_value = c(
				HW_test_1$pval, 
				HW_test_2$pval, 
				HW_test_3$pval, 
				HW_test_4$pval, 
				HW_test_5$pval, 
				HW_test_6$pval,
				HW_test_7$pval
			)
		)
		# colnames(HardyWeinberg_data) <- c("Test", "Statistic", "p-value")
		
		
		

# Compute the observed allele count and frequency table
		allele_obs_df <- data.frame(
			allele = allele_name, 
			count = allele_count,
			frequency = allele_freq
		)
		
# Compute the observed genotype count and frequency table
		genotype_obs_df <- data.frame(
			genotype = genotype_name, 
			count = genotype_count,
			frequency = genotype_freq
		)
		
# Compute the chi-square table
		hwe_df <- data.frame(
			genotype = genotype_name, 
			observed = hwe_obs_count,
			expected = hwe_exp_count,
			residuals = hwe_chisq_resid
		)
		
# Compute the diversity table
		diversity_df <- data.frame(
			n = nb_ind, 
			D = (hwe_obs_count[2] - hwe_exp_count[2]) / 2, 
			Ho = genotype_freq[2],
			He = hwe_exp_freq[2],
			F = 1 - genotype_freq[2] / hwe_exp_freq[2]
		)
		
		pop_data <- list(
			genotype_name = genotype_name,
			allele_name = allele_name,
			
			nb_ind = nb_ind,
			genotype_count = genotype_count,
			genotype_freq = genotype_freq,
			
			nb_copy = nb_copy,
			allele_count = allele_count,
			allele_freq = allele_freq,
			
			hwe_exp_count = hwe_exp_count,
			hwe_obs_count = hwe_obs_count,
			hwe_chisq_resid = hwe_chisq_resid,
			hwe_chisq_stat = hwe_chisq_stat,
			hwe_chisq_df = hwe_chisq_df,
			hwe_chisq_pvalue = hwe_chisq_pvalue,
			
			pegas_perm_flag = pegas_perm_flag,
			pegas_data = pegas_data,
			HardyWeinberg_perm_flag = HardyWeinberg_perm_flag,
			HardyWeinberg_data = HardyWeinberg_data,

			genotype_obs_df = genotype_obs_df,
			allele_obs_df = allele_obs_df,
			diversity_df = diversity_df,

			hwe_df = hwe_df,

			param = param_checking
		)
#------------------------------------------------------------------------------
	} else {
		pop_data <- list(param = param_checking)
	}
#------------------------------------------------------------------------------

	return(pop_data)
}



#==============================================================================
# plot_error
#==============================================================================
plot_error <- function(msg)
{
	plot(c(0, 1), c(0, 1), type = "n", xlab = "", ylab = "", main = "", 
		xaxt = "n", yaxt = "n", bty = "n")
	err_msg <- sprintf("ERROR: check parameter values%s", msg)
	text(0.5, 0.5, err_msg, col = "red", adj = 0.5)
}



#==============================================================================
# plot_triangle
#==============================================================================
plot_triangle <- function(pop_data) 
{
	if (pop_data$param$flag) {

		## Retreive allele frequency, heterozygous frequency and sample size.
		p <- pop_data$allele_freq[1]
		H <- pop_data$genotype_freq[2]
		n <- pop_data$nb_ind

		## Save parameters.
		par_bak <- par(no.readonly = TRUE)

		## Adjust margins.
		par(mar = c(0, 0, 0, 0))

		## Prepare graphic area.
		inc_x <- 0.05
		inc_y <- 0.05
		box_x <- c(0 - inc_x, 1 + inc_x)
		box_y <- c(0 - inc_y, 1 + inc_y)
		plot(box_x, box_y, type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1)

		## Invalid area
		p_lwr <- min(c(sqrt(5 / n), 0.5))
		val_x <- c(0, p_lwr, p_lwr, 0)
		val_y <- c(0, 2 * p_lwr, 0, 0)
		polygon(val_x, val_y, col = "lightgray", lty = 3)

		p_upr <- 1 - p_lwr
		val_x <- c(1, p_upr, p_upr, 1)
		val_y <- c(0, 2 - 2 * p_upr, 0, 0)
		polygon(val_x, val_y, col = "lightgray", lty = 3)


		## Draw triangle.
		val_x <- c(0, 1/2, 1, 0)
		val_y <- c(0, 1, 0, 0)
		lines(val_x, val_y)

		## Draw genotype labels.
		text(0, 0, pop_data$genotype_name[2], pos = 2)
		text(1, 0, pop_data$genotype_name[1], pos = 4)
		text(1/2, 1, pop_data$genotype_name[3], pos = 3)

		## Draw HWE parabola.
		val_x <- seq(from = 0, to = 1, length = 200)
		val_y <- 2 * val_x * (1 - val_x)
		lines(val_x, val_y, lty = 2)

		## Draw genotypic structure point.
		val_x <- p
		val_y <- H
		alpha <- 0.05
		if (pop_data$hwe_chisq_pvalue <= alpha) {
			points(val_x, val_y, col = "red", pch = 19)
		} else {
			points(val_x, val_y, col = "blue", pch = 19)
		}
		
		## Chi-square limits
		alpha <- 0.05
		chisq_lim <- qchisq(1 - alpha, df = 1)

		F_lwr <- sqrt(chisq_lim / n)
		val_x <- seq(from = 0, to = 1, length = 200)
		val_y <- 2 * (val_x) * (1 - val_x) * (1 - F_lwr)
		lines(val_x, val_y, col = "gray")
		
		F_upr <- -sqrt(chisq_lim / n)
		val_x <- seq(from = -F_upr / (1 - F_upr), to = 1 / (1 - F_upr), length = 200)
		val_y <- 2 * (val_x) * (1 - val_x) * (1 - F_upr)
		lines(val_x, val_y, col = "gray")
		
		
		

		## Restore parameters.
		par(par_bak)
#------------------------------------------------------------------------------
	} else {
		plot_error(pop_data$param$msg)
	}
}



#==============================================================================
# non_interactive
#==============================================================================
non_interactive <- function()
{
	pop_data <- compute_pop() 
	print(pop_data)
	print(pop_data$obs_df)
	print(pop_data$hwe_df)
	print(pop_data$pegas_data)
	print(str(pop_data$pegas_data))
	
	return(pop_data)
}




