library(tidyverse)
library(AnnotationDbi)
library(org.Hs.eg.db)

data_folder <- '/cellfile/projects/rna-kinetics/rna_decay_cafe/input_data'
output_folder <- '/cellfile/projects/rna-kinetics/rna_decay_cafe/output_figures'

rna_decay_cafe <- read_csv(file.path(data_folder, 'AvgKdegs_genes_v1.csv'))

sd_kdeg_by_gene <- rna_decay_cafe |>
  dplyr::group_by(feature_ID) |>
  dplyr::summarise(
    sd_avg_donorm_log_kdeg = stats::sd(avg_donorm_log_kdeg, na.rm = TRUE),
    .groups = "drop"
  ) |>
  drop_na()

symbol_to_ensg <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(sd_kdeg_by_gene$feature_ID),
  keytype = "SYMBOL",
  columns = c("ENSEMBL")
) |>
  as_tibble() |>
  filter(!is.na(ENSEMBL)) |>
  distinct(SYMBOL, ENSEMBL)

sd_kdeg_by_gene <- sd_kdeg_by_gene |>
  left_join(
    symbol_to_ensg,
    by = c("feature_ID" = "SYMBOL")
  ) |>
  drop_na()

test_results <- read_csv(file.path(data_folder, 'test_results.csv'))


for (kinetic_rate in c('elongation_speed', 'splicing_speed')) {
  for (treatment_group in c('replicative_senescent', 'ICM')) {

    print(paste("Kinetic rate:", kinetic_rate, "Treatment group:", treatment_group))
    test_results_subset <- test_results |>
      dplyr::filter(
        parameter_type == kinetic_rate,
        group_1 == treatment_group,
        group_2 == "proliferating") |>
      left_join(sd_kdeg_by_gene, join_by(gene_name == ENSEMBL)) |>
      drop_na(sd_avg_donorm_log_kdeg) |>
      mutate(abs_lfc_regularized = abs(lfc_regularized)) |>
      dplyr::filter(p_value < 0.05)

    df_cor <- test_results_subset |>
      dplyr::select(sd_avg_donorm_log_kdeg, abs_lfc_regularized) |>
      tidyr::drop_na()

    lm_fit <- lm(abs_lfc_regularized ~ sd_avg_donorm_log_kdeg, data = df_cor)

    lm_summary <- summary(lm_fit)
    r2 <- lm_summary$r.squared
    p_lm <- lm_summary$coefficients[2, 4]  # slope p-value
    print(lm_summary)

    scatter_plot <- ggplot(
      df_cor,
      aes(
        x = sd_avg_donorm_log_kdeg,
        y = abs_lfc_regularized
      )
    ) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = "lm", se = TRUE, color = "firebrick") +
      theme_minimal() +
      labs(
        x = "Sd. of log(k_deg)",
        y = paste("|log2FC| kinetic rate", kinetic_rate),
        title = paste0("LFC ", kinetic_rate, "(", treatment_group, " vs. proliferating) vs. sd. of (log(k_deg))")
      ) +
      annotate(
        "text",
        x = Inf, y = Inf,
        hjust = 1.1, vjust = 1.5,
        label = sprintf(
          "R² = %.2f\np = %.2g",
          r2, p_lm
        ),
        size = 4
      )
    print(scatter_plot)

    ggsave(file.path(output_folder, paste0(kinetic_rate, "_", treatment_group, '.png')),
           bg = 'white')

  }
}


