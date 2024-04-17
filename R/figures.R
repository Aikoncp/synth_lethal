library(tidyverse)
library(scales)


#Expression vs essentiality
gene_A <- "EGFR"
gene_B <- "STMN1"

x <- gene_effect |> select("EGFR") |> pull()
y <- expression |> select("STMN1") |> pull()


plot_tibble <- tibble(Essentiality = x, Expression = y)

ggplot(data = plot_tibble, aes(x=Essentiality, y=Expression)) + geom_point(color = "#73a9c9") + labs(title="Expression of STMN1 vs Essentiality of EGFR",
                                                                                   x ="Essentiality of EGFR", y = "Expression of STMN1")





########comparison of p-values


ggplot(
  data = SL_pairs |> dplyr::filter(gene1 == "ERBB2"),
  aes(
    x = q_value,
    y = depletion_q_value,
    colour = "cadetblue"
  )
) + geom_point() + labs(title = "Comparison of p-values of EGFR",
                        x =
                          "Wilcoxon p-value", y = "Hypergeometric Depletion p-value")







pearson <-
  cor(
    SL_pairs |> dplyr::filter(gene1 == "EGFR") |> select(p_value),
    SL_pairs |> dplyr::filter(gene1 == "EGFR") |> select(depletion_p_value),
    method = "pearson",
    use = "complete.obs"
  )
ggplot(
  data = SL_pairs |> dplyr::filter(gene1 == "EGFR"),
  aes(x = p_value + 1e-30,
      y = depletion_p_value + 1e-25)
) + scale_x_log10(
  breaks = trans_breaks("log10", function(x)
    10 ^ x),
  labels = trans_format("log10", math_format(10 ^
                                               .x))
) +
  scale_y_log10(
    breaks = trans_breaks("log10", function(x)
      10 ^ x),
    labels = trans_format("log10", math_format(10 ^ .x))
  ) + geom_point(size = 1, color = "#73a9c9")  + labs(
    title = str_c(
      "Plot of p-values of the Wilcoxon and Hypergeometric Depletion Test of EGFR, Pearson Correlation = ",
      round(pearson, 3)
    ),
    x = "Wilcoxon p-value",
    y = "Hypergeometric Depletion p-value"
  )


pearson <-
  cor(
    SL_pairs |> dplyr::filter(gene1 == "EGFR") |> select(p_value),
    SL_pairs |> dplyr::filter(gene1 == "EGFR") |> select(survival_coef),
    method = "pearson",
    use = "complete.obs"
  )
ggplot(data = SL_pairs |> dplyr::filter(gene1 == "EGFR"),
       aes(x = p_value + 1e-20, y = survival_coef)) + geom_point(size = 1, color = "#73a9c9") + scale_x_log10(
         breaks = trans_breaks("log10", function(x)
           10 ^ x),
         labels = trans_format("log10", math_format(10 ^
                                                      .x))
       ) + labs(
         title = str_c(
           "Plot of p-values of the Wilcoxon Test and coefficients of the Cox model of EGFR, Pearson Correlation = ",
           round(pearson, 3)
         ),
         x =
           "Wilcoxon p-value",
         y = "Cox model coefficient"
       )

pearson <-
  cor(
    SL_pairs |> dplyr::filter(gene1 == "EGFR") |> select(survival_coef),
    SL_pairs |> dplyr::filter(gene1 == "EGFR") |> select(depletion_p_value),
    method = "pearson",
    use = "complete.obs"
  )
ggplot(data = SL_pairs |> dplyr::filter(gene1 == "EGFR"),
       aes(x = depletion_p_value + 1e-15, y = survival_coef)) + geom_point(size = 1, color = "#73a9c9") + scale_x_log10(
         breaks = trans_breaks("log10", function(x)
           10 ^ x),
         labels = trans_format("log10", math_format(10 ^
                                                      .x))
       ) + labs(
         title = str_c(
           "Plot of p-values of the Hypergeometric Depletion Test and coefficients of the Cox model of EGFR, Pearson Correlation = ",
           round(pearson, 3)
         ),
         x =
           "Hypergeometric Depletion p-value",
         y = "Cox model coefficient"
       )

####Histogram of p-value and q-value

ggplot(data = SL_pairs, aes(x=p_value)) + geom_histogram(binwidth = 0.005, colour =  "black", fill =  "#73a9c9") + labs(title="Histogram of Wilcoxon p-values",
                                                                                                                       x ="Wilcoxon p-value", y = "Frequency")

ggplot(data = SL_pairs, aes(x=q_value)) + geom_histogram(binwidth = 0.005, colour =  "black", fill =  "#73a9c9") + labs(title="Histogram of Wilcoxon q-values",
                                                                                                                        x ="Wilcoxon q-value", y = "Frequency")