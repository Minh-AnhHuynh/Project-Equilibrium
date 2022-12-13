
# My helper function
make_mean_comparison <-
  function(data,
           quantvar,
           groupvar,
           levene = TRUE,
           do_kruskal = FALSE,
           agostino = FALSE,
           bartlett = FALSE,
           return_table = FALSE,
           return_tukey = FALSE) {
    librarian::shelf(rstatix, glue, moments, car)
    # moments::agostino()
    # car::levene()

    quantvar_sym <- sym(quantvar)

    if (str_detect(quantvar, " ") == TRUE) {
      quantvar <- glue("`{quantvar}`")
    }
    if (str_detect(groupvar, " ") == TRUE) {
      groupvar <- glue("`{groupvar}`")
    }
    formula_expression <- formula(glue("{quantvar} ~ {groupvar}"))

    df_name <- deparse(substitute(data))


    message("Test normality with Shapiro-Wilk test")

    # Test normality for each group ----
    shapiro <- data %>%
      group_by(!!sym(groupvar)) %>%
      group_nest() %>%
      mutate(shapiro = map(.data$data, ~ shapiro_test(.x, !!quantvar_sym))) %>%
      select(-data) %>%
      unnest(cols = shapiro) %>%
      print()


    if (all(shapiro$p <= 0.05)) { # Every group needs to satisfy the condition
      message(
        glue(
          "Null hypothesis : Assume that the data is normally distributed.
               p<= 0.05, so H0 is rejected. {quantvar} is **NOT normally distributed**"
        )
      )
    } else {
      message(
        glue(
          "Null hypothesis : Assume that the data is normally distributed.
               p>= 0.05, so H0 is not rejected. {quantvar} is **normally distributed**."
        )
      )
    }



    if (agostino == TRUE) {
      message("Test normality with D'Agostino")
      ag <-
        moments::agostino.test(pull(data, {{ quantvar }})) %>% print()
      if (ag[["p.value"]] <= 0.05) {
        message("H0 is rejected. Data doesn't have a normal distribution")
      } else {
        message("H0 is not rejected. Data has a normal distribution")
      }
    }

    # Perform Kruskal-Wallis non-parametric test ----
    if (any(shapiro$p <= 0.05) | do_kruskal == TRUE) {
      message("Kruskal-Wallis test can be performed on non parametric data")

      kruskal <- data %>%
        kruskal_test(formula_expression) %>%
        print()
      if (kruskal$p >= 0.05) {
        message(
          "P is not significant. There is no evidence of stochastic dominance between the samples."
        )
      }

      if (kruskal$p <= 0.05) {
        message(
          glue(
            "Kruskral-Wallis p <= 0.05. \n
              At least one sample stochastically dominates another sample. \n
              Post-hoc testing is allowed."
          )
        )


        # We can do wilcoxon test but dunn_test seems to be better for kruskal
        # wallis post-hoc analysis as Dunn's test takes into account the ranking
        # used by the Kruskal-Wallis test. It also makes adjustments for
        # exaequos.
        message("Wilcoxon post hoc test")
        wilcoxon <- data %>%
          wilcox_test(formula_expression, p.adjust.method = "holm") %>%
          dplyr::select(
            "Group 1" = group1,
            "Group 2" = group2,
            n1,
            n2,
            p,
            p.adj,
            p.adj.signif
          ) %>%
          print()
        # holm method for correction is by default but seems better than
        # bonferonni


        ## Dunn post-hoc testing ----
        message("Dunn post hoc test")
        dunn <- data %>%
          dunn_test(formula_expression, p.adjust.method = "holm") %>%
          print()

        pk <- kruskal$p
        dunn_tab <- dplyr::select(dunn,
                                  "Group 1" = group1,
                                  "Group 2" = group2,
                                  n1,
                                  n2,
                                  p,
                                  p.adj,
                                  p.adj.signif
        )
        stat_table <-
          ggtexttable(dunn_tab, rows = NULL, theme = ttheme("light")) %>%
          tab_add_title(glue("Kruskal-Wallis : p = {pk}
                       Dunn test"), face = "bold") %>%
          tab_add_footnote("p.adj using Holm method.",
                           size = 10,
                           face = "italic"
          ) %>%
          print()
        file_path_save <-
          glue("04_figures/ggplot/{df_name}-Kruskal-dunn-table-{quantvar}.png")
        ggsave(file_path_save)
        print(glue("Table saved in {file_path_save}"))
        if (return_table == TRUE) {
          return(
            list(
              "kruskal" = kruskal,
              "dunn" = dunn,
              "wilcoxon" = wilcoxon,
              "bartlett" = bartlett,
              "shapiro" = shapiro
            )
          )
        }
        return(stat_table)
      } else {
        message("Kruskal-Wallis p >= 0.05, post-hoc testing is not allowed")
        return(kruskal)
      }
    }


    # ANOVA ----
    if (bartlett == TRUE & all(shapiro$p >= 0.05)) {
      # Test homogeneity of variances in order to perform anova
      bartlett <-
        bartlett.test(formula_expression, data) %>% print()
      if (bartlett[["p.value"]] <= 0.05) {
        message(
          "The hypothesis that the assumption of equal variances between groups is rejected. Anova cannot be performed."
        )
      }

      if (bartlett[["p.value"]] >= 0.05) {
        message(
          "The assumption of equal variances is valid as H0 is not rejected.
              Hypothesis of homoscedasticity is confirmed.
              ANOVA can be performed."
        )
      }
    }

    ## Test homoscedasticity (levene) ====
    # TODO : Replace car::leveneTest by rstatix::levene_test()
    if (levene == TRUE & all(shapiro$p >= 0.05)) {
      levene_stat <- leveneTest(formula_expression, data) %>% print()

      if (levene_stat$`Pr(>F)`[1] <= 0.05) {
        message("H0 is rejected. Variances are NOT equal between groups.")
      } else {
        message("H0 is not rejected. Variances are equal between groups.")
      }
    }


    if (all(shapiro$p >= 0.05) & levene_stat$`Pr(>F)`[1] >= 0.05) {
      message("Perform ANOVA test")


      anova <- aov(formula_expression, data = data) %>%
        print()

      anova %>%
        report::report() %>%
        print()
      anova_data <- anova_test(data, formula_expression) %>% print()
    } else {
      message("ANOVA cannot be performed")
    }

    ## Post-hoc Tukey's test ====
    if (anova_data$p <= 0.05) {
      message("Post-hoc one-way ANOVA testing can be done.
              Performing Tukey's test")
      # The reason we do not use tukey_hsd(anova) is because add_xy_position()
      # doesn't work on it while it works on an independent tukey with data and
      # formula expression
      tukey <- tukey_hsd(data, formula_expression) %>% print()
    } else {
      message("Post-hoc analysis cannot be done.")
    }
    if (return_tukey == TRUE) {
      return(tukey)
    }
  }

