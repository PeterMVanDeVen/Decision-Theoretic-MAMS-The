#' R6 Trial Class
#'
#' @description
#' This is the class used to define a clinical trial object.
#'
#' @export
Trial <- R6::R6Class(

    # Class Trial ----
    classname = 'Trial',

    # * Public Members ----
    public = list(

        #' @description
        #' Create a new trial object.
        initialize = function() {

        },

        solve = function() {

            n_stage_per_arm <- 20

            private$data %>%
                mutate(id_stage = ceiling(id_patient / n_per_stage_per_arm)) %>%
                group_by(id_trial, id_arm, id_stage) %>%
                summarise(y = sum(y), n = n(), .groups = 'drop') %>%
                tidyr::pivot_wider(names_from = 'id_arm', values_from = c('y', 'n')) %>%
                tidyr::nest(Y = starts_with('y'), N = starts_with('n')) %>%
                mutate(
                    opt_dec =
                        purrr::map2_int(
                            Y, N,
                            function(x, y) {
                                optimalDecision(unlist(x), unlist(y), prior = c(1, 1), decisions = 1:4, loss = binaryLoss)
                            }
                        )
                )


        }


    ),

    # * Private Members ----
    private = list(

        # Private attributes
        data = NULL,

        # Private methods

        create = function(resp_rate, n_trials, n_max) {
            private$data <-
                tibble(
                    y = runif(n = length(resp_rate) * n_trials * n_max),
                    id_trial = rep(1:n_trials, times = n_max * length(resp_rate)),
                    id_arm = rep(1:length(resp_rate), times = n_trials * n_max)
                ) %>%
                left_join(
                    enframe(resp_rate, name = 'id_arm', value = 'resp_rate'),
                    by = 'id_arm'
                ) %>%
                mutate(y = as.integer(y < resp_rate)) %>%
                group_by(id_arm, id_trial) %>%
                mutate(id_patient = row_number()) %>%
                ungroup() %>%
                arrange(id_trial, id_patient, id_arm) %>%
                select(id_trial, id_arm, id_patient, y)
        },

        binaryLoss = function(decision, theta, delta = 0) {

            # decision is an integer vector
            # theta is a tibble with n rows and k + 1 columns (k = nr. experimental arms)

            # Decisions:
            # 1. no experiental arm is better than control
            # 2. arm 2 is better than control
            # 3. arm 3 is better than control
            # 4. both arms 2 and 3 are better than control

            # Loss:
            # 0, if the decision is the correct one
            # 1, if the decision is NOT the correct one

            # Create a tibble with names theta_1, theta_2, etc.
            tbl_theta <-
                as_tibble(theta) %>%
                magrittr::set_colnames(paste0('theta_', 1:ncol(theta)))

            tbl_theta %>%
                mutate(
                    # Subtract the first col (theta_1) to all others
                    across(-1, ~.x - theta_1),
                    # Check whether the difference is > delta
                    across(-1, ~.x > delta),
                    # Proposed decision
                    dec = as.integer(decision),
                    # Optimal decision
                    opt_dec = case_when(
                        theta_2 & theta_3 ~ 4L,
                        theta_3 ~ 3L,
                        theta_2 ~ 2L,
                        TRUE ~ 1L
                    ),
                    # Zero-one loss function
                    loss = as.integer(dec != opt_dec)
                ) %>%
                pull(loss)

        },

        optimalDecision = function(y, n, prior, decisions = 1:4, loss, ...) {

            # Number of arms
            K <- length(y)

            # Sample size for posterior draws per arm
            G <- 1e4

            # Prior parameters
            a <- prior[1]
            b <- prior[2]

            # Matrix of posterior draws (G x K)
            theta_post <-
                rbeta(
                    n = G * K,
                    shape1 = a + y,
                    shape2 = b + n - y
                ) %>%
                matrix(ncol = K, byrow = TRUE)

            # Compute the loss for each sample, for each decision.
            # Then the mean is computed, and the decision with the
            # lowest mean loss is the optimal.
            loss_val <-
                lapply(
                    decisions,
                    loss,
                    theta = theta_post,
                    ... = ...
                ) %>%
                lapply(mean) %>%
                unlist()

            # Optimal decision: argmin(E[Loss])
            opt_dec <- decisions[which.min(loss_val)]

            return(opt_dec)

        }

    )

)
