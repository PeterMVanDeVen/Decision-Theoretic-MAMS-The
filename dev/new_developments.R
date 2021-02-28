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

        }

    )

)
