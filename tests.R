evalnodrop <- function(y, n1, n2, prior, delta = 0, gamma){
    if(delta < 0){
        cat("ERROR: delta must be a non-negative number \n")
        return(NULL)
    }
    J <- length(y)
    if(gamma>0) {G <- ceiling(2.5*(1/gamma))} else {G<-10000} # number of draws per arm
    a <- prior[1]
    b <- prior[2]
    benefit <- NA
    theta <- matrix(rbeta(G*J, shape1 = a + y, shape2 = b + n1 - y), ncol = J, byrow = T)
    Tt <- apply(theta, 1, Ttheta, delta)
    Ty <- Tdata(y = y, n = n1, prior, delta = delta)
    Ynew <- matrix(rbinom(n = G*J, size = n2, prob = t(theta)), ncol = J, byrow = T)
    Tynew <- apply( t(t(Ynew) + y), 1, Tdata, n1+n2, prior, delta)
    b0 <- sum(Ty == Tt)/G
    benefit[1] <- sum(Tynew == Tt)/G
    return(list(ben.stop = b0, ben.cont = benefit))}


###########################


n_per_stage_per_arm <- 20

tbx <-
    private$data %>%
    mutate(id_stage = ceiling(id_patient / n_per_stage_per_arm)) %>%
    filter(id_stage == 1) %>%
    group_by(id_trial, id_arm, id_stage) %>%
    summarise(y = sum(y), n = n(), .groups = 'drop') %>%
    tidyr::pivot_wider(names_from = 'id_arm', values_from = c('y', 'n')) %>%
    tidyr::nest(Y = starts_with('y'), N = starts_with('n')) %>%
    mutate(
        opt_dec =
            purrr::map2(
                Y, N,
                function(x, y) {
                    optimalDecision(unlist(x), unlist(y), prior = c(1, 1), decisions = 1:4, loss = binaryLoss)
                }
            )
    )

# Sample new data for an extra stage, using the posterior draws for response
# rate (theta).
tbl_y_new <-
    tbx %>%
    mutate(N_new = N) %>% # TODO da rimuovere dopo
    tidyr::unnest(N_new) %>%
    tidyr::unnest(opt_dec) %>%
    tidyr::unnest(theta) %>%
    select(id_trial, id_stage, starts_with('n_'), starts_with('theta'))

for (i in 1:K) {
    tbl_y_new %<>%
        mutate(
            !! paste0('y_new_', i) :=
                rbinom(
                    n = nrow(.),
                    size = get(paste0('n_', i)),
                    prob = get(paste0('theta_', i))
                )
        )
}

# Add new data (samples)
tbx %<>%
    left_join(
        tbl_y_new %>%
            tidyr::nest(Y_new = starts_with('y_new')) %>%
            select(id_trial, id_stage, Y_new),
        # tbl_y_new %>%
        #     select(id_trial, id_stage, starts_with('y_new')) %>%
        #     group_by(id_trial, id_stage) %>%
        #     tidyr::nest(Y_new = starts_with('y_new')),
        by = c('id_trial', 'id_stage')
    )

# New step evaluation
tbx %<>%
    mutate(
        opt_dec_new =
            purrr::map2(
                Y_new, N,
                function(x, y) {
                    optimalDecision(unlist(x), unlist(y), prior = c(1, 1), decisions = 1:4, loss = binaryLoss)
                }
            )
    )

tbx %<>%
    select(id_trial, id_stage, opt_dec, opt_dec_new) %>%
    tidyr::unnest(opt_dec) %>%
    select(id_trial, id_stage, opt_dec_0 = opt_dec, opt_loss_0 = opt_loss, opt_dec_new) %>%
    tidyr::unnest(opt_dec_new) %>%
    select(-theta) %>%
    group_by(id_trial, id_stage) %>%
    summarise(
        opt_loss_0 = min(opt_loss_0),
        opt_dec_0 = min(opt_dec_0),
        across(opt_loss, mean),
        .groups = 'drop'
    )

##########################################
# new test, easier but requires a lot of memory
##########################################

K <- 3 # nr. arms
n_per_stage_per_arm <- 20

tbx <-
    private$data %>%
    mutate(id_stage = ceiling(id_patient / n_per_stage_per_arm)) %>%
    group_by(id_trial, id_arm, id_stage) %>%
    summarise(y = sum(y), n = n(), .groups = 'drop') %>%
    tidyr::pivot_wider(names_from = 'id_arm', values_from = c('y', 'n'))

G <- 1e2
tbl_theta <- tibble(
    id_trial = rep(tbx$id_trial, each = G),
    id_stage = rep(tbx$id_stage, each = G)
) %>%
    left_join(
        tbx,
        by = c('id_trial', 'id_stage')
    )
for (i in 1:K) {
    tbl_theta %<>%
        mutate(
            !! paste0('theta_', i) :=
                rbeta(
                    n = nrow(.),
                    shape1 = 1 + get(paste0('y_', i)),
                    shape2 = 1 + get(paste0('n_', i)) - get(paste0('y_', i))
                )
        )
}

Q <- 1e2
tbl_y_new <-
    tibble(
        id_trial = rep(tbl_theta$id_trial, each = Q),
        id_stage = rep(tbl_theta$id_stage, each = Q)
    ) %>%
    left_join(
        tbl_theta,
        by = c('id_trial', 'id_stage')
    )
for (i in 1:K) {
    tbl_y_new %<>%
        mutate(
            !! paste0('y_new_', i) :=
                rbinom(
                    n = nrow(.),
                    size = get(paste0('n_', i)),
                    prob = get(paste0('theta_', i))
                )
        )
}


