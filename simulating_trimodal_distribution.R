# GR 0-1 occurrs about 200 times
# GR 1-2 occurrs about 100 times
# GR 6-8 occurrs about 100 times

ts_length <- 350
draws_a <- runif(ts_length*.25, min=.12, max=.25)
draws_b <- runif(ts_length*.5, min=.75, max=1.25)
draws_c <- runif(ts_length*.125, min=6.5, max=6.75)
draws_d <- runif(ts_length*.125, min=7.75, max=8)
draws <- c(draws_a, draws_b, draws_c, draws_d)

hist(draws, breaks = seq(0,8, by=0.25))
hist(log(draws))
mean(log(draws))

# recreating our GR recreates the trimodal distribution of log(GR)
# but why does that persist across runs?

# what if we do it across a bunch of runs
ts_length <- 350
runs <- 20
results <- data.frame(GR=double(), run=integer())
for (xx in 1:runs) {
  draws_a <- runif(ts_length*.25, min=.12, max=.25)
  draws_b <- runif(ts_length*.5, min=.75, max=1.25)
  draws_c <- runif(ts_length*.125, min=6.5, max=6.75)
  draws_d <- runif(ts_length*.125, min=7.75, max=8)
  draws <- c(draws_a, draws_b, draws_c, draws_d)
  run_num <- rep(xx, length(draws))
  current <- data.frame(draws, run_num)
  results <- rbind(results, current)

}

hist(results$draws, breaks = seq(0,8, by=0.25))
hist(log(results$draws)) # is tri-modal

run.gr <- results %>% group_by(run_num) %>%
  dplyr::summarise(
    gr.an.log.mean = mean(log(draws)),
    gr.an.mean = mean(draws)
  )

hist(run.gr$gr.an.mean) # not tri-modal if we just take the means
hist(run.gr$gr.an.log.mean) # not tri-modal...hmmm
mean(run.gr$gr.an.log.mean)
