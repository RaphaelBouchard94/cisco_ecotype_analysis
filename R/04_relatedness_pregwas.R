library(tidyverse)


res <- read_tsv("data/ngsrelate/notsummer.rel")   # use a/b instead of ida/idb if you didn't run with -z

flagged <- res |>
  filter(KING > 0.177 | rab > 0.354) |>
  select(a, b, nSites, rab, theta, KING, R0, R1, Fa, Fb, J9, J8, J7) |>
  arrange(desc(KING))

per_ind <- bind_rows(
  res |> transmute(id = a, F = Fa, rab, KING),
  res |> transmute(id = b, F = Fb, rab, KING)
) |>
  group_by(id) |>
  summarise(median_F = median(F), mean_rab = mean(rab),
            mean_KING = mean(KING), n_flagged = sum(rab > 0.177)) |>
  arrange(desc(median_F))

ggplot(res, aes(KING, rab)) + geom_point(alpha = 0.5)


res <- res |>
  mutate(rel = case_when(
    KING > 0.354              ~ "duplicate",
    KING > 0.177 & R0 < 0.02  ~ "PO",
    KING > 0.177              ~ "FS",
    KING > 0.0884             ~ "2nd degree",
    KING > 0.0442             ~ "3rd degree",
    TRUE                      ~ "unrelated"),
    R0_above_0.5  = R0 > 0.5,                      # possible different origin
    freq_vs_free  = abs(theta - KING) > 0.05)      # rab-based and frequency-free estimates disagree

count(res, rel)

ggplot(res, aes(R1, R0, colour = rel)) + geom_point() +
  annotate("point", x = 0.5, y = 0, shape = 4, size = 4)      # expected PO position
ggplot(res, aes(R1, KING, colour = rel)) + geom_point() +
  annotate("point", x = 0.5, y = 0.25, shape = 4, size = 4)

pairs <- res |> filter(KING > 0.177) |> select(a, b)


removed <- character()
while (nrow(pairs) > 0) {
  n <- table(c(pairs$a, pairs$b))
  candidates <- names(n)[n == max(n)]
  drop <- if (exists("depth")) candidates[which.min(depth[candidates])] else candidates[1]
  removed <- c(removed, drop)
  pairs <- filter(pairs, a != drop, b != drop)
}

removed




