# ============================================================
# GO enrichment of the candidate length-at-maturity regions
#
# Input: the candidate regions of 04_gwas_withineco_fst_overlap.R
# (section 3: >= min_adjacent contiguous outlier windows, at least one
# with median p < 0.01; Montejo-Kovacevich et al. 2021), from the runs
# listed in go_runs.
#
# Two gene sets (parameter gene_set):
#   "candidate"  the 1-2 candidate genes of each region: the gene that
#                contains the region's core interval (strongest to second
#                strongest SNP), or the two genes closest to it
#   "region"     every gene overlapping a candidate region (the gene set
#                Montejo-Kovacevich et al. 2021 used for their enrichment)
#
# The unit of analysis is the region cluster, not the gene:
#   - overlapping regions, and regions closer than merge_gap_bp, are
#     merged into one unit, so a cluster of nearby regions (or a region
#     found in both runs) counts once
#   - each GO term is scored by the number of units with at least one
#     gene annotated with it
#   - null = the same regions circularly shifted along their chromosome
#     (one random shift per chromosome), with the genes then chosen
#     exactly as for the real regions. For "candidate", the gene-picking
#     rule is re-applied to each shifted core interval, so the rule's own
#     preferences (long genes contain intervals more often; genes next to
#     gene-poor stretches are picked as "closest" more often) are part of
#     the null and cannot create enrichment.
#
# Multiple testing: every term within the size limits (min_term_size,
# max_term_size, which do not depend on the data) is tested and enters
# the BH correction, whether or not the candidate regions carry it.
# min_regions only decides which terms are REPORTED, after the
# correction. Testing only the terms seen in >= min_regions units and
# then correcting would pick terms by their own outcome: in a simulation
# with no enrichment at all (3,000 terms, ~0.4 units expected per term),
# that gave at least one false discovery at FDR < 0.1 in 94% of
# datasets, against 3.5% when all size-filtered terms are corrected
# together.
#
# GO terms come from the GAWN annotation table (column GeneGo = GO IDs of
# each transcript's UniProt best hit). Only Biological Process terms are
# kept, and each is propagated to its ancestors with GO.db (is_a, part_of
# and the regulates links: "positive regulation of growth" also counts
# for "growth").
#
# GO.db: BiocManager::install("GO.db") if it is not installed.
# ============================================================

library(tidyverse)
library(data.table)
library(GenomicRanges)
library(GO.db)

# ------------------------------------------------------------
# Parameters
# ------------------------------------------------------------
regions_file <- "res/gwas_length_candidate_regions.tsv"        # outputs of 04_gwas_withineco_fst_overlap.R
genes_file   <- "res/gwas_length_candidate_region_genes.tsv"   # (only used to check the gene picking)
annot_file   <- "data/genome_annotation_table.tsv"             # GAWN genome annotation table
gwas_files   <- c("data/gwas_length/summergwas_ALL_CHR.lrt0.gz",
                  "data/gwas_length/fallgwas_ALL_CHR.lrt0.gz")  # only used for chromosome lengths

gene_set      <- "candidate"   # "candidate" or "region" (see header)
go_runs       <- "fall"        # runs whose regions are included: "fall", "summer" or c("summer", "fall")
go_classes    <- c("shared_fst", "shared", "run_fst", "run")   # region classes included
search_bp     <- 1e6     # candidate genes are looked for within this distance of the core interval
merge_gap_bp  <- 1e5     # candidate regions closer than this (bp) are merged into one unit
n_shift       <- 50000   # circular shifts; see the note printed at the end
seed          <- 2026    # random seed for the shifts
min_term_size <- 10      # a term is tested only if it annotates >= this many genes genome-wide
max_term_size <- 500     # ... and <= this many (very broad terms say little)
min_regions   <- 3       # a term is REPORTED only if found in >= this many units (not a test filter)
fdr_report    <- 0.1     # FDR threshold used in the summary

stopifnot(gene_set %in% c("candidate", "region"))

# The run(s) go into the file names, so a fall-only analysis does not
# overwrite the results of the pooled one
run_tag          <- paste(go_runs, collapse = "_")
out_file         <- sprintf("res/gwas_length_go_enrichment_%s_%s.tsv", gene_set, run_tag)
out_genes_file   <- sprintf("res/gwas_length_go_enrichment_%s_%s_genes.tsv", gene_set, run_tag)
out_summary_file <- sprintf("res/gwas_length_go_summary_%s_%s.txt", gene_set, run_tag)

# ============================================================
# Chromosome lengths and candidate regions
# ============================================================

# Chromosome length = furthest SNP position in either GWAS.
# lapply() calls fread() on each file; `select =` is passed on to fread(),
# so only the two columns needed are read.
chr_len <- bind_rows(lapply(gwas_files, fread, select = c("Chromosome", "Position"))) %>%
  mutate(Chromosome = as.character(Chromosome)) %>%
  group_by(Chromosome) %>%
  summarize(len = max(Position), .groups = "drop") %>%
  { setNames(.$len, .$Chromosome) }

reg <- fread(regions_file) %>%
  mutate(Chromosome = as.character(Chromosome)) %>%
  filter(run %in% go_runs, class %in% go_classes, Chromosome %in% names(chr_len))
if (nrow(reg) == 0) stop("No candidate region left after filtering on go_runs and go_classes")

# Units: reduce() merges overlapping regions into one range. With
# min.gapwidth = merge_gap_bp, it also merges regions separated by a gap
# shorter than merge_gap_bp, so a cluster of nearby regions (e.g. several
# regions pointing to the same genes) counts as one locus.
# findOverlaps(..., select = "first") then returns, for each region, the
# number of the merged range it lies in, i.e. one integer per region
# instead of a Hits object. The unit of each region is fixed here and
# reused in every shift: all regions of a chromosome move together, so a
# cluster also moves as one block in the null.
reg_gr   <- GRanges(reg$Chromosome, IRanges(reg$start, reg$end))
merged   <- reduce(reg_gr, min.gapwidth = merge_gap_bp)
reg$unit <- findOverlaps(reg_gr, merged, select = "first")
n_units  <- length(merged)

cat(sprintf("%d candidate regions (%s) -> %d units after merging regions closer than %g kb\n",
            nrow(reg), paste(names(table(reg$run)), table(reg$run), sep = ": ", collapse = ", "),
            n_units, merge_gap_bp / 1e3))

# ============================================================
# GO annotation: Biological Process terms + their ancestors
# ============================================================

annot <- fread(annot_file) %>%
  mutate(ScaffoldName = as.character(ScaffoldName)) %>%
  filter(ScaffoldName %in% names(chr_len))

# Ancestor table from GO.db: one row per (term, ancestor) pair, plus each
# term paired with itself. A gene annotated with "regulation of growth"
# then also counts for "growth" and every other parent up the hierarchy.
bp_anc   <- as.list(GOBPANCESTOR)
anc_long <- tibble(go   = rep(names(bp_anc), lengths(bp_anc)),
                   term = unlist(bp_anc, use.names = FALSE)) %>%
  filter(term != "all") %>%
  bind_rows(tibble(go = names(bp_anc), term = names(bp_anc)))

# Many transcripts share the same GeneGo string (same UniProt hit), so
# the terms are worked out once per distinct string, then looked up.
# str_extract_all() pulls every "GO:" + 7 digits out of a string,
# whatever the separators.
go_strings <- unique(annot$GeneGo)
ids        <- str_extract_all(go_strings, "GO:\\d{7}")

string_term <- tibble(string = rep(seq_along(go_strings), lengths(ids)),
                      go     = unlist(ids)) %>%
  # inner_join keeps only IDs found in the BP ancestor table, which drops
  # Cellular Component / Molecular Function IDs and obsolete terms
  inner_join(anc_long, by = "go", relationship = "many-to-many") %>%
  distinct(string, term)

# Terms are handled as integers 1..n (positions in `terms`): faster than
# strings, and tabulate() can count them directly
terms <- sort(unique(string_term$term))
string_term$term_i <- match(string_term$term, terms)

# terms_of_string[[s]] = integer terms of GO string s. Splitting by a factor
# whose levels are ALL strings keeps an empty entry for strings without BP terms.
terms_of_string <- split(string_term$term_i,
                         factor(string_term$string, levels = seq_along(go_strings)))

# terms_of_tx[[i]] = integer terms of row (transcript) i of annot
terms_of_tx <- unname(terms_of_string[match(annot$GeneGo, go_strings)])

# Term size = number of distinct UniProt hits annotated with the term,
# a proxy for the number of genes
term_size <- tibble(acc = annot$GeneAccession,
                    string = match(annot$GeneGo, go_strings)) %>%
  filter(acc != "-", acc != "") %>%
  distinct() %>%
  inner_join(string_term, by = "string", relationship = "many-to-many") %>%
  distinct(acc, term_i) %>%                       # one row per (UniProt hit, term)
  { tabulate(.$term_i, nbins = length(terms)) }   # rows per term = hits per term

# Terms tested: the size limits only, which do not depend on the data
size_ok <- which(term_size >= min_term_size & term_size <= max_term_size)
if (length(size_ok) == 0) stop("No GO term passes the size limits")

tx_gr <- GRanges(annot$ScaffoldName, IRanges(annot$FromPosition, annot$ToPosition))

# findOverlaps() builds a search index of its subject on every call. With
# ~n_shift calls, it is faster to build it once: GNCList() is that index.
tx_index <- GNCList(tx_gr)

# Candidate genes are picked among named genes only, as in the overlap
# script. `named` = row numbers of those transcripts in annot.
named     <- which(annot$GeneName != "-" & annot$GeneName != "")
named_idx <- GNCList(tx_gr[named])

# ============================================================
# Genes of the (shifted) regions
# ============================================================
#
# Both functions take `sh`, the shift (bp) of each region, and return the
# genes that count for each unit as a list(unit, tx, region): one element
# per (unit, transcript) pair, tx = row of annot. The real data are
# simply the case sh = 0, so observed and null counts go through exactly
# the same code.

reg_L <- unname(chr_len[reg$Chromosome])            # length of each region's chromosome
chr_i <- match(reg$Chromosome, unique(reg$Chromosome))
chr_L <- unname(chr_len[unique(reg$Chromosome)])

# One random shift per chromosome (0 to L - 1), given to every region on
# it, so the regions of a chromosome keep their spacing
draw_shift <- function() floor(runif(length(chr_L)) * chr_L)[chr_i]

# --- "region": every transcript overlapping the region -------
hits_region <- function(sh) {
  s    <- (reg$start - 1 + sh) %% reg_L + 1          # shifted start, wrapping past the end
  e    <- s + (reg$end - reg$start)
  wrap <- e > reg_L                                  # region now runs past the chromosome end:
  gr   <- GRanges(c(reg$Chromosome, reg$Chromosome[wrap]),   # keep the part up to the end and
                  IRanges(c(s, rep(1, sum(wrap))),           # add the overflow at the start
                          c(pmin(e, reg_L), e[wrap] - reg_L[wrap])))
  piece_region <- c(seq_len(nrow(reg)), which(wrap)) # region of each piece
  h <- findOverlaps(gr, tx_index)
  list(unit   = reg$unit[piece_region[queryHits(h)]],
       tx     = subjectHits(h),
       region = piece_region[queryHits(h)])
}

# --- "candidate": the 1-2 candidate genes of the core interval ---
# Same rule as the overlap script: the gene(s) containing the core
# interval, otherwise the two genes closest to it (distance 0 if they
# overlap it), ties broken by distance to the strongest SNP.
hits_candidate <- function(sh) {
  w   <- reg$core_end - reg$core_start
  s   <- (reg$core_start - 1 + sh) %% reg_L + 1
  s   <- pmin(s, reg_L - w)       # a core interval pushed past the chromosome end is moved
  e   <- s + w                    # back so that it ends there (rare: cores are a few kb)
  top <- s + (reg$top_pos - reg$core_start)          # the strongest SNP moves with its core
  
  # Genes within search_bp of each core interval (`+ search_bp` widens a
  # GRanges on both sides). q = which core interval, tx = row of annot.
  h  <- findOverlaps(GRanges(reg$Chromosome, IRanges(s, e)) + search_bp, named_idx)
  dt <- data.table(q = queryHits(h), tx = named[subjectHits(h)])
  
  # data.table syntax: dt[rows, columns, by]. `:=` adds columns to dt
  # in place (no copy, unlike mutate()). s[q] looks up, for each row, the
  # value of s for that row's core interval. This runs n_shift times, and
  # data.table is much faster than dplyr on many small grouped tables.
  dt[, `:=`(gene_start = annot$FromPosition[tx],
            gene_end   = annot$ToPosition[tx],
            name       = annot$GeneName[tx])]
  dt[, `:=`(inside = gene_start <= s[q] & gene_end >= e[q],
            dist   = pmax(0, gene_start - e[q], s[q] - gene_end),   # 0 = overlaps the core
            dtop   = pmax(0, gene_start - top[q], top[q] - gene_end))]
  
  # Closest genes first within each core, with exact ties broken by the
  # same fixed rules as the overlap script (a copy containing the core,
  # then position, then name). findOverlaps() does not return genes in
  # file order, so without these rules tied genes could be kept in a
  # different order from the overlap script.
  dt <- dt[order(q, dist, dtop, !inside, gene_start, gene_end, name)]
  dt <- unique(dt, by = c("q", "name"))      # first (closest) copy of each gene name
  
  # any_inside: does any gene contain this core interval? rowid(q) numbers
  # the rows within each core (1, 2, 3, ...), so rowid(q) <= 2 keeps the
  # two closest genes. Keep the containing gene(s) when there is one,
  # otherwise the two closest.
  dt[, any_inside := any(inside), by = q]
  dt <- dt[(any_inside & inside) | (!any_inside & rowid(q) <= 2)]
  
  list(unit = reg$unit[dt$q], tx = dt$tx, region = dt$q)
}

hit_fun <- switch(gene_set, candidate = hits_candidate, region = hits_region)

# Returns an integer vector, element t = number of units with >= 1 gene
# annotated with term t
count_terms <- function(hits) {
  tx_by_unit   <- split(hits$tx, hits$unit)                     # genes of each unit
  term_by_unit <- lapply(tx_by_unit, function(tx)
    unique(unlist(terms_of_tx[tx], use.names = FALSE)))         # each term once per unit
  # as.integer() turns NULL (no unit hit any gene) into integer(0)
  tabulate(as.integer(unlist(term_by_unit, use.names = FALSE)), nbins = length(terms))
}

# ============================================================
# Observed counts
# ============================================================

obs_hits <- hit_fun(rep(0, nrow(reg)))
obs      <- count_terms(obs_hits)

obs_genes <- tibble(region_id = reg$region_id[obs_hits$region],
                    unit      = obs_hits$unit,
                    tx        = obs_hits$tx,
                    gene      = annot$GeneName[obs_hits$tx]) %>%
  distinct()

cat(sprintf("Gene set '%s': %d distinct named genes in %d of %d units\n",
            gene_set, n_distinct(obs_genes$gene[!obs_genes$gene %in% c("-", "")]),
            n_distinct(obs_genes$unit), n_units))

# A unit without any gene cannot carry a term. In the null, shifted
# units that land in gene deserts (no gene within search_bp, or no gene
# overlapping the region) are counted below, so you can see whether that
# happens much more often than in the real data.
n_empty_obs <- n_units - n_distinct(obs_hits$unit)
if (n_empty_obs > 0)
  warning(sprintf("%d units have no gene (no named gene within %g kb for 'candidate')",
                  n_empty_obs, search_bp / 1e3))

# Gene names found in more than one unit. Unnamed transcripts ("-") are
# left out: they all share that name without being one gene.
#   same chromosome      the units are probably one locus that
#                        merge_gap_bp did not join, so a term carried by
#                        that gene is counted twice: increase merge_gap_bp
#   different chromosomes  most likely two paralogues with the same best
#                        UniProt hit (salmonid whole-genome duplication):
#                        two separate loci, correctly counted twice
shared_genes <- obs_genes %>%
  filter(!gene %in% c("-", "")) %>%
  mutate(Chromosome = annot$ScaffoldName[tx]) %>%
  distinct(unit, region_id, gene, Chromosome) %>%
  group_by(gene) %>%
  filter(n_distinct(unit) > 1) %>%
  summarize(units       = paste(sort(unique(unit)), collapse = ", "),
            region_ids  = paste(sort(unique(region_id)), collapse = ", "),
            chromosomes = paste(sort(unique(Chromosome)), collapse = ", "),
            same_chromosome = n_distinct(Chromosome) == 1,
            .groups = "drop")
if (any(shared_genes$same_chromosome))
  warning(sprintf("%d gene names occur in more than one unit on the same chromosome (printed below): if those units are one locus, increase merge_gap_bp.",
                  sum(shared_genes$same_chromosome)))
if (nrow(shared_genes) > 0) {
  cat("Gene names found in more than one unit:\n")
  print(shared_genes, n = Inf, width = Inf)
} else {
  cat("No gene name is shared by two units\n")
}

# The candidate genes picked here should be those of the overlap script.
# A difference means the two scripts read different annotation or regions
# files, or a core interval is more than search_bp from any gene.
if (gene_set == "candidate" && file.exists(genes_file)) {
  # as_tibble(): fread() returns a data.table, and dplyr verbs keep that
  # class, so print(diff, n = Inf) below would go to data.table's print
  # method, which has no `n` argument
  file_genes <- fread(genes_file) %>%
    as_tibble() %>%
    filter(region_id %in% reg$region_id, !is.na(gene_name)) %>%
    distinct(region_id, gene = gene_name)
  mine <- distinct(obs_genes, region_id, gene)
  # anti_join(x, y) = rows of x with no match in y
  diff <- bind_rows(
    anti_join(file_genes, mine, by = c("region_id", "gene")) %>% mutate(only_in = "overlap script"),
    anti_join(mine, file_genes, by = c("region_id", "gene")) %>% mutate(only_in = "this script")
  ) %>%
    arrange(region_id)
  if (nrow(diff) > 0) {
    warning(sprintf("%d region-gene pairs differ from %s (printed below). Re-run the overlap script if it predates the tie-breaking rules.",
                    nrow(diff), genes_file))
    print(diff, n = Inf)
  } else {
    cat("Candidate genes identical to those of the overlap script\n")
  }
}

obs_ok <- obs[size_ok]
cat(sprintf("%d BP terms tested (size %d-%d genes); %d of them found in >= %d units\n",
            length(size_ok), min_term_size, max_term_size, sum(obs_ok >= min_regions), min_regions))

# ============================================================
# Null: circular shifts
# ============================================================
#
# Only two running totals are kept per term, so memory does not grow
# with n_shift (a terms x shifts matrix would reach several GB):
#   n_ge    number of shifts whose count is >= the observed count
#   sum_nl  sum of the null counts; / n_shift = expected count
# `x >= obs_ok` is a logical vector. Adding it to an integer vector counts
# TRUE as 1 and FALSE as 0, so n_ge goes up by one for every term whose
# shifted count reaches the observed one.
#
# The seed is set right before the only random step of the script.

set.seed(seed)
n_ge      <- integer(length(size_ok))
sum_nl    <- numeric(length(size_ok))
n_empty   <- 0                         # shifted units without any gene, summed over shifts
every     <- max(1, n_shift %/% 10)    # progress message every 10% (at least every shift)
t0        <- Sys.time()

cat(sprintf("Running %s circular shifts...\n", format(n_shift, big.mark = ",")))
for (i in seq_len(n_shift)) {
  hits    <- hit_fun(draw_shift())
  x       <- count_terms(hits)[size_ok]
  n_ge    <- n_ge + (x >= obs_ok)
  sum_nl  <- sum_nl + x
  n_empty <- n_empty + (n_units - length(unique(hits$unit)))
  if (i %% every == 0)
    cat(sprintf("  %s / %s shifts (%.1f min)\n", format(i, big.mark = ","),
                format(n_shift, big.mark = ","),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

cat(sprintf("Units without any gene: %d of %d in the data; %.2f of %d on average in the shifts\n",
            n_empty_obs, n_units, n_empty / n_shift, n_units))

# ============================================================
# Results
# ============================================================

# Genes and regions behind each term in the real data
unit_info <- reg %>%
  group_by(unit) %>%
  summarize(Chromosome = Chromosome[1],
            unit_start = min(start),
            unit_end   = max(end),
            region_ids = paste(region_id, collapse = ", "),
            .groups = "drop")

term_gene_region <- obs_genes %>%
  mutate(term_i = terms_of_tx[tx]) %>%             # list column: the terms of each gene
  unnest(term_i) %>%                               # one row per (gene, term)
  filter(term_i %in% size_ok) %>%
  transmute(go_id = terms[term_i], unit, region_id, gene,
            accession = annot$GeneAccession[tx]) %>%
  distinct() %>%
  left_join(unit_info, by = "unit")

term_genes <- term_gene_region %>%
  group_by(go_id) %>%
  summarize(genes = paste(sort(unique(gene)), collapse = ", "), .groups = "drop")

# One row per tested term. tibble() can use a column defined earlier in
# the same call (fold uses n_units and expected). fold is NA when no
# shift ever gave the term a unit (expected = 0), instead of Inf or NaN.
# The BH correction runs over ALL tested terms; `reported` only flags the
# terms found in >= min_regions units, for the summary.
res <- tibble(
  go_id          = terms[size_ok],
  term           = unname(Term(terms[size_ok])),
  n_genes_genome = term_size[size_ok],
  n_units        = obs_ok,                         # units (merged regions) with the term
  expected       = sum_nl / n_shift,
  fold           = if_else(expected > 0, n_units / expected, NA_real_),
  emp_p          = (1 + n_ge) / (n_shift + 1),
  fdr            = p.adjust(emp_p, method = "BH"),
  reported       = n_units >= min_regions
) %>%
  left_join(term_genes, by = "go_id") %>%
  arrange(emp_p, desc(fold))

term_gene_region <- term_gene_region %>%
  left_join(res %>% dplyr::select(go_id, term, emp_p, fdr), by = "go_id") %>%
  arrange(emp_p, go_id, Chromosome, unit_start)

fwrite(res, out_file, sep = "\t")
fwrite(term_gene_region, out_genes_file, sep = "\t")

# ============================================================
# Summary (console + text file)
# ============================================================
# The smallest possible empirical p is 1 / (n_shift + 1). BH multiplies
# p by (number of terms tested / rank), so with too few shifts no term
# can reach a low FDR: a single top term can at best reach
# n_tested / (n_shift + 1). For it to reach FDR = 0.05, n_shift must be
# at least 20 x n_tested.
#
# sink(file, split = TRUE) sends everything printed to the console to
# the file as well; on.exit(sink()) undoes it when the function ends,
# even after an error.
write_summary <- function(file) {
  sink(file, split = TRUE)
  on.exit(sink(), add = TRUE)
  old_opt <- options(width = 200)
  on.exit(options(old_opt), add = TRUE)
  
  n_tested <- nrow(res)
  cat(sprintf("\nGO enrichment (Biological Process), gene set '%s', runs: %s\n", gene_set, run_tag))
  cat(sprintf("%d candidate regions -> %d units; %d BP terms tested (size %d-%d), BH over all of them\n",
              nrow(reg), n_units, n_tested, min_term_size, max_term_size))
  cat(sprintf("%s circular shifts: smallest possible P = %.1e; smallest possible FDR for a single top term = %.3f%s\n",
              format(n_shift, big.mark = ","), 1 / (n_shift + 1), n_tested / (n_shift + 1),
              if (n_tested / (n_shift + 1) > 0.05)
                sprintf(" (raise n_shift to >= %s to allow FDR 0.05)", format(20 * n_tested, big.mark = ","))
              else ""))
  cat(sprintf("Terms with FDR < %g: %d (%d of them found in >= %d units)\n",
              fdr_report, sum(res$fdr < fdr_report),
              sum(res$fdr < fdr_report & res$reported), min_regions))
  
  cat(sprintf("\nTop terms found in >= %d units:\n", min_regions))
  res %>%
    filter(reported) %>%
    dplyr::select(go_id, term, n_units, expected, fold, emp_p, fdr, genes) %>%
    mutate(across(c(expected, fold), ~ round(.x, 2)),
           across(c(emp_p, fdr), ~ signif(.x, 2))) %>%
    head(20) %>%
    print(n = Inf, width = Inf)
  invisible(NULL)
}

write_summary(out_summary_file)
cat("\nWritten:", out_file, out_genes_file, out_summary_file, sep = "\n  ")
