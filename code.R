# Packages ---------------------------------------
load_or_install <- function(pkgs) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) install.packages(missing_pkgs, dependencies = TRUE)
  invisible(lapply(pkgs, \(p) suppressPackageStartupMessages(library(p, character.only = TRUE))))
}

pkgs <- c(
  "readxl", "tidyverse", "vegan", "MASS", "performance", "cowplot",
  "broom", "spdep", "RColorBrewer", "ggdist", "purrr",
  "sjPlot", "patchwork"
)

load_or_install(pkgs)

# Getting data ready -------------------------------------------
data <- read_excel("data.xlsx")

assemblages <- data[, 13:32]
env <- data[, c(1, 4, 9:12)]

thermal_traits <- read_excel("thermal.xlsx") %>%
  rename(species = Species) %>%
  mutate(across(2:last_col(), as.numeric))

morphology <- read_excel("morphology.xlsx")

morphology_clean <- morphology %>%
  mutate(across(2:last_col(), as.numeric)) %>%
  filter(!if_any(where(is.numeric), is.na))

# Calculating and combining traits ------------------------------
morphology_clean <- morphology_clean %>%
  mutate(
    compression_index = Bd / Bw,
    relative_depth = Bd / SL,
    index_ventral_flat = Mmd / Bd,
    relative_eye_posit = Eh / Hd,
    finenness = (SL / sqrt(Bd)) * Bw,
    relative_mouth_width = Mw / SL,
    relative_eye_size = Ed / Hd
  )

morphology_traits <- morphology_clean %>%
  group_by(species) %>%
  summarise(
    across(
      c(
        compression_index,
        relative_depth,
        index_ventral_flat,
        relative_eye_posit,
        finenness,
        relative_mouth_width,
        relative_eye_size
      ),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

traits_combined <- morphology_traits %>%
  left_join(thermal_traits, by = "species") %>%
  dplyr::select(-any_of(c("CTMax_CV")))

# Species richness ----------------------------------------------
env$S <- specnumber(assemblages)

# Creating predictors -------------------------------------------
data$maximum_length_m <- as.numeric(data$maximum_length_m)
data$maximum_width_m <- as.numeric(data$maximum_width_m)
data$maximum_depth_cm <- as.numeric(data$maximum_depth_cm)

data$depth_m <- data$maximum_depth_cm / 100

env$Volume <- (2/3) * pi *
  (data$maximum_length_m/2) *
  (data$maximum_width_m/2) *
  data$depth_m

env <- env %>%
  mutate(across(3:last_col(), as.numeric))

# Environmental stress indexes ---------------------------------
pca_env <- prcomp(
  env %>% dplyr::select(Temp, DO, pH),
  center = TRUE,
  scale. = TRUE
)

scores <- as.data.frame(pca_env$x)

env <- env %>%
  mutate(
    Stress_index_1 = scores$PC1,
    Stress_index_2 = scores$PC2
  )

loadings <- as.data.frame(pca_env$rotation)
loadings
summary(pca_env)

# Trait similarity ----------------------------------------------
traits_matrix <- traits_combined %>%
  column_to_rownames("species") %>%
  as.matrix()

traits_scaled <- scale(traits_matrix)

dist_matrix <- dist(traits_scaled, method = "euclidean")
dist_mat <- as.matrix(dist_matrix)

traits_similarity <- apply(assemblages, 1, function(abund_row) {
  sp_present <- names(abund_row)[abund_row > 0]
  sp_present <- intersect(sp_present, rownames(dist_mat))
  
  if (length(sp_present) <= 1) return(0)
  
  dsub <- dist_mat[sp_present, sp_present, drop = FALSE]
  dvals <- dsub[upper.tri(dsub)]
  
  mean(1 / (1 + dvals), na.rm = TRUE)
})

env$traits_similarity <- traits_similarity

# Terrestrial dispersal capacity --------------------------------
amphibious_vec <- c(
  "Atlantirivulus peruibensis",
  "Leptopanchax itanhaensis",
  "Callichthys callichthys",
  "Synbranchus marmoratus",
  "Hoplosternum littorale"
)

env$terrestrial_dispersal_capacity <- apply(assemblages, 1, function(abund_row) {
  sp_present <- names(abund_row)[abund_row > 0]
  as.integer(any(sp_present %in% amphibious_vec))
})

# Aquatic dispersal capacity ------------------------------------
fin_map <- setNames(traits_combined$finenness, traits_combined$species)

env$aquatic_dispersal_capacity <- apply(assemblages, 1, function(abund_row) {
  sp_present <- names(abund_row)[abund_row > 0]
  sp_present <- intersect(sp_present, names(fin_map))
  
  if (length(sp_present) == 0) return(NA_real_)
  
  mean(fin_map[sp_present], na.rm = TRUE)
})

# Scale predictors ----------------------------------------------
env_scaled <- env %>%
  mutate(
    across(
      where(is.numeric) & !any_of(c("S", "terrestrial_dispersal_capacity")),
      ~ scale(.x)[, 1]
    )
  )

env_filt <- env %>%
  filter(S > 1)

env_filt_scaled <- env_filt %>%
  mutate(
    across(
      where(is.numeric) & !any_of(c("S", "terrestrial_dispersal_capacity")),
      ~ scale(.x)[, 1]
    )
  )

# Models --------------------------------------------------------
m1 <- glm.nb(S ~ stream_distance + Volume + Period, data = env_scaled)

m2 <- glm.nb(S ~ Stress_index_1 + Stress_index_2 + Period, data = env_scaled)

m3 <- glm(
  S ~ traits_similarity + aquatic_dispersal_capacity +
    terrestrial_dispersal_capacity + Period,
  data = env_filt_scaled,
  family = poisson(link = "log")
)

summary(m1); r2(m1)
summary(m2); r2(m2)
summary(m3); r2(m3)

# Figure 1A -----------------------------------------------------
model_list <- list(
  "Island biogeography" = m1,
  "Environmental filtering" = m2,
  "Metacommunity" = m3
)

coef_df <- purrr::imap_dfr(model_list, ~{
  broom::tidy(.x, conf.int = TRUE) %>%
    filter(term != "(Intercept)",
           term != "PeriodWet Period") %>%
    mutate(model = .y)
})

coef_df$model <- factor(coef_df$model, levels = names(model_list))

coef_df <- coef_df %>%
  mutate(
    term = dplyr::recode(
      term,
      stream_distance = "Stream distance",
      Volume = "Volume",
      Stress_index_1 = "Stress index 1",
      Stress_index_2 = "Stress index 2",
      traits_similarity = "Trait similarity",
      aquatic_dispersal_capacity = "Aquatic dispersal",
      terrestrial_dispersal_capacity = "Terrestrial dispersal"
    )
  ) %>%
  mutate(term_id = paste(model, term, sep = "___"))

term_order <- coef_df %>%
  mutate(model_num = as.numeric(model)) %>%
  arrange(model_num) %>%
  pull(term_id)

coef_df$term_id <- factor(coef_df$term_id, levels = rev(term_order))

bg_df <- coef_df %>%
  distinct(model, term_id) %>%
  mutate(y = as.numeric(term_id)) %>%
  group_by(model) %>%
  summarise(
    ymin = min(y) - 0.5,
    ymax = max(y) + 0.5,
    .groups = "drop"
  ) %>%
  mutate(xmin = -Inf, xmax = Inf)

label_df <- coef_df %>%
  distinct(term_id, term)

framework_cols <- c(
  "Island biogeography"      = "#9C6644",
  "Environmental filtering" = "#DDA15E", 
  "Metacommunity"           = "#283618"
)

pal <- framework_cols

fig1_A <- ggplot(coef_df, aes(x = estimate, y = term_id, color = model)) +
  geom_rect(
    data = bg_df,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = model
    ),
    inherit.aes = FALSE,
    alpha = 0.18,
    color = NA
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  geom_errorbar(
    aes(
      xmin = conf.low,
      xmax = conf.high
    ),
    orientation = "y",
    width = 0.18,
    linewidth = 1.6,
    color = "black"
  ) +
  geom_errorbar(
    aes(
      xmin = conf.low,
      xmax = conf.high,
      color = model
    ),
    orientation = "y",
    width = 0.12,
    linewidth = 1
  ) +
  geom_point(
    aes(fill = model),
    shape = 21,
    size = 4.5,
    stroke = 1.1,
    color = "black"
  ) +
  scale_y_discrete(
    labels = setNames(label_df$term, label_df$term_id)
  ) +
  scale_x_continuous(
    limits = c(-1.2, 1.2)
  ) +
  scale_color_manual(
    values = framework_cols
  ) +
  scale_fill_manual(
    values = framework_cols
  ) +
  theme_classic(base_size = 20) +
  labs(
    x = "Coefficient estimate (± 95% CI)",
    y = NULL,
    fill = "Framework"
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(
        shape = 21,
        size = 4,
        color = "black"
      )
    ),
    color = "none"
  ) +
  theme(
    legend.title = element_text(
      face = "bold",
      size = 16
    ),
    legend.text = element_text(
      face = "bold",
      size = 14
    ),
    legend.position = "right",
    legend.justification = "top",
    legend.box.just = "top"
  )

fig1_A

# Figure 1B -----------------------------------------------------
r2_df <- purrr::imap_dfr(model_list, ~{
  tibble(
    model = .y,
    R2 = unname(performance::r2(.x)$R2)
  )
}) %>%
  arrange(desc(R2)) %>%
  mutate(model = factor(model, levels = model))

fig1_B <- ggplot(r2_df, aes(x = model, y = R2, fill = model)) +
  geom_col(
    color = "black",
    width = 0.93,
    linewidth = 1,
    alpha = 0.7
  ) +
  scale_fill_manual(values = framework_cols) +
  coord_flip() +
  theme_classic(base_size = 20) +
  labs(x = NULL, y = expression(R^2)) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 0.50),
    breaks = seq(0, 0.50, by = 0.1)
  ) +
  theme(
    axis.title.x = element_text(face = "bold"),
    legend.position = "none"
  )

fig1_B

# Spatial autocorrelation ---------------------------------------
# Geographic coordinates were converted to decimal degrees, 
# and duplicate points were slightly jittered (±0.00001°) to prevent 
# overlap while preserving spatial relationships

env$coord <- data$coordinates

dms_to_decimal <- function(x) {
  x <- str_trim(x)
  parts <- str_split_fixed(x, "\\s+", 2)
  lat_s <- parts[, 1]
  lon_s <- parts[, 2]
  
  parse_one <- function(s) {
    s <- str_replace_all(s, "\\s+", "")
    m <- str_match(s, "^(\\d+)°(\\d+)'([0-9.]+)\"?([NSEW])$")
    
    deg <- as.numeric(m[, 2])
    min <- as.numeric(m[, 3])
    sec <- as.numeric(m[, 4])
    hem <- m[, 5]
    
    dec <- deg + min / 60 + sec / 3600
    dec[hem %in% c("S", "W")] <- -dec[hem %in% c("S", "W")]
    dec
  }
  
  tibble(
    coord = x,
    lat = parse_one(lat_s),
    lon = parse_one(lon_s)
  )
}

coords_df <- dms_to_decimal(env$coord)

env_spatial <- env %>%
  mutate(row_id = row_number()) %>%
  bind_cols(coords_df %>% dplyr::select(lon, lat)) %>%
  group_by(lon, lat) %>%
  mutate(
    n_dup = n(),
    dup_id = row_number(),
    lon_jit = if (n() > 1) lon + seq(-0.00001, 0.00001, length.out = n())[dup_id] else lon,
    lat_jit = if (n() > 1) lat + seq(-0.00001, 0.00001, length.out = n())[dup_id] else lat
  ) %>%
  ungroup()

env_spatial_scaled <- env_spatial %>%
  mutate(
    across(
      where(is.numeric) & !any_of(c(
        "S", "terrestrial_dispersal_capacity",
        "row_id", "n_dup", "dup_id",
        "lon", "lat", "lon_jit", "lat_jit"
      )),
      ~ scale(.x)[, 1]
    )
  )

env_spatial_filt_scaled <- env_spatial %>%
  filter(S > 1) %>%
  mutate(
    across(
      where(is.numeric) & !any_of(c(
        "S", "terrestrial_dispersal_capacity",
        "row_id", "n_dup", "dup_id",
        "lon", "lat", "lon_jit", "lat_jit"
      )),
      ~ scale(.x)[, 1]
    )
  )

m1_sp <- glm.nb(
  S ~ stream_distance + Volume + Period,
  data = env_spatial_scaled
)

m2_sp <- glm.nb(
  S ~ Stress_index_1 + Stress_index_2 + Period,
  data = env_spatial_scaled
)

m3_sp <- glm(
  S ~ traits_similarity +
    aquatic_dispersal_capacity +
    terrestrial_dispersal_capacity +
    Period,
  data = env_spatial_filt_scaled,
  family = poisson(link = "log")
)

moran_residuals <- function(mod, dat, k) {
  
  dat2 <- dat %>%
    filter(
      is.finite(lon_jit),
      is.finite(lat_jit)
    )
  
  r <- residuals(mod, type = "pearson")
  
  coords <- as.matrix(
    dat2[, c("lon_jit", "lat_jit")]
  )
  
  nb <- knn2nb(
    knearneigh(coords, k = k)
  )
  
  lw <- nb2listw(
    nb,
    style = "W",
    zero.policy = TRUE
  )
  
  mi <- moran.test(
    r,
    lw,
    zero.policy = TRUE
  )
  
  tibble(
    k = k,
    Moran_I = unname(
      mi$estimate[["Moran I statistic"]]
    ),
    p_value = mi$p.value,
    n = length(r)
  )
}

model_list_spatial <- list(
  "Island biogeography" = list(
    model = m1_sp,
    data = env_spatial_scaled
  ),
  
  "Environmental filtering" = list(
    model = m2_sp,
    data = env_spatial_scaled
  ),
  
  "Metacommunity" = list(
    model = m3_sp,
    data = env_spatial_filt_scaled
  )
)

k_values <- c(4, 7, 12, 15)

spatial_tests <- purrr::imap_dfr(
  model_list_spatial,
  function(x, model_name) {
    
    purrr::map_dfr(
      k_values,
      ~ moran_residuals(
        mod = x$model,
        dat = x$data,
        k = .x
      )
    ) %>%
      mutate(
        Model = model_name,
        .before = 1
      )
  }
)

spatial_tests

spatial_tests %>% arrange(p_value)

range(spatial_tests$Moran_I)
range(spatial_tests$p_value)

# Figure S3 -----------------------------------------------------
m1_s <- glm.nb(
  S ~ stream_distance + Volume,
  data = env_filt_scaled
)

m2_s <- glm.nb(
  S ~ Stress_index_1 + Stress_index_2,
  data = env_filt_scaled
)

m3_s <- glm(
  S ~ traits_similarity + aquatic_dispersal_capacity +
    terrestrial_dispersal_capacity,
  data = env_filt_scaled,
  family = poisson(link = "log")
)

model_list_s <- list(
  "Island biogeography" = m1_s,
  "Environmental filtering" = m2_s,
  "Metacommunity" = m3_s
)

coef_df_s <- purrr::imap_dfr(model_list_s, ~{
  broom::tidy(.x, conf.int = TRUE) %>%
    filter(term != "(Intercept)") %>%
    mutate(model = .y)
})

coef_df_s$model <- factor(coef_df_s$model, levels = names(model_list_s))

coef_df_s <- coef_df_s %>%
  mutate(
    term = dplyr::recode(
      term,
      stream_distance = "Stream distance",
      Volume = "Volume",
      Stress_index_1 = "Stress index 1",
      Stress_index_2 = "Stress index 2",
      traits_similarity = "Trait similarity",
      aquatic_dispersal_capacity = "Aquatic dispersal",
      terrestrial_dispersal_capacity = "Terrestrial dispersal"
    )
  ) %>%
  mutate(term_id = paste(model, term, sep = "___"))

term_order_s <- coef_df_s %>%
  mutate(model_num = as.numeric(model)) %>%
  arrange(model_num) %>%
  pull(term_id)

coef_df_s$term_id <- factor(coef_df_s$term_id, levels = rev(term_order_s))

bg_df_s <- coef_df_s %>%
  distinct(model, term_id) %>%
  mutate(y = as.numeric(term_id)) %>%
  group_by(model) %>%
  summarise(
    ymin = min(y) - 0.5,
    ymax = max(y) + 0.5,
    .groups = "drop"
  ) %>%
  mutate(xmin = -Inf, xmax = Inf)

label_df_s <- coef_df_s %>%
  distinct(term_id, term)

fig_s3_A <- ggplot(coef_df_s, aes(x = estimate, y = term_id, color = model)) +
  geom_rect(
    data = bg_df_s,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = model
    ),
    inherit.aes = FALSE,
    alpha = 0.18,
    color = NA
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  geom_errorbar(
    aes(
      xmin = conf.low,
      xmax = conf.high
    ),
    orientation = "y",
    width = 0.18,
    linewidth = 1.6,
    color = "black"
  ) +
  geom_errorbar(
    aes(
      xmin = conf.low,
      xmax = conf.high,
      color = model
    ),
    orientation = "y",
    width = 0.12,
    linewidth = 1
  ) +
  geom_point(
    aes(fill = model),
    shape = 21,
    size = 4.5,
    stroke = 1.1,
    color = "black"
  ) +
  scale_y_discrete(
    labels = setNames(label_df_s$term, label_df_s$term_id)
  ) +
  scale_x_continuous(
    limits = c(-1.2, 1.2)
  ) +
  scale_color_manual(
    values = framework_cols
  ) +
  scale_fill_manual(
    values = framework_cols
  ) +
  theme_classic(base_size = 20) +
  labs(
    x = "Coefficient estimate (± 95% CI)",
    y = NULL,
    fill = "Framework"
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(
        shape = 21,
        size = 4,
        color = "black"
      )
    ),
    color = "none"
  ) +
  theme(
    legend.title = element_text(face = "bold", size = 16),
    legend.text = element_text(face = "bold", size = 14),
    legend.position = "right",
    legend.justification = "top",
    legend.box.just = "top"
  )

fig_s3_A


# Figure S3B ----------------------------------------------------
r2_df_s <- purrr::imap_dfr(model_list_s, ~{
  tibble(
    model = .y,
    R2 = unname(performance::r2(.x)$R2)
  )
}) %>%
  arrange(desc(R2)) %>%
  mutate(model = factor(model, levels = model))

fig_s3_B <- ggplot(r2_df_s, aes(x = model, y = R2, fill = model)) +
  geom_col(
    color = "black",
    width = 0.93,
    linewidth = 1,
    alpha = 0.7
  ) +
  scale_fill_manual(values = framework_cols) +
  coord_flip() +
  theme_classic(base_size = 20) +
  labs(
    x = NULL,
    y = expression(R^2)
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 0.50),
    breaks = seq(0, 0.50, by = 0.1)
  ) +
  theme(
    axis.title.x = element_text(face = "bold"),
    legend.position = "none"
  )

fig_s3_B


# Figure S3C ----------------------------------------------------

aic_df_s <- purrr::imap_dfr(model_list_s, ~{
  tibble(
    model = .y,
    AIC = AIC(.x)
  )
}) %>%
  mutate(delta_AIC = AIC - min(AIC)) %>%
  arrange(delta_AIC) %>%
  mutate(model = factor(model, levels = model))

fig_s3_C <- ggplot(aic_df_s, aes(x = model, y = delta_AIC, fill = model)) +
  geom_col(
    color = "black",
    width = 0.93,
    linewidth = 1,
    alpha = 0.7
  ) +
  scale_fill_manual(values = framework_cols) +
  coord_flip() +
  theme_classic(base_size = 20) +
  labs(
    x = NULL,
    y = expression(Delta*AIC)
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0,20)
  ) +
  theme(
    axis.title.x = element_text(face = "bold"),
    legend.position = "none"
  )

fig_s3_C


# Combine and save ----------------------------------------------

fig_s3 <- fig_s3_A + (fig_s3_B / fig_s3_C) +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")"
  )


fig_s3

ggsave(
  "Figure_S3.jpg",
  fig_s3,
  width = 15,
  height = 7.5
)

# Figure 1 Complete ------------------------
fig1 <- fig1_A + (fig1_B / fig_s3_C) +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")"
  )

fig1

ggsave("Figure_1.jpg", fig1, width = 15, height = 7.5)


# Null models: only metacommunity -------------------------------
only_scale <- function(x) scale(x)[, 1]

calc_traits_similarity_sp <- function(sp_present) {
  sp_present <- intersect(sp_present, rownames(dist_mat))
  if (length(sp_present) <= 1) return(0)
  
  dsub <- dist_mat[sp_present, sp_present, drop = FALSE]
  dvals <- dsub[upper.tri(dsub)]
  
  mean(1 / (1 + dvals), na.rm = TRUE)
}

calc_aquatic_disp <- function(sp_present) {
  sp_present <- intersect(sp_present, names(fin_map))
  if (length(sp_present) == 0) return(NA_real_)
  mean(fin_map[sp_present], na.rm = TRUE)
}

calc_terrestrial_disp <- function(sp_present) {
  if (length(sp_present) == 0) return(0)
  as.integer(any(sp_present %in% amphibious_vec))
}

species_pool <- colnames(assemblages)
site_richness <- env$S
n_reps <- 1000

generate_null_assemblage <- function(assemblages, site_richness, species_pool) {
  null_mat <- matrix(
    0,
    nrow = nrow(assemblages),
    ncol = length(species_pool),
    dimnames = list(rownames(assemblages), species_pool)
  )
  
  for (i in seq_len(nrow(assemblages))) {
    S_i <- site_richness[i]
    if (is.na(S_i) || S_i == 0) next
    
    sp_i <- sample(species_pool, size = S_i, replace = FALSE)
    null_mat[i, sp_i] <- 1
  }
  
  as.data.frame(null_mat)
}

build_null_env <- function(null_assemblage, env) {
  sp_list_by_site_null <- apply(null_assemblage, 1, function(abund_row) {
    names(abund_row)[abund_row > 0]
  })
  
  env_null <- env %>%
    dplyr::select(S)
  
  env_null$traits_similarity <- vapply(
    sp_list_by_site_null,
    calc_traits_similarity_sp,
    numeric(1)
  )
  
  env_null$aquatic_dispersal_capacity <- vapply(
    sp_list_by_site_null,
    calc_aquatic_disp,
    numeric(1)
  )
  
  env_null$terrestrial_dispersal_capacity <- vapply(
    sp_list_by_site_null,
    calc_terrestrial_disp,
    numeric(1)
  )
  
  env_null
}

prepare_env_null_scaled_filt <- function(env_null) {
  env_null %>%
    filter(S > 1) %>%
    mutate(
      across(
        where(is.numeric) & !any_of(c("S", "terrestrial_dispersal_capacity")),
        ~ only_scale(.x)
      )
    )
}

fit_null_model <- function(env_null_scaled_filt) {
  glm(
    S ~ traits_similarity + aquatic_dispersal_capacity +
      terrestrial_dispersal_capacity,
    data = env_null_scaled_filt,
    family = poisson(link = "log")
  )
}

set.seed(123)

null_results <- vector("list", n_reps)

for (r in seq_len(n_reps)) {
  null_assemblage <- generate_null_assemblage(
    assemblages = assemblages,
    site_richness = site_richness,
    species_pool = species_pool
  )
  
  env_null <- build_null_env(null_assemblage, env)
  env_null_scaled_filt <- prepare_env_null_scaled_filt(env_null)
  
  null_results[[r]] <- suppressWarnings(
    fit_null_model(env_null_scaled_filt)
  )
}

coef_null_df <- purrr::imap_dfr(null_results, ~{
  broom::tidy(.x, conf.int = TRUE) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      replicate = .y,
      model = "Metacommunity"
    )
})

coef_null_sum <- coef_null_df %>%
  group_by(model, term) %>%
  summarise(
    mean_est = mean(estimate, na.rm = TRUE),
    ci_low = quantile(estimate, 0.025, na.rm = TRUE),
    ci_high = quantile(estimate, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

coef_emp_df <- broom::tidy(m3) %>%
  filter(term != "(Intercept)",
         term != "PeriodWet Period") %>%
  mutate(model = "Metacommunity")

# Figure 2 ------------------------------------------------------
coef_keep <- c(
  "traits_similarity",
  "aquatic_dispersal_capacity",
  "terrestrial_dispersal_capacity"
)

coef_labels <- c(
  traits_similarity = "Trait similarity",
  aquatic_dispersal_capacity = "Aquatic dispersal",
  terrestrial_dispersal_capacity = "Terrestrial dispersal"
)

coef_null_df2 <- coef_null_df %>%
  filter(term %in% coef_keep) %>%
  mutate(term = dplyr::recode(term, !!!coef_labels))

coef_null_sum2 <- coef_null_sum %>%
  filter(term %in% coef_keep) %>%
  mutate(term = dplyr::recode(term, !!!coef_labels))

coef_emp_df2 <- coef_emp_df %>%
  filter(term %in% coef_keep) %>%
  mutate(term = dplyr::recode(term, !!!coef_labels))

term_order <- c(
  "Aquatic dispersal",
  "Terrestrial dispersal",
  "Trait similarity"
)

coef_null_df2$term <- factor(coef_null_df2$term, levels = term_order)
coef_null_sum2$term <- factor(coef_null_sum2$term, levels = term_order)
coef_emp_df2$term <- factor(coef_emp_df2$term, levels = term_order)

coef_null_df2$dummy <- " "
coef_null_sum2$dummy <- " "
coef_emp_df2$dummy <- " "

coef_null_df2 <- coef_null_df2 %>%
  mutate(framework = "Metacommunity")

coef_null_sum2 <- coef_null_sum2 %>%
  mutate(framework = "Metacommunity")

coef_emp_df2 <- coef_emp_df2 %>%
  mutate(framework = "Metacommunity")

figure_2 <- ggplot(coef_null_df2, aes(x = estimate, y = term)) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.8,
    color = "gray40"
  ) +
  ggdist::stat_halfeye(
    aes(fill = framework),
    adjust = 0.7,
    width = 0.65,
    .width = 0,
    justification = -0.15,
    alpha = 0.35,
    color = NA,
    slab_color = NA,
    point_colour = NA
  ) +
  geom_segment(
    data = coef_null_sum2,
    aes(
      x = ci_low,
      xend = ci_high,
      y = term,
      yend = term,
      color = framework
    ),
    linewidth = 1.5,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = coef_null_sum2,
    aes(
      x = mean_est,
      y = term,
      fill = framework
    ),
    shape = 21,
    size = 4,
    color = "black",
    inherit.aes = FALSE
  ) +
  geom_point(
    data = coef_emp_df2,
    aes(
      x = estimate,
      y = term,
      fill = framework
    ),
    shape = 23,
    size = 3,
    color = "black",
    stroke = 1,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = framework_cols) +
  scale_color_manual(values = framework_cols) +
  scale_x_continuous(limits = c(-0.5, 0.95)) +
  labs(
    x = "Coefficient estimate",
    y = NULL
  ) +
  theme_classic(base_size = 15) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())

figure_2

ggsave("Figure_2.jpg", figure_2, width = 6, height = 4)

# Interactions --------------------------------------------------
m_int_1 <- glm(
  S ~ traits_similarity * stream_distance,
  data = env_filt_scaled,
  family = poisson(link = "log")
)

m_int_2 <- glm(
  S ~ traits_similarity * aquatic_dispersal_capacity,
  data = env_filt_scaled,
  family = poisson(link = "log")
)

m_int_3 <- glm(
  S ~ aquatic_dispersal_capacity * stream_distance,
  data = env_filt_scaled,
  family = poisson(link = "log")
)

summary(m_int_1); r2(m_int_1)
summary(m_int_2); r2(m_int_2)
summary(m_int_3); r2(m_int_3)

# Figure 3 -------------------------------------
fig_int_C <- plot_model(
  m_int_1,
  type = "pred",
  terms = c("traits_similarity", "stream_distance")
) +
  labs(
    x = "Trait similarity",
    y = "Predicted richness",
    color = "Stream distance",
    fill = "Stream distance",
    title = NULL
  ) +
  scale_color_brewer(
    palette = "Dark2",
    labels = c("-1" = "Low", "0" = "Mean", "1" = "High")
  ) +
  scale_fill_brewer(
    palette = "Dark2",
    labels = c("-1" = "Low", "0" = "Mean", "1" = "High")
  ) +
  scale_y_continuous(limits = c(0, 17), breaks = seq(0, 17, by = 4)) +
  theme_classic(base_size = 18)+
  theme(
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank()
  )

fig_int_B <- plot_model(
  m_int_3,
  type = "pred",
  terms = c("stream_distance", "aquatic_dispersal_capacity")
) +
  labs(
    color = "Aquatic dispersal",
    fill = "Aquatic dispersal",
    y = "Predicted richness",
    x = "Stream distance",
    title = NULL
  ) +
  scale_color_brewer(
    palette = "Dark2",
    labels = c("-1" = "Low", "0" = "Mean", "1" = "High")
  ) +
  scale_fill_brewer(
    palette = "Dark2",
    labels = c("-1" = "Low", "0" = "Mean", "1" = "High")
  ) +
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, by = 2)) +
  theme_classic(base_size = 18)+
  theme(
    legend.position = c(0.35, 0.99),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank()
  )

fig_int_A <- plot_model(
  m_int_2,
  type = "pred",
  terms = c("traits_similarity", "aquatic_dispersal_capacity")
) +
  labs(
    x = "Trait similarity",
    y = "Predicted richness",
    color = "Aquatic dispersal",
    fill = "Aquatic dispersal",
    title = NULL
  ) +
  scale_color_brewer(
    palette = "Dark2",
    labels = c("-1" = "Low", "0" = "Mean", "1" = "High")
  ) +
  scale_fill_brewer(
    palette = "Dark2",
    labels = c("-1" = "Low", "0" = "Mean", "1" = "High")
  ) +
  scale_y_continuous(limits = c(0, 230), breaks = seq(0, 230, by = 50)) +
  theme_classic(base_size = 18)+
  theme(
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank()
  )

fig_3_new <- fig_int_A + fig_int_B + fig_int_C +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")"
  )


fig_3_new

ggsave("Figure_3.jpg", fig_3_new, width = 13, height = 5)
