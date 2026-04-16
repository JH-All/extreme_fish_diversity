# Packages ---------------------------------------
load_or_install <- function(pkgs) {
  if (!requireNamespace("utils", quietly = TRUE)) stop("utils not available.")
  
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    message("Installing: ", paste(missing_pkgs, collapse = ", "))
    install.packages(missing_pkgs, dependencies = TRUE)
  }
  
  invisible(lapply(pkgs, function(p) {
    suppressPackageStartupMessages(
      library(p, character.only = TRUE)
    )
  }))
  
  message("Packages read.")
  invisible(TRUE)
}


pkgs <- c(
  "readxl", "tidyverse", "vegan", "lme4", "DHARMa", "glmmTMB", "MASS",
  "performance", "cowplot", "car", "mgcv", "broom", "mgcViz", "spdep",
  "RColorBrewer", "ggdist", "purrr", "piecewiseSEM", "ggeffects", 
  "DiagrammeR", "sjPlot", "sjmisc", "ade4", "ape", "patchwork"
)

load_or_install(pkgs)

# Getting data ready -------------------------------------------
data = read_excel("data.xlsx")
assemblages = data[,13:32]
env = data[,c(1,4, 9:12)]

thermal_traits = read_excel("thermal.xlsx")
thermal_traits <- thermal_traits %>%
  mutate(across(2:last_col(), as.numeric))

morphology  = read_excel("morphology.xlsx")
morphology_clean <- morphology %>%
  mutate(across(2:last_col(), as.numeric)) %>%  
  filter(!if_any(where(is.numeric), is.na))

## Calculating and combining traits ------------------
morphology_clean$compression_index = morphology_clean$Bd / morphology_clean$Bw
morphology_clean$relative_depth = morphology_clean$Bd / morphology_clean$SL
morphology_clean$index_ventral_flat = morphology_clean$Mmd / morphology_clean$Bd
morphology_clean$relative_eye_posit = morphology_clean$Eh / morphology_clean$Hd
morphology_clean$finenness = (morphology_clean$SL / sqrt(morphology_clean$Bd)) * morphology_clean$Bw
morphology_clean$relative_mouth_width = morphology_clean$Mw / morphology_clean$SL
morphology_clean$relative_eye_size = morphology_clean$Ed / morphology_clean$Hd


morphology_traits <- morphology_clean %>%
  group_by(species) %>%
  summarise(
    across(
      c(compression_index,
        relative_depth,
        index_ventral_flat,
        relative_eye_posit,
        finenness,
        relative_mouth_width,
        relative_eye_size),
      ~ mean(.x, na.rm = TRUE)
    )
  )

thermal_traits <- thermal_traits %>%
  rename(species = Species)

traits_combined <- morphology_traits %>%
  left_join(thermal_traits, by = "species") %>%
  dplyr::select(-CTMax_CV)

## Species richness ----------------------------------
env$S = specnumber(assemblages)

## Creating some predictors -----------------------------------
### Water volume --------------------
data$maximum_length_m<- as.numeric(data$maximum_length_m)
data$maximum_width_m <- as.numeric(data$maximum_width_m)
data$maximum_depth_cm <- as.numeric(data$maximum_depth_cm)
data$depth_m <- data$maximum_depth_cm / 100
V = (2/3) * pi * data$maximum_length_m * data$maximum_width_m  * data$depth_m
env$Volume = V
env <- env %>%
  mutate(across(3:last_col(), as.numeric))

### Environmental stress indexes ----------------------
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

### Traits similarity ------------------------------
traits_matrix <- traits_combined %>%
  column_to_rownames("species") %>%  
  as.matrix()

traits_scaled <- scale(traits_matrix)
dist_matrix <- dist(traits_scaled, method = "euclidean")
dist_mat <- as.matrix(dist_matrix)

species_in_traits <- rownames(dist_mat)

traits_similarity <- apply(assemblages, 1, function(abund_row) {
  sp_present <- names(abund_row)[abund_row > 0]
  sp_present <- intersect(sp_present, rownames(dist_mat))
  if (length(sp_present) <= 1) return(0)
  
  dsub <- dist_mat[sp_present, sp_present, drop = FALSE]
  dvals <- dsub[upper.tri(dsub)]
  
  mean(1 / (1 + dvals), na.rm = TRUE)
})

env$traits_similarity <- traits_similarity

### Terrestrial dispersal capacity presence  ----------------------
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

### Aquatic dispersal capacity ---------------------
fin_map <- setNames(traits_combined$finenness, traits_combined$species)
species_with_traits <- names(fin_map)

env$aquatic_dispersal_capacity <- apply(assemblages, 1, function(abund_row) {
  sp_present <- names(abund_row)[abund_row > 0]
  
  sp_present <- intersect(sp_present, species_with_traits)
  
  if (length(sp_present) == 0) return(NA_real_)  
  mean(fin_map[sp_present], na.rm = TRUE)
})

### Site average ARR -----------------
arr_map <- setNames(traits_combined$ARR, traits_combined$species)
valid_species <- traits_combined$species

env$ARR <- apply(assemblages, 1, function(abund_row) {
  
  sp_present <- names(abund_row)[abund_row > 0]
  sp_present <- intersect(sp_present, valid_species)
  if (length(sp_present) == 0) return(NA_real_)
  
  mean(arr_map[sp_present], na.rm = TRUE)
})

### Scale all predictors ----------------------
env_scaled <- env %>%
  mutate(
    across(
      where(is.numeric) & !any_of(c("S", "terrestrial_dispersal_capacity")),
      ~ scale(.x)[,1]
    )
  )

# Separated for models with trait similarity
env_filt <- env %>%
  filter(S > 1)

env_filt_scaled <- env_filt %>%
  mutate(
    across(
      where(is.numeric) & !any_of(c("S", "terrestrial_dispersal_capacity")),
      ~ scale(.x)[,1]
    )
  )

# Model 1 -----------------------------------------
m1 <- glm.nb(S ~ stream_distance + Volume + Period, data = env_scaled)

summary(m1)
r2(m1)

# Model 2 ---------------------------------------
m2 <- glm.nb(S ~ Temp + DO + pH + Volume + Period, data = env_scaled)

summary(m2)
r2(m2)

# Model 3 ---------------------------------------
m3 <- glm.nb(S ~ Stress_index_1 + Stress_index_2 + Period, data = env_scaled)

summary(m3)
r2(m3)

# Model 4 ----------------------------------------
m4 <- glm(
  S ~ traits_similarity + aquatic_dispersal_capacity +
    terrestrial_dispersal_capacity + Period,
  data = env_filt_scaled,
  family = poisson(link = "log")
)

summary(m4)
r2(m4)

# Model 5 --------------------------------------------
m5 <- glm.nb(S ~ ARR + Period, data = env_scaled)
summary(m5)
r2(m5)

# Figure 1A ------------------------------------------
model_list <- list(
  "Island biogeography" = m1,
  "Environment" = m2,
  "Environmental filtering" = m3,
  "Metacommunity" = m4,
  "Acclimation response" = m5
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
      Temp = "Temperature",
      DO = "Dissolved oxygen",
      Stress_index_1 = "Stress index 1",
      Stress_index_2 = "Stress index 2",
      traits_similarity = "Trait similarity",
      aquatic_dispersal_capacity = "Aquatic dispersal",
      terrestrial_dispersal_capacity = "Terrestrial dispersal",
      ARR = "Acclimation response ratio"
    )
  )


coef_df <- coef_df %>%
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
  mutate(
    xmin = -Inf,
    xmax = Inf
  )

label_df <- coef_df %>%
  distinct(term_id, term)

pal <- brewer.pal(5, "Dark2")
names(pal) <- levels(coef_df$model)

fig1_A = ggplot(coef_df, aes(x = estimate, y = term_id, color = model)) +
  geom_rect(
    data = bg_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = model),
    inherit.aes = FALSE,
    alpha = 0.18,
    color = NA
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.8) +
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high),
    width = 0.18,
    linewidth = 1.6,
    color = "black"
  ) +
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high, color = model),
    width = 0.12,
    linewidth = 1.0
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
  scale_x_continuous(limits = c(-1.2, 1.2))+
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
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
  )+
  theme(
    legend.title = element_text(face = "bold", size = 16),
    legend.text  = element_text(face = "bold", size = 14),
    
    legend.position = "right",
    legend.justification = "top",
    legend.box.just = "top"
  )

fig1_A

# Figure 1B -------------------------------
model_list <- list(
  "Island biogeography" = m1,
  "Environment" = m2,
  "Environmental filtering" = m3,
  "Metacommunity" = m4,
  "Acclimation response" = m5
)

pal <- brewer.pal(5, "Dark2")
names(pal) <- names(model_list)


r2_df <- purrr::imap_dfr(model_list, ~{
  tibble(
    model = .y,
    R2 = unname(performance::r2(.x)$R2)
  )
})


r2_df <- r2_df %>%
  arrange(desc(R2)) %>%
  mutate(
    model = factor(model, levels = model)
  )

fig1_B = ggplot(r2_df, aes(x = model, y = R2)) +
  geom_col(
    color = "black",
    width = 0.93,
    linewidth = 1,
    alpha = 0.85,
    fill = "lightgray"
  ) +
  coord_flip() +
  theme_classic(base_size = 20) +
  labs(
    x = NULL,
    y = expression(R^2)
  ) +
  scale_y_continuous(
    expand = c(0,0),
    limits = c(0, 0.50),
    breaks = seq(0, 0.50, by = .1)
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(face = "bold"),
    legend.title = element_text(face = "bold", size = 16),
    legend.text = element_text(face = "bold", size = 14),
    legend.position = "none"
  )

fig1_B

# Figure 1C ---------------------------
aic_df <- purrr::imap_dfr(model_list, ~{
  tibble(
    model = .y,
    AIC = AIC(.x)
  )
})

aic_df <- aic_df %>%
  mutate(delta_AIC = AIC - min(AIC)) %>%
  arrange(delta_AIC) %>%
  mutate(model = factor(model, levels = model))

fig1_C = ggplot(aic_df, aes(x = model, y = delta_AIC)) +
  geom_col(
    color = "black",
    width = 0.93,
    linewidth = 1,
    show.legend = FALSE,
    fill = "lightgray"
  ) +
  coord_flip() +
  theme_classic(base_size = 20) +
  labs(
    x = NULL,
    y = expression(Delta*AIC)
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(face = "bold"),
    legend.position = "none"
  ) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0,165),
                     breaks = seq(0,165, by = 30))

fig1_C

# Figure Complete ---------------------------
fig1 = fig1_A + (fig1_B / fig1_C) +
  plot_annotation(tag_levels = "A")

fig1

ggsave("Figure_1.jpg", fig1, width = 15, height = 6.5)

# Spatial autocorrelation -----------------------------
env$coord <- data$coordinates

dms_to_decimal <- function(x){
  x <- str_trim(x)
  
  parts <- str_split_fixed(x, "\\s+", 2)
  lat_s <- parts[,1]
  lon_s <- parts[,2]
  
  parse_one <- function(s){
    s <- str_replace_all(s, "\\s+", "")
    
    m <- str_match(s, "^(\\d+)°(\\d+)'([0-9.]+)\"?([NSEW])$")
    if(any(is.na(m))) {
      stop("Falha ao parsear coordenada: ", s[which(is.na(m))[1]])
    }
    
    deg <- as.numeric(m[,2])
    min <- as.numeric(m[,3])
    sec <- as.numeric(m[,4])
    hem <- m[,5]
    
    dec <- deg + min/60 + sec/3600
    dec[hem %in% c("S","W")] <- -dec[hem %in% c("S","W")]
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
  bind_cols(coords_df %>% dplyr::select(lon, lat))

env_spatial <- env_spatial %>%
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
      ~ scale(.x)[,1]
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
      ~ scale(.x)[,1]
    )
  )

m1_sp <- glm.nb(S ~ stream_distance + Volume, data = env_spatial_scaled)

m2_sp <- glm.nb(S ~ Temp + DO + pH + Volume, data = env_spatial_scaled)

m3_sp <- glm.nb(S ~ Stress_index_1 + Stress_index_2, data = env_spatial_scaled)

m4_sp <- glm(
  S ~ traits_similarity + aquatic_dispersal_capacity +
    terrestrial_dispersal_capacity,
  data = env_spatial_filt_scaled,
  family = poisson(link = "log")
)

m5_sp <- glm.nb(S ~ ARR, data = env_spatial_scaled)

moran_residuals <- function(mod, dat, k = 4){
  
  dat2 <- dat %>%
    filter(is.finite(lon_jit), is.finite(lat_jit))
  
  r <- residuals(mod, type = "pearson")
  
  if(length(r) != nrow(dat2)){
    stop("Mismatch entre resíduos e coordenadas: resíduos = ", length(r),
         ", coordenadas = ", nrow(dat2))
  }
  
  if(nrow(dat2) < (k + 2)){
    return(tibble(
      Moran_I = NA_real_,
      p_value = NA_real_,
      n = nrow(dat2)
    ))
  }
  
  coords <- as.matrix(dat2[, c("lon_jit", "lat_jit")])
  
  nb <- knn2nb(knearneigh(coords, k = k))
  lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
  
  mi <- moran.test(r, lw, zero.policy = TRUE)
  
  tibble(
    Moran_I = unname(mi$estimate[["Moran I statistic"]]),
    p_value = mi$p.value,
    n = length(r)
  )
}

model_list_spatial <- list(
  "Island biogeography" = list(model = m1_sp, data = env_spatial_scaled),
  "Environment" = list(model = m2_sp, data = env_spatial_scaled),
  "Environmental stress" = list(model = m3_sp, data = env_spatial_scaled),
  "Metacommunity" = list(model = m4_sp, data = env_spatial_filt_scaled),
  "Acclimation Response" = list(model = m5_sp, data = env_spatial_scaled)
)

spatial_tests <- purrr::imap_dfr(model_list_spatial, function(x, model_name){
  
  moran_residuals(
    mod = x$model,
    dat = x$data,
    k = 4
  ) %>%
    mutate(Model = model_name, .before = 1)
})

spatial_tests %>%
  arrange(p_value)

# Figure S1A ----------------------------------------
## Redoing models with site richness > 1
m1_s <- glm.nb(S ~ stream_distance + Volume,
             data = env_filt_scaled)
m2_s <- glm.nb(S ~ Temp + DO + pH + Volume,
               data = env_filt_scaled)
m3_s <- glm.nb(S ~ Stress_index_1 + Stress_index_2, 
               data = env_filt_scaled)
m5_s <- glm.nb(S ~ ARR, data = env_filt_scaled)

model_list_s <- list(
  "Island biogeography" = m1_s,
  "Environment" = m2_s,
  "Environmental filtering" = m3_s,
  "Metacommunity" = m4,
  "Acclimation response" = m5_s
)

coef_df_s <- purrr::imap_dfr(model_list_s, ~{
  broom::tidy(.x, conf.int = TRUE) %>%
    filter(term != "(Intercept)",
           term != "PeriodWet Period") %>%
    mutate(model = .y)
})

coef_df_s$model <- factor(coef_df_s$model, levels = names(model_list_s))

coef_df_s <- coef_df_s %>%
  mutate(
    term = dplyr::recode(
      term,
      stream_distance = "Stream distance",
      Temp = "Temperature",
      DO = "Dissolved oxygen",
      Stress_index_1 = "Stress index 1",
      Stress_index_2 = "Stress index 2",
      traits_similarity = "Trait similarity",
      aquatic_dispersal_capacity = "Aquatic dispersal",
      terrestrial_dispersal_capacity = "Terrestrial dispersal",
      ARR = "Acclimation response ratio"
    )
  )


coef_df_s <- coef_df_s %>%
  mutate(term_id = paste(model, term, sep = "___"))

term_order_s <- coef_df_s %>%
  mutate(model_num = as.numeric(model)) %>%
  arrange(model_num) %>%
  pull(term_id)

coef_df_s$term_id <- factor(coef_df_s$term_id, levels = rev(term_order))

bg_df_s <- coef_df_s %>%
  distinct(model, term_id) %>%
  mutate(y = as.numeric(term_id)) %>%
  group_by(model) %>%
  summarise(
    ymin = min(y) - 0.5,
    ymax = max(y) + 0.5,
    .groups = "drop"
  ) %>%
  mutate(
    xmin = -Inf,
    xmax = Inf
  )

label_df_s <- coef_df_s %>%
  distinct(term_id, term)

pal <- brewer.pal(5, "Dark2")
names(pal) <- levels(coef_df_s$model)

fig_s1_A = ggplot(coef_df_s, aes(x = estimate, y = term_id, color = model)) +
  geom_rect(
    data = bg_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = model),
    inherit.aes = FALSE,
    alpha = 0.18,
    color = NA
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.8) +
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high),
    width = 0.18,
    linewidth = 1.6,
    color = "black"
  ) +
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high, color = model),
    width = 0.12,
    linewidth = 1.0
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
  scale_x_continuous(limits = c(-1.2, 1.2))+
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
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
  )+
  theme(
    legend.title = element_text(face = "bold", size = 16),
    legend.text  = element_text(face = "bold", size = 14),
    
    legend.position = "right",
    legend.justification = "top",
    legend.box.just = "top"
  )

fig_s1_A

# Figure S1B ----------------------------------------
model_list_s <- list(
  "Island biogeography" = m1_s,
  "Environment" = m2_s,
  "Environmental filtering" = m3_s,
  "Metacommunity" = m4,
  "Acclimation response" = m5_s
)

pal <- brewer.pal(5, "Dark2")
names(pal) <- names(model_list_s)

r2_df_s <- purrr::imap_dfr(model_list_s, ~{
  tibble(
    model = .y,
    R2 = unname(performance::r2(.x)$R2)
  )
})


r2_df_s <- r2_df_s %>%
  arrange(desc(R2)) %>%
  mutate(
    model = factor(model, levels = model)
  )

fig_s1_B = ggplot(r2_df_s, aes(x = model, y = R2)) +
  geom_col(
    color = "black",
    width = 0.93,
    linewidth = 1,
    alpha = 0.85,
    fill = "lightgray"
  ) +
  coord_flip() +
  theme_classic(base_size = 20) +
  labs(
    x = NULL,
    y = expression(R^2)
  ) +
  scale_y_continuous(
    expand = c(0,0),
    limits = c(0, 0.50),
    breaks = seq(0, 0.50, by = .1)
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(face = "bold"),
    legend.title = element_text(face = "bold", size = 16),
    legend.text = element_text(face = "bold", size = 14),
    legend.position = "none"
  )

fig_s1_B

# Figure S1C -----------------------------------------
aic_df_s <- purrr::imap_dfr(model_list_s, ~{
  tibble(
    model = .y,
    AIC = AIC(.x)
  )
})

aic_df_s <- aic_df_s %>%
  mutate(delta_AIC = AIC - min(AIC)) %>%
  arrange(delta_AIC) %>%
  mutate(model = factor(model, levels = model))

fig_s1_C = ggplot(aic_df_s, aes(x = model, y = delta_AIC)) +
  geom_col(
    color = "black",
    width = 0.93,
    linewidth = 1,
    show.legend = FALSE,
    fill = "lightgray"
  ) +
  coord_flip() +
  theme_classic(base_size = 20) +
  labs(
    x = NULL,
    y = expression(Delta*AIC)
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(face = "bold"),
    legend.position = "none"
  ) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0,25),
                     breaks = seq(0,25, by = 5))

fig_s1_C

# Figure S1 Complete ---------------------------
fig_s1 = fig_s1_A + (fig_s1_B / fig_s1_C) +
  plot_annotation(tag_levels = "A")

fig_s1

ggsave("Figure_S1.jpg", fig_s1, width = 15, height = 6.5)

# Null models -----------------------------------
## Functions ----------------------
only_scale <- function(x){
  scale(x)[,1]
}

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

calc_mean_from_map <- function(sp_present, map_obj, valid_species) {
  sp_present <- intersect(sp_present, valid_species)
  if (length(sp_present) == 0) return(NA_real_)
  mean(map_obj[sp_present], na.rm = TRUE)
}

## Species pool ----------------------
species_pool <- colnames(assemblages)
site_richness <- env$S
n_reps <- 1000

## Null assemblages ------------------
generate_null_assemblage <- function(assemblages, site_richness, species_pool) {
  
  null_mat <- matrix(
    0,
    nrow = nrow(assemblages),
    ncol = length(species_pool),
    dimnames = list(rownames(assemblages), species_pool)
  )
  
  for(i in seq_len(nrow(assemblages))) {
    S_i <- site_richness[i]
    if (is.na(S_i) || S_i == 0) next
    sp_i <- sample(species_pool, size = S_i, replace = FALSE)
    null_mat[i, sp_i] <- 1
  }
  
  as.data.frame(null_mat)
}

## Build null predictors dataframe ------------------
build_null_env <- function(null_assemblage, env) {
  
  sp_list_by_site_null <- apply(null_assemblage, 1, function(abund_row) {
    sp_present <- names(abund_row)[abund_row > 0]
    sp_present <- intersect(sp_present, species_pool)
    sp_present
  })
  
  env_null <- env %>%
    dplyr::select(S)
  
  env_null$traits_similarity <- vapply(
    sp_list_by_site_null, calc_traits_similarity_sp, numeric(1)
  )
  
  env_null$aquatic_dispersal_capacity <- vapply(
    sp_list_by_site_null, calc_aquatic_disp, numeric(1)
  )
  
  env_null$terrestrial_dispersal_capacity <- vapply(
    sp_list_by_site_null, calc_terrestrial_disp, numeric(1)
  )
  
  env_null$ARR <- vapply(
    sp_list_by_site_null, calc_mean_from_map, numeric(1),
    map_obj = arr_map, valid_species = valid_species
  )
  
  env_null
}

## Scaled versions -----------------------
# for m5
prepare_env_null_scaled <- function(env_null) {
  env_null %>%
    mutate(
      across(
        where(is.numeric) & !any_of(c("S", "terrestrial_dispersal_capacity")),
        ~ only_scale(.x)
      )
    )
}

# for m4:  S > 1 
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

## Fit null models -----------------------------
fit_null_models <- function(env_null_scaled, env_null_scaled_filt) {
  
  null_metacom <- glm(
    S ~ traits_similarity + aquatic_dispersal_capacity +
      terrestrial_dispersal_capacity,
    data = env_null_scaled_filt,
    family = poisson(link = "log")
  )
  
  null_acclimation <- glm.nb(
    S ~ ARR,
    data = env_null_scaled
  )
  
  list(
    "Metacommunity" = null_metacom,
    "Acclimation response" = null_acclimation
  )
}

## Run null models --------------------------------
set.seed(123)
null_results <- vector("list", n_reps)

for(r in 1:n_reps) {
  
  null_assemblage <- generate_null_assemblage(
    assemblages = assemblages,
    site_richness = site_richness,
    species_pool = species_pool
  )
  
  env_null <- build_null_env(null_assemblage, env)
  env_null_scaled <- prepare_env_null_scaled(env_null)
  env_null_scaled_filt <- prepare_env_null_scaled_filt(env_null)
  
  null_results[[r]] <- suppressWarnings(
    fit_null_models(env_null_scaled, env_null_scaled_filt)
  )
}

length(null_results)
names(null_results[[1]])

## Extract null coefficients ----------------
coef_null_df <- purrr::imap_dfr(null_results, ~{
  
  rep_id <- .y
  model_list_rep <- .x
  
  purrr::imap_dfr(model_list_rep, ~{
    broom::tidy(.x, conf.int = TRUE) %>%
      dplyr::filter(term != "(Intercept)") %>%
      dplyr::mutate(
        replicate = rep_id,
        model = .y
      )
  })
})

coef_null_sum <- coef_null_df %>%
  group_by(model, term) %>%
  summarise(
    mean_est = mean(estimate, na.rm = TRUE),
    ci_low = quantile(estimate, 0.025, na.rm = TRUE),
    ci_high = quantile(estimate, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

## Empirical models values -----------------------------
emp_models <- list(
  "Metacommunity" = m4,
  "Acclimation response" = m5
)

coef_emp_df <- purrr::imap_dfr(emp_models, ~{
  broom::tidy(.x) %>%
    dplyr::filter(term != "(Intercept)") %>%
    dplyr::mutate(model = .y)
})

## Figure 2 ------------------------
coef_keep <- c(
  "traits_similarity",
  "aquatic_dispersal_capacity",
  "terrestrial_dispersal_capacity",
  "ARR"
)

coef_labels <- c(
  traits_similarity = "Trait similarity",
  aquatic_dispersal_capacity = "Aquatic dispersal",
  terrestrial_dispersal_capacity = "Terrestrial dispersal",
  ARR = "Acclimation response ratio"
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
  "Trait similarity",
  "Aquatic dispersal",
  "Terrestrial dispersal",
  "Acclimation response ratio"
)

coef_null_df2$term  <- factor(coef_null_df2$term, levels = term_order)
coef_null_sum2$term <- factor(coef_null_sum2$term, levels = term_order)
coef_emp_df2$term   <- factor(coef_emp_df2$term, levels = term_order)

coef_null_df2$dummy <- " "
coef_null_sum2$dummy <- " "
coef_emp_df2$dummy  <- " "

figure_2 = ggplot(coef_null_df2, aes(x = estimate, y = dummy)) +
  geom_errorbarh(
    data = coef_null_sum2,
    aes(xmin = ci_low, xmax = ci_high, y = dummy),
    width = 0.05,
    linewidth = 2.4,
    color = "black",
    inherit.aes = FALSE
  )+
  geom_errorbarh(
    data = coef_null_sum2,
    aes(xmin = ci_low, xmax = ci_high, y = dummy),
    width = 0.03,
    linewidth = 1.7,
    color = "gray40",
    inherit.aes = FALSE
  ) +
  
  geom_point(
    data = coef_null_sum2,
    aes(x = mean_est, y = dummy),
    color = "grey20",
    fill = "gray40",
    shape = 21,
    size = 4.5,
    stroke = 1.1,
    inherit.aes = FALSE
  ) +
  
  geom_vline(
    data = coef_emp_df2,
    aes(xintercept = estimate),
    color = "darkred",
    linetype = "dashed",
    linewidth = 1,
    alpha = 0.8
  ) +
  ggdist::stat_halfeye(
    adjust = 0.7,
    width = 0.6,
    .width = 0,
    justification = -0.2,
    alpha = 0.8,
    point_colour = NA,
    fill = "gray70",
    color = "black"
  ) +
  facet_wrap(~ term, ncol = 2, scales = "fixed") +
  scale_y_discrete(expand = expansion(mult = c(0.10, 0.35))) +
  scale_x_continuous(limits = c(-0.95, 0.95)) +
  theme_classic(base_size = 18) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  labs(
    x = "Coefficient estimate",
    y = NULL
  )

figure_2

ggsave("Figure_2.jpg", figure_2, width = 10, height = 7)

# Interactions -----------------------------------------
## Model Int 1 -----------------------------
m_int_1 <- glm(
  S ~ traits_similarity * stream_distance,
  data = env_filt_scaled,
  family = poisson(link = "log")
)

summary(m_int_1)
r2(m_int_1)

## Model Int 2 -----------------------------
m_int_2 <- glm(
  S ~ traits_similarity * aquatic_dispersal_capacity,
  data = env_filt_scaled,
  family = poisson(link = "log")
)

summary(m_int_2)
r2(m_int_2)

## Model Int 3 -----------------------------
m_int_3 <- glm(
  S ~ aquatic_dispersal_capacity * stream_distance,
  data = env_filt_scaled,
  family = poisson(link = "log")
)

summary(m_int_3)
r2(m_int_3)

## Figure 3A ----------------------
int_A = plot_model(m_int_1, type = "pred",
                   terms = c("traits_similarity",  "stream_distance"))
int_A

fig_int_A = int_A +
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
  scale_y_continuous(
    limits = c(1, 12),
    breaks = seq(2, 12, by = 2)
  ) +
  theme_classic(base_size = 18) +
  theme(
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank()
  )

fig_int_A

## Figure 3B  ----------------------
int_B = plot_model(m_int_3, type = "pred",
                   terms = c("stream_distance",
                             "aquatic_dispersal_capacity"))
int_B

fig_int_B = int_B +
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
  scale_y_continuous(
    limits = c(0, 8),
    breaks = seq(0, 8, by = 2)
  ) +
  theme_classic(base_size = 18) +
  theme(
    legend.position = c(0.35, 0.99),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank()
  )

fig_int_B

## Figure 3C ----------------------
int_C = plot_model(m_int_2, type = "pred",
                   terms = c("traits_similarity",  "aquatic_dispersal_capacity"))
int_C

fig_int_C = int_C +
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
  scale_y_continuous(
    limits = c(0, 40),
    breaks = seq(0, 40, by = 10)
  ) +
  theme_classic(base_size = 18) +
  theme(
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank()
  )

fig_int_C

## Figure 3 Complete ----------------------------------

fig_3 = fig_int_A + fig_int_B + fig_int_C

fig_3_new = fig_3 + 
  plot_annotation(tag_levels = 'A')

fig_3_new

ggsave("Figure_3.jpg", fig_3_new, width = 13, height = 5)