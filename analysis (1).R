# Crop yield, climate and management in Italy (2010-2025)
# Author: Prabhath Samarahewa
# October 2026
#
# What this script does:
#   - builds a small farm database (SQLite) from simulated data
#   - joins the tables and prepares the data
#   - fits a mixed model for crop yield
#   - compares a few climate and adaptation scenarios
#
# The data are simulated, not real observations. Needs R >= 4.1.


# 1. Setup ---------------------------------------------------------------

pkgs <- c("tidyverse", "DBI", "RSQLite", "lme4", "lmerTest",
          "broom", "broom.mixed", "patchwork", "scales")
missing_pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing_pkgs)) install.packages(missing_pkgs)
invisible(lapply(pkgs, library, character.only = TRUE))

set.seed(2026)
dir.create("data", showWarnings = FALSE)
dir.create("outputs/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/tables", showWarnings = FALSE, recursive = TRUE)

# assumptions used in the simulation and in the analysis
years          <- 2010:2025
TEMP_CENTRE    <- 20    # deg C, temperature is centred on this value in the model
IRRIG_WATER_MM <- 150   # extra water from irrigation (mm)
COST_IRRIG     <- 120   # EUR/ha
COST_OTHER     <- 550   # EUR/ha, other variable costs


# 2. Simulate the data and build the database ----------------------------

# regions: crop-season climate, soil quality and share of irrigated farms
regions <- tibble(
  region_id   = 1:6,
  region      = c("Lombardia", "Veneto", "Emilia-Romagna",
                  "Toscana", "Puglia", "Sicilia"),
  country     = "IT",
  area        = c("North", "North", "North", "Centre", "South", "South"),
  base_temp   = c(19.0, 19.5, 20.0, 20.5, 22.0, 22.5),   # deg C
  base_precip = c(480, 500, 330, 280, 200, 190),         # mm
  soil_mean   = c(0.78, 0.75, 0.75, 0.60, 0.50, 0.48),   # index 0-1
  irrig_prob  = c(0.60, 0.45, 0.55, 0.20, 0.50, 0.45)
)

# climate per region and year: slow warming, slight drying, and
# year-to-year anomalies shared by all regions
shared_anom <- tibble(year = years,
                      t_anom = rnorm(length(years), 0, 0.5),
                      p_anom = rnorm(length(years), 0, 40))

climate <- expand_grid(region_id = regions$region_id, year = years) |>
  left_join(regions, by = "region_id") |>
  left_join(shared_anom, by = "year") |>
  mutate(
    temp_gs   = base_temp + 0.05 * (year - 2010) + t_anom + rnorm(n(), 0, 0.4),
    precip_gs = pmax(100, base_precip - 1.5 * (year - 2010) + p_anom + rnorm(n(), 0, 40))
  ) |>
  select(region_id, year, temp_gs, precip_gs)

# farms (farm_re is a hidden farm effect, it is not stored in the database)
farms_truth <- tibble(
  farm_id   = sprintf("F%03d", 1:150),
  region_id = sample(regions$region_id, 150, replace = TRUE,
                     prob = c(.20, .18, .20, .12, .18, .12))
) |>
  left_join(select(regions, region_id, soil_mean, irrig_prob), by = "region_id") |>
  mutate(
    area_ha      = round(rlnorm(n(), log(40), 0.8), 1),
    soil_quality = pmin(pmax(rnorm(n(), soil_mean, 0.10), 0.10), 1),
    irrigated    = rbinom(n(), 1, irrig_prob),
    practice     = sample(c("conventional", "conservation", "organic"),
                          n(), replace = TRUE, prob = c(.60, .25, .15)),
    farm_re      = rnorm(n(), 0, 0.08)
  )

farms_db <- select(farms_truth, farm_id, region_id, area_ha,
                   soil_quality, irrigated, practice)

# crops, prices and nitrogen price
crop_params <- tibble(
  crop       = c("wheat", "maize", "sunflower"),
  base_yield = c(4.8, 7.0, 2.6),      # t/ha at reference conditions
  t_opt      = c(18.5, 21.0, 21.0),   # best crop-season temperature (deg C)
  n_mean     = c(150, 200, 60),       # usual nitrogen rate (kg/ha)
  base_price = c(190, 180, 380)       # EUR/t
)

prices <- expand_grid(crop = crop_params$crop, year = years) |>
  left_join(select(crop_params, crop, base_price), by = "crop") |>
  mutate(price_eur_t = round(base_price * (1 + 0.02 * (year - 2010)) *
                               exp(rnorm(n(), 0, 0.08)), 1)) |>
  select(crop, year, price_eur_t)

input_costs <- tibble(year = years,
                      n_cost_eur_kg = round(0.9 * exp(rnorm(length(years), 0, 0.20)), 3))

# field observations: farm x year panel with 8% of rows missing.
# Yield is generated on the log scale from the terms below.
field_obs <- expand_grid(farm_id = farms_truth$farm_id, year = years) |>
  left_join(farms_truth, by = "farm_id") |>
  left_join(climate, by = c("region_id", "year")) |>
  slice_sample(prop = 0.92) |>
  mutate(crop = sample(crop_params$crop, n(), replace = TRUE, prob = c(.5, .3, .2))) |>
  left_join(crop_params, by = "crop") |>
  mutate(
    n_mult   = case_when(practice == "organic"      ~ 0.35,
                         practice == "conservation" ~ 0.90,
                         TRUE                       ~ 1.00),
    n_kg_ha  = round(pmax(5, rnorm(n(), n_mean * n_mult, 25)), 1),
    water_mm = precip_gs + IRRIG_WATER_MM * irrigated,
    log_y =  log(base_yield) +
      0.30 * (n_kg_ha / 100) - 0.06 * (n_kg_ha / 100)^2 +            # nitrogen
      -0.012 * (temp_gs - t_opt)^2 +                                 # temperature
      0.0009 * (water_mm - 400) - 0.0000015 * (water_mm - 400)^2 +   # water
      0.30 * (soil_quality - 0.6) +                                  # soil
      0.006 * (year - 2010) +                                        # time trend
      case_when(practice == "organic" ~ -0.12,
                practice == "conservation" ~ -0.03, TRUE ~ 0) +
      farm_re + rnorm(n(), 0, 0.10),
    yield_t_ha = round(exp(log_y), 2)
  ) |>
  select(farm_id, year, crop, n_kg_ha, yield_t_ha)

# write the tables to SQLite
con <- dbConnect(RSQLite::SQLite(), "data/agri_synthetic.sqlite")
dbWriteTable(con, "regions",     select(regions, region_id, region, country, area), overwrite = TRUE)
dbWriteTable(con, "climate",     climate,     overwrite = TRUE)
dbWriteTable(con, "farms",       farms_db,    overwrite = TRUE)
dbWriteTable(con, "field_obs",   field_obs,   overwrite = TRUE)
dbWriteTable(con, "prices",      prices,      overwrite = TRUE)
dbWriteTable(con, "input_costs", input_costs, overwrite = TRUE)
print(dbListTables(con))


# 3. Query and prepare the data ------------------------------------------

query <- "
SELECT o.farm_id, o.year, o.crop, o.n_kg_ha, o.yield_t_ha,
       f.region_id, r.region, r.country,
       f.area_ha, f.soil_quality, f.irrigated, f.practice,
       c.temp_gs, c.precip_gs,
       p.price_eur_t, ic.n_cost_eur_kg
FROM field_obs o
JOIN farms       f  ON f.farm_id   = o.farm_id
JOIN regions     r  ON r.region_id = f.region_id
JOIN climate     c  ON c.region_id = f.region_id AND c.year = o.year
JOIN prices      p  ON p.crop      = o.crop      AND p.year = o.year
JOIN input_costs ic ON ic.year     = o.year
"
raw <- as_tibble(dbGetQuery(con, query))
dbDisconnect(con)

dat <- raw |>
  mutate(
    crop      = factor(crop, levels = c("wheat", "maize", "sunflower")),
    practice  = factor(practice, levels = c("conventional", "conservation", "organic")),
    n_100     = n_kg_ha / 100,                     # nitrogen in 100 kg/ha
    temp_c    = temp_gs - TEMP_CENTRE,
    water_c   = (precip_gs + IRRIG_WATER_MM * irrigated - 400) / 100,   # per 100 mm
    trend     = year - 2010,
    log_yield = log(yield_t_ha),
    rel_yield = yield_t_ha / ave(yield_t_ha, crop),  # yield relative to the crop mean
    pfp       = yield_t_ha * 1000 / n_kg_ha,         # kg grain per kg N
    gm_eur_ha = yield_t_ha * price_eur_t - n_kg_ha * n_cost_eur_kg -
                COST_IRRIG * irrigated - COST_OTHER
  )

stopifnot(!anyNA(dat))
cat(nrow(dat), "observations from", n_distinct(dat$farm_id), "farms\n")


# 4. Summary tables ------------------------------------------------------

tab_region_crop <- dat |>
  group_by(region, crop) |>
  summarise(n_obs = n(), mean_yield = mean(yield_t_ha), sd_yield = sd(yield_t_ha),
            mean_N = mean(n_kg_ha), .groups = "drop")

tab_practice <- dat |>
  group_by(practice) |>
  summarise(n_obs = n(), rel_yield = mean(rel_yield), mean_N = mean(n_kg_ha),
            PFP_kg_per_kgN = mean(pfp), gross_margin_eur_ha = mean(gm_eur_ha),
            .groups = "drop")

# yield trend per region (% per year), adjusted for crop
tab_trend <- dat |>
  group_by(region) |>
  group_modify(~ tidy(lm(log_yield ~ year + crop, data = .x), conf.int = TRUE) |>
                 filter(term == "year")) |>
  ungroup() |>
  transmute(region,
            trend_pct_per_year = 100 * (exp(estimate) - 1),
            ci_low  = 100 * (exp(conf.low)  - 1),
            ci_high = 100 * (exp(conf.high) - 1),
            p_value = p.value)

print(tab_practice)
print(tab_trend)


# 5. Yield model ---------------------------------------------------------

# log yield with a crop-specific temperature effect, curved water and
# nitrogen responses, soil, practice and a time trend.
# Random intercepts for region and for farm within region.
m_yield <- lmer(
  log_yield ~ crop * temp_c + I(temp_c^2) +
    water_c + I(water_c^2) +
    n_100 + I(n_100^2) +
    soil_quality + practice + trend +
    (1 | region / farm_id),
  data = dat
)
summary(m_yield)

# back-transform from log scale (with the usual variance correction)
back_transform <- function(m, eta) exp(eta + sigma(m)^2 / 2)

# marginal and conditional R2
vf <- var(as.vector(model.matrix(m_yield) %*% fixef(m_yield)))
vc <- as.data.frame(VarCorr(m_yield))
v_rand  <- sum(vc$vcov[vc$grp != "Residual"])
v_resid <- vc$vcov[vc$grp == "Residual"]
r2_marg <- vf / (vf + v_rand + v_resid)
r2_cond <- (vf + v_rand) / (vf + v_rand + v_resid)
cat(sprintf("R2 marginal %.2f, conditional %.2f\n", r2_marg, r2_cond))

tab_coef <- tidy(m_yield, effects = "fixed", conf.int = TRUE)

# nitrogen rate that maximises yield according to the model
fe <- fixef(m_yield)
n_opt <- -fe["n_100"] / (2 * fe["I(n_100^2)"]) * 100
cat(sprintf("Yield-maximising N rate: %.0f kg/ha\n", n_opt))

# compare estimates with the values used in the simulation
truth <- c(n_100 = 0.30, `I(n_100^2)` = -0.06, soil_quality = 0.30,
           practiceconservation = -0.03, practiceorganic = -0.12, trend = 0.006)
tab_recovery <- tibble(term = names(truth), true = unname(truth),
                       estimated = unname(fe[names(truth)])) |>
  mutate(abs_error = abs(estimated - true))
print(tab_recovery)

# cross-validation, leaving out whole farms (5 folds)
farm_folds <- tibble(farm_id = unique(dat$farm_id)) |>
  mutate(fold = sample(rep(1:5, length.out = n())))
dat_cv <- left_join(dat, farm_folds, by = "farm_id")

tab_cv <- map_dfr(1:5, function(k) {
  train <- filter(dat_cv, fold != k)
  test  <- filter(dat_cv, fold == k)
  m_k   <- update(m_yield, data = train)
  pred  <- back_transform(m_k, predict(m_k, newdata = test, re.form = NA))
  tibble(fold = k,
         RMSE_t_ha = sqrt(mean((test$yield_t_ha - pred)^2)),
         R2        = cor(test$yield_t_ha, pred)^2)
})
print(tab_cv)


# 6. Scenarios -----------------------------------------------------------

# the adaptation scenario irrigates every field and caps nitrogen at 180 kg/ha
scenarios <- tribble(
  ~scenario,                          ~dT,  ~dP_pct, ~irrigate_all, ~n_cap,
  "Baseline",                          0.0,   0,      FALSE,         Inf,
  "+1.5 C",                            1.5,   0,      FALSE,         Inf,
  "+3 C",                              3.0,   0,      FALSE,         Inf,
  "+3 C, -10% rain",                   3.0, -10,      FALSE,         Inf,
  "+3 C, -10% rain + adaptation",      3.0, -10,      TRUE,          180
)

make_scenario <- function(data, dT, dP_pct, irrigate_all, n_cap) {
  data |>
    mutate(
      irrigated_s = if (irrigate_all) 1L else irrigated,
      n_kg_ha_s   = pmin(n_kg_ha, n_cap),
      temp_c      = (temp_gs + dT) - TEMP_CENTRE,
      water_c     = (precip_gs * (1 + dP_pct / 100) +
                       IRRIG_WATER_MM * irrigated_s - 400) / 100,
      n_100       = n_kg_ha_s / 100
    )
}

scen_results <- pmap_dfr(scenarios, function(scenario, dT, dP_pct, irrigate_all, n_cap) {
  nd <- make_scenario(dat, dT, dP_pct, irrigate_all, n_cap)
  nd$yield_pred <- back_transform(m_yield, predict(m_yield, newdata = nd, re.form = NA))
  nd |>
    mutate(scenario = scenario,
           gm_pred  = yield_pred * price_eur_t - n_kg_ha_s * n_cost_eur_kg -
                      COST_IRRIG * irrigated_s - COST_OTHER) |>
    select(scenario, region, crop, yield_pred, gm_pred, n_kg_ha_s)
}) |>
  mutate(scenario = factor(scenario, levels = scenarios$scenario))

tab_scen <- scen_results |>
  group_by(scenario, region) |>
  summarise(yield = mean(yield_pred), gross_margin = mean(gm_pred),
            N_use = mean(n_kg_ha_s), .groups = "drop") |>
  group_by(region) |>
  mutate(yield_change_pct = 100 * (yield / yield[scenario == "Baseline"] - 1),
         gm_change_eur_ha = gross_margin - gross_margin[scenario == "Baseline"]) |>
  ungroup()

print(filter(tab_scen, scenario != "Baseline"), n = 30)


# 7. Figures -------------------------------------------------------------

th <- theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

p1 <- ggplot(dat, aes(crop, yield_t_ha, fill = practice)) +
  geom_boxplot(outlier.alpha = 0.3) +
  scale_fill_viridis_d(end = 0.85) +
  labs(title = "A. Yield by crop and practice",
       x = NULL, y = "Yield (t/ha)", fill = NULL) + th

# nitrogen response curves, other variables held at their mean
n_grid <- expand_grid(crop = levels(dat$crop), n_100 = seq(0, 3, length.out = 60)) |>
  mutate(crop         = factor(crop, levels = levels(dat$crop)),
         temp_c       = mean(dat$temp_c),
         water_c      = mean(dat$water_c),
         soil_quality = mean(dat$soil_quality),
         practice     = factor("conventional", levels = levels(dat$practice)),
         trend        = mean(dat$trend))
n_grid$yield <- back_transform(m_yield, predict(m_yield, newdata = n_grid, re.form = NA))

p2 <- ggplot(n_grid, aes(n_100 * 100, yield, colour = crop)) +
  geom_line(linewidth = 1) +
  scale_colour_viridis_d(end = 0.85) +
  labs(title = "B. Modelled nitrogen response",
       x = "Nitrogen (kg/ha)", y = "Predicted yield (t/ha)", colour = NULL) + th

p3 <- dat |>
  group_by(region, year) |>
  summarise(rel_yield = mean(rel_yield), .groups = "drop") |>
  ggplot(aes(year, rel_yield)) +
  geom_line(colour = "grey60") +
  geom_smooth(method = "lm", se = TRUE, colour = "#2166ac", fill = "#2166ac", alpha = 0.15) +
  facet_wrap(~ region, nrow = 2) +
  labs(title = "C. Yield trend by region (relative to crop mean)",
       x = NULL, y = "Relative yield") + th

p4 <- tab_scen |>
  filter(scenario != "Baseline") |>
  ggplot(aes(scenario, region, fill = yield_change_pct)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%+.1f%%", yield_change_pct)), size = 3) +
  scale_fill_gradient2(low = "#b2182b", mid = "white", high = "#2166ac", midpoint = 0) +
  scale_x_discrete(labels = function(x) str_wrap(x, 14)) +
  labs(title = "D. Yield change relative to baseline",
       x = NULL, y = NULL, fill = "% change") + th

panel <- (p1 | p2) / (p3 | p4)
ggsave("outputs/figures/fig_panel.png",              panel, width = 14, height = 10, dpi = 300)
ggsave("outputs/figures/fig1_yield_by_practice.png", p1, width = 7, height = 5, dpi = 300)
ggsave("outputs/figures/fig2_nitrogen_response.png", p2, width = 7, height = 5, dpi = 300)
ggsave("outputs/figures/fig3_regional_trends.png",   p3, width = 9, height = 5, dpi = 300)
ggsave("outputs/figures/fig4_scenarios.png",         p4, width = 8, height = 5, dpi = 300)


# 8. Save tables and print main results ----------------------------------

write_csv(tab_region_crop, "outputs/tables/summary_region_crop.csv")
write_csv(tab_practice,    "outputs/tables/summary_practice.csv")
write_csv(tab_trend,       "outputs/tables/yield_trend_by_region.csv")
write_csv(tab_coef,        "outputs/tables/model_fixed_effects.csv")
write_csv(tab_recovery,    "outputs/tables/parameter_check.csv")
write_csv(tab_cv,          "outputs/tables/cross_validation.csv")
write_csv(tab_scen,        "outputs/tables/scenario_results.csv")

worst <- tab_scen |> filter(scenario == "+3 C, -10% rain") |> slice_min(yield_change_pct, n = 1)
adapt <- tab_scen |> filter(scenario == "+3 C, -10% rain + adaptation") |>
  summarise(mean_change = mean(yield_change_pct), mean_N = mean(N_use))

cat("\nMain results\n")
cat(sprintf("R2 marginal %.2f, conditional %.2f; cross-validated RMSE %.2f t/ha\n",
            r2_marg, r2_cond, mean(tab_cv$RMSE_t_ha)))
cat(sprintf("Yield-maximising N rate: about %.0f kg/ha\n", n_opt))
cat(sprintf("Largest loss under +3 C and -10%% rain: %s (%.1f%%)\n",
            worst$region, worst$yield_change_pct))
cat(sprintf("Adaptation package: mean yield change %+.1f%% with %.0f kg N/ha\n",
            adapt$mean_change, adapt$mean_N))

sessionInfo()
