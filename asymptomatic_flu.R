# Start with a clear environment
rm(list = ls())

# Asymptomatic vs symptomatic influenza cases
flu <- readRDS('./tmp/flu_p2_dat_proc.rds')
flu_pos <- flu |> filter(flu == 1)

# Drop full data
remove(flu)

# Create symptomatic variable
flu_pos$asymp <- ifelse(flu_pos$symp_count > 0, 0, 1)

# Asymptomatic percent over time
asymp_time <- flu_pos |> 
  group_by(week_date, asymp) |>
  mutate(ppl_sum = sum(patient_count_imp)) |>
  distinct(week_date, asymp, ppl_sum, .keep_all = FALSE) |>
  ungroup() |> group_by(week_date) |>
  mutate(total = sum(ppl_sum),
         percent = ppl_sum / total) |>
  filter(asymp == 1)

# Plot
coeff <- max(asymp_time$total) / max(asymp_time$percent)
ggplot(asymp_time, aes(x=week_date)) +
  geom_line(aes(y=percent), color = 'red') + 
  geom_line(aes(y=total / coeff), color = 'blue') + # Divide by 10 to get the same range than the temperatur
  scale_y_continuous(
    # Features of the first axis
    name = "Percent Asymptomatic",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff, name="All Flu Cases")
  )

# Asymptomatic percent by county
asymp_place <- flu_pos |> 
  group_by(county_fips, asymp) |>
  mutate(ppl_sum = sum(patient_count_imp)) |>
  distinct(county_fips, asymp, ppl_sum, .keep_all = FALSE) |>
  ungroup() |> group_by(county_fips) |>
  mutate(total = sum(ppl_sum),
         percent = ppl_sum / total) |>
  filter(asymp == 1)

# Plot
ggplot(asymp_place, aes(x=percent)) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8)

# Create a map
# Load maps
usa_albers_state <- st_as_sf(readRDS('tmp/usa_albers_state.rds'))   # convert to sf
usa_albers_county <- st_as_sf(readRDS('tmp/usa_albers_county.rds')) # convert to sf

# Create county to county group x-walk
x_walk <- flu_pos |> separate(county_fips, sep = "_", remove = FALSE, 
                                into = c('fips_1', 'fips_2', 'fips_3', 'fips_4',
                                         'fips_5', 'fips_6', 'fips_7', 'fips_8')) |>
  pivot_longer(cols=c('fips_1', 'fips_2', 'fips_3', 'fips_4', 'fips_5', 'fips_6', 
                      'fips_7', 'fips_8'),
               names_to='fips_num',
               values_to='fips') |>
  filter(!is.na(fips)) |> select(c(county_fips, fips))
x_walk$state_fips <- substr(x_walk$county_fips, 1, 2)
x_walk <- x_walk |> distinct(fips, county_fips)

# Collapse county map to county group
usa_albers_county_x_walk <- left_join(x_walk, usa_albers_county[, c(5, 6)], by = c('fips' = 'GEOID'))
usa_albers_county_group <- usa_albers_county_x_walk %>% 
  group_by(county_fips) %>%
  summarise(geometry = sf::st_union(geometry)) %>%
  ungroup()

# Merge on data
usa_albers_county_group <- left_join(usa_albers_county_group, asymp_place,
                                     by = c('county_fips' = 'county_fips'))
# Create maps
map <- ggplot() +
  geom_sf(data = st_as_sf(usa_albers_county_group), aes(fill = percent, group = county_fips), color= 'black', linewidth = 0.15) +
  geom_sf(data = usa_albers_state, aes(group = STATEFP), fill = '#FFFFFF00', color= 'black', linewidth = 0.4) +
  scale_fill_viridis_c('Percent', direction = -1) + ggtitle('Asymptomatic Percent') + 
  theme_void() + theme(legend.position = 'right',
                       plot.title = element_text(size = 20, hjust = 0.5),
                       panel.border = element_rect(fill=NA, linewidth = 0.8, color = 'white'),
                       legend.title = element_text(size = 12),
                       legend.text = element_text(size = 12),
                       legend.key.size = unit(0.6, 'cm')) 
map

# Differences by demographics
# Age
asymp_age <- flu_pos |> 
  group_by(age_grp, asymp) |>
  mutate(ppl_sum = sum(patient_count_imp)) |>
  distinct(age_grp, asymp, ppl_sum, .keep_all = FALSE) |>
  ungroup() |> group_by(age_grp) |>
  mutate(total = sum(ppl_sum),
         percent = ppl_sum / total) |>
  filter(asymp == 1)
asymp_age
# Sex
asymp_sex <- flu_pos |> 
  group_by(patient_gender_code, asymp) |>
  mutate(ppl_sum = sum(patient_count_imp)) |>
  distinct(patient_gender_code, asymp, ppl_sum, .keep_all = FALSE) |>
  ungroup() |> group_by(patient_gender_code) |>
  mutate(total = sum(ppl_sum),
         percent = ppl_sum / total) |>
  filter(asymp == 1)
asymp_sex
# Cases over time
asymp_timeseries <- flu_pos |> 
  group_by(week_date, asymp) |>
  mutate(ppl_sum = sum(patient_count_imp)) |>
  distinct(week_date, asymp, ppl_sum, .keep_all = FALSE) |>
  ungroup() |> group_by(week_date) |>
  mutate(total = sum(ppl_sum),
         percent = ppl_sum / total)

asymp_timeseries_county <- flu_pos |> 
  group_by(week_date, county_fips, asymp) |>
  mutate(ppl_sum = sum(patient_count_imp)) |>
  distinct(week_date, county_fips, asymp, ppl_sum, .keep_all = FALSE) |>
  ungroup() |> group_by(week_date, county_fips) |>
  mutate(total = sum(ppl_sum),
         percent = ppl_sum / total)

asymp_timeseries_county_lim <- asymp_timeseries_county |> 
  filter(county_fips == "17031" | county_fips == "53033" | 
           county_fips == "06037" | county_fips == "04013" |
           county_fips == "36061" | county_fips == "48113" |
           county_fips == "11001" | county_fips == "13121" |
           county_fips == "29510")

# Plot
ggplot(asymp_timeseries, aes(x=week_date)) +
  geom_line(aes(y=ppl_sum, color = as.factor(asymp), group = as.factor(asymp)))

ggplot(asymp_timeseries_county_lim) +
  geom_line(aes(x = week_date, y = ppl_sum, color = as.factor(asymp), group = as.factor(asymp))) +
  facet_wrap(vars(county_fips), scales = "free", nrow = 5) 

pop <- read.csv('raw/uscounties.csv')
pop$county_fips <- as.character(pop$county_fips)
pop$county_fips <- ifelse(nchar(pop$county_fips) == 4, paste0('0', pop$county_fips), pop$county_fips)
pop <- left_join(x_walk, pop, by = c('fips' = 'county_fips'))
pop <- pop |> group_by(county_fips) |>
  mutate(pop_sum = sum(population)) |>
  distinct(county_fips, pop_sum, .keep_all = FALSE)

asymp_place <- left_join(asymp_place, pop, by = c('county_fips' = 'county_fips'))
ggplot(asymp_place) + geom_point(aes(y = percent, x = log(pop_sum)))
