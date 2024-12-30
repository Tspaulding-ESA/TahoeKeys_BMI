#################### TAHOE KEYS BMI ANALYSIS ###############################

# Setup ####################################################################
source(file.path("2_Code","0_Setup.R"))

trt_color = read.csv(file.path("1_Data", "TreatmentColors.csv")) %>%
  mutate(trt = factor(trt, levels = trt))

# Get Data #################################################################
file_list <- list.files(file.path("Z:","Shared","Projects","2022",
                                  "D202200318.00 - TRPA Tahoe Keys Lagoons Aquatic Weeds CMT",
                                  "10 Data Collection and Management","2_BMI","BMI_Master"),
                        pattern = "CMT BMI", full.names = TRUE)

abundance <- data.frame()
traits <- data.frame()
for(i in 1:length(file_list)){
  abund <- readxl::read_excel(file_list[i], sheet = "Long output") |>
    fix_names() |> 
    rename("site" = replicate) |>
    mutate( date = lubridate::ymd(date),
            site = as.numeric(site))
  trait <- readxl::read_excel(file_list[i], sheet = "Traits") |> 
    fix_names()
  abundance <- bind_rows(abundance, abund)
  traits <- bind_rows(traits, trait) |> distinct()
  abundance
  traits
}

treatment_schedule <- read.csv(file.path("1_Data","treatments_rev.csv"))

dat <- 
  abundance |>
  mutate(year = year(date),) |>
  right_join(treatment_schedule) |>
  mutate(season = case_when(
           month(date) %in% c(9,10,11) ~ "Fall",
           month(date) %in% c(4,5,6) ~ "Spring",
           TRUE ~ "NA"),
         period = case_when(
           year <= 2021 ~ "Before",
           year == 2022 & season == "Spring" ~ "Before",
           year == 2022 & season == "Fall" ~ "After",
           year == 2023 ~ "After",
           year == 2024 ~ "After")
         ) |>
  left_join(traits)

# Make some adjustments per Toni
dat <- 
  dat |>
  # Lets just use the 2022 and 23 data for now. 
  # Our meeting with Dr. Lars Anderson (science advisor to the project) 
  # will better inform whether we use previous data or not. 
  filter(year >= 2022) |>
  # In 2022, the "before" samples were collected in Oct in sites 7, 25, and 27.
  mutate(period = case_when(
    site %in% c(7,25,27) & year == 2022 & season == "Fall" ~ "Before",
    site == 26 ~ "After",
    TRUE ~ period
  )) |>
  filter(site != 21) %>%
  filter(!is.na(habitat))

sample_dates <- dat %>%
  select(site,date) %>%
  distinct() %>%
  mutate(year = year(date)) %>%
  pivot_wider(names_from = "year", values_from = date) %>%
  mutate(across(`2022`:last_col(), .fns = function(x) as_date(as.numeric(as.character(x)),
                                                              origin = origin)))

write.csv(sample_dates, "SampleDates.csv")

# in 2023 Sites 5, 10, 12, 14, 21 are all zero counts
zero_counts <- dat %>%
  filter(year(date) == 2023, site %in% c(5,10,12,14,21)) %>%
  select(waterbody:date, location:period) %>%
  mutate(habitat = "mid-channel") %>%
  distinct()

dat <- dat %>%
  bind_rows(zero_counts) %>%
  mutate(inverts = case_when(
    year(date) %in% c(2021,2022,2024) ~ TRUE,
    year(date) %in% c(2023) & 
      site %in% c(5,10,12,13,21) ~ FALSE,
    TRUE ~ TRUE
  )) %>%
  mutate(period = factor(period, levels = c("Before","After"
  )),
  habitat = factor(habitat, levels = c("near-shore","mid-channel")))

# Summary #################################################################
summary <- dat |>
  group_by(waterbody, location, year, area, period, site, treatment, habitat) |>
  summarise(all_rich = length(unique(taxon[!is.na(taxon)])),
            ept_rich = sum(ifelse(order %in% c("Ephemeroptera","Plecoptera",
                                               "Trichoptera"),1,0)),
            perc_dom = max(abundance)/sum(abundance),
            dom_taxa = paste(list(taxon[abundance == max(abundance)])),
            avg_tol = sum((abundance[which(wy.hbi != 11)]/
                             sum(abundance[which(wy.hbi != 11)]))*
                            wy.hbi[which(wy.hbi != 11)]),
            perc_scraper = sum(ifelse(feeding.group == "SC",1,0)*abundance)/sum(abundance),
            perc_predator = sum(ifelse(feeding.group == "PR",1,0)*abundance)/sum(abundance),
            perc_filter = sum(ifelse(feeding.group == "CF",1,0)*abundance)/sum(abundance),
            perc_gather = sum(ifelse(feeding.group == "CG",1,0)*abundance)/sum(abundance),
            perc_other = sum(ifelse(feeding.group %notin% c("SC","PR","CF","CG"),
                                    1,0)*abundance)/sum(abundance)
  )

summary_table <- dat |>
  ungroup() |>
  group_by(area, treatment, habitat, period, year) |>
  summarise(sites = as.character(list(unique(site))),
            all_rich = length(unique(taxon[!is.na(taxon)])),
            ept_rich = sum(ifelse(order %in% c("Ephemeroptera","Plecoptera",
                                               "Trichoptera"),1,0)),
            perc_dom = max(abundance, na.rm = TRUE)/sum(abundance, na.rm = TRUE),
            dom_taxa = paste(list(taxon[abundance == max(abundance, na.rm = TRUE) &
                                          !is.na(taxon)])),
            avg_tol = sum((abundance[which(wy.hbi != 11)]/
                             sum(abundance[which(wy.hbi != 11)]))*
                            wy.hbi[which(wy.hbi != 11)]),
            perc_scraper = sum(ifelse(feeding.group == "SC",1,0)*abundance, 
                               na.rm = TRUE)/sum(abundance, na.rm = TRUE),
            perc_predator = sum(ifelse(feeding.group == "PR",1,0)*abundance, na.rm = TRUE)/sum(abundance, na.rm = TRUE),
            perc_filter = sum(ifelse(feeding.group == "CF",1,0)*abundance, na.rm = TRUE)/sum(abundance, na.rm = TRUE),
            perc_gather = sum(ifelse(feeding.group == "CG",1,0)*abundance, na.rm = TRUE)/sum(abundance, na.rm = TRUE),
            perc_other = sum(ifelse(feeding.group %notin% c("SC","PR","CF","CG"),
                                    1,0)*abundance, na.rm = TRUE)/sum(abundance, na.rm = TRUE)
  )
  

write.csv(summary_table, file.path("1_Data","Output","Summary_Table.csv"))

# Plots ###################################################################
## Richness ===============================================================
richness_hab <- summary |>
  dplyr::select(waterbody:ept_rich) |>
  pivot_longer(cols = all_rich:ept_rich, names_to = "richness", 
               values_to = "value") |>
  mutate(richness = case_when(
    richness == "all_rich" ~ "Overall",
    richness == "ept_rich" ~ "EPT"
    ),
    ) %>%
  # ungroup() |>
  # group_by(waterbody,location,area,treatment) |>
  # mutate(site_list = list(unique(site))) |>
  mutate(treatment_area = paste0(area,"\n",
                                 #site_list,"\n",
                                 treatment)) |>
  mutate(year = stringr::str_sub(as.character(year),3,4)) %>%
  mutate(pd_yr = paste(period,year, sep = "-'"))

ggpubr::ggbarplot(data = richness_hab |>
                    filter(richness == "Overall"),
                  x = "pd_yr", y = "value",
                  add = "mean_se_",
                  fill = "treatment",
                  facet.by = c("habitat","treatment_area"))+
  stat_compare_means(comparisons = list(c("Before-'22","After-'23"),
                                        c("Before-'22","After-'24")),
                     label = "p.format",
                     size = 3)+
  labs(title = "Overall Richness", x = "Treatment", 
       y = "Average Richness \u00B1 se")+
  scale_fill_manual(breaks = trt_color$trt,
                    values = trt_color$color_cmt)+
  theme_classic()+
  theme(
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5)
  )+
  coord_cartesian(ylim = c(0,45))

ggsave(file.path("1_Data","Output","Overall_Richness_by_Habitat_Plot_yearsplit.png"), device = "png",
       width = 15, height = 7, units = "in", dpi = 300)

ggpubr::ggbarplot(data = richness_hab |>
                    filter(richness == "Overall"),
                  x = "period", y = "value",
                  add = "mean_se_",
                  fill = "treatment",
                  facet.by = c("habitat","treatment_area"))+
  stat_compare_means(comparisons = list(c("Before","After")),
                     label = "p.format",
                     size = 3)+
  labs(title = "Overall Richness", x = "Treatment", 
       y = "Average Richness \u00B1 se")+
  scale_fill_manual(breaks = trt_color$trt,
                    values = trt_color$color_cmt)+
  theme_classic()+
  theme(
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5)
  )+
  coord_cartesian(ylim = c(0,45))

ggsave(file.path("1_Data","Output","Overall_Richness_by_Habitat_Plot_periodsplit.png"), device = "png",
       width = 15, height = 7, units = "in", dpi = 300)

ggpubr::ggbarplot(data = richness_hab |>
                    filter(richness == "EPT"),
                  x = "pd_yr", y = "value",
                  fill = "treatment",
                  add = "mean_se_",
                  facet.by = c("habitat","treatment_area"))+
  stat_compare_means(comparisons = list(c("Before-'22","After-'23"),
                                       c("Before-'22","After-'24")),
                    label = "p.format",
                    size = 3)+
  labs(title = "EPT Richness", x = "Treatment", 
       y = "Average Richness \u00B1 se")+
  scale_fill_manual(breaks = trt_color$trt,
                    values = trt_color$color_cmt)+
  scale_y_continuous(breaks = seq(0,10,2))+
  theme_classic()+
  theme(
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5)
  )+
  coord_cartesian(ylim = c(0,9))

ggsave(file.path("1_Data","Output","EPT_Richness_by_Habitat_Plot_yearsplit.png"), device = "png",
       width = 15, height = 7, units = "in", dpi = 300)

ggpubr::ggbarplot(data = richness_hab |>
                    filter(richness == "EPT"),
                  x = "period", y = "value",
                  fill = "treatment",
                  add = "mean_se_",
                  facet.by = c("habitat","treatment_area"))+
  stat_compare_means(comparisons = list(c("Before","After")),
                     label = "p.format",
                     size = 3)+
  labs(title = "EPT Richness", x = "Treatment", 
       y = "Average Richness \u00B1 se")+
  scale_fill_manual(breaks = trt_color$trt,
                    values = trt_color$color_cmt)+
  theme_classic()+
  theme(
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5)
  )+
  coord_cartesian(ylim = c(0,8.5))

ggsave(file.path("1_Data","Output","EPT_Richness_by_Habitat_Plot_periodsplit.png"), device = "png",
       width = 15, height = 7, units = "in", dpi = 300)


## ANOVA ====================================================================
richness_anova_hab <- aov(cbind(all_rich,ept_rich)  ~  area * period * treatment * habitat, summary)
summary(richness_anova_hab)

## Dominant Taxa ===============================================================
dominant_hab <- summary |>
  dplyr::select(waterbody:habitat,perc_dom) |>
  mutate(treatment_area = paste0(area,"\n",
                                 #site_list,"\n",
                                 treatment)) |>
  mutate(year = stringr::str_sub(as.character(year),3,4)) %>%
  mutate(pd_yr = paste(period,year, sep = "-'"))

ggpubr::ggbarplot(data = dominant_hab,
                  x = "pd_yr", y = "perc_dom",
                  add = "mean_se_",
                  fill = "treatment",
                  facet.by = c("habitat","treatment_area"))+
  stat_compare_means(comparisons = list(c("Before-'22","After-'23"),
                                        c("Before-'22","After-'24")),
                     label = "p.format",
                     size = 3)+
  labs(title = "Dominant Taxa", x = "Treatment", 
       y = "Average Dominant Taxa %-Abundance \u00B1 se")+
  scale_fill_manual(breaks = trt_color$trt,
                    values = trt_color$color_cmt)+
  scale_y_continuous(labels = scales::percent)+
  theme_classic()+
  theme(
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5)
  )

ggsave(file.path("1_Data","Output","Dominant_Taxa_by_Habitat_Plot_yearsplit.png"), device = "png",
       width = 15, height = 7, units = "in", dpi = 300)

ggpubr::ggbarplot(data = dominant_hab,
                  x = "period", y = "perc_dom",
                  add = "mean_se_",
                  fill = "treatment",
                  facet.by = c("habitat","treatment_area"))+
  stat_compare_means(comparisons = list(c("Before","After")),
                     label = "p.format",
                     size = 3)+
  labs(title = "Overall Dominant Taxa", x = "Treatment", 
       y = "Average Dominant Taxa \u00B1 se")+
  scale_fill_manual(breaks = trt_color$trt,
                    values = trt_color$color_cmt)+
  theme_classic()+
  theme(
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5)
  )+
  scale_y_continuous(labels = scales::percent)

ggsave(file.path("1_Data","Output","Dominant_Taxa_by_Habitat_Plot_periodsplit.png"), device = "png",
       width = 15, height = 7, units = "in", dpi = 300)

dominant_anova_hab <- aov(perc_dom  ~  area * period * treatment * habitat, summary)
summary(dominant_anova_hab)

dominance_summary <- dominant_hab %>%
  group_by(area, treatment, period) %>%
  summarise(mean = mean(perc_dom, na.rm = TRUE),
            min = min(perc_dom, na.rm = TRUE),
            max = max(perc_dom, na.rm = TRUE)) %>%
  print(n = 50)

# Feeding Garea# Feeding Group Analysis ###################################################
feed_groups <- summary |>
  dplyr::select(waterbody:habitat, perc_scraper:perc_other) |>
  pivot_longer(cols = perc_scraper:perc_other, names_to = "feed", 
               names_prefix = "perc_", values_to = "percent") |>
  ungroup() |>
  group_by(waterbody,location,area,treatment) |>
  mutate(site_list = list(unique(site))) |>
  mutate(treatment_area = paste0(area,"\n",
                                 #"Sites: ",site_list,"\n",
                                 treatment))
  # group_by(area, period, treatment, habitat, feed) |>
  # summarise(mean = mean(percent, na.rm = TRUE),
  #           min = min(percent, na.rm = TRUE),
  #           max = max(percent, na.rm = TRUE)) %>%
  # filter(mean > 0)

ggpubr::ggbarplot(data = feed_groups,
                  x = "period", y = "percent",
                  add = "mean", fill = "feed",
                  facet.by = c("habitat","treatment_area"))+
  stat_compare_means(comparisons = list(c("Before","After")),
                     label = "p.format", method = "wilcox.test",
                     size = 3)+
  labs(title = "Feeding Groups", x = "Treatment", 
       y = "Average Percent Abundance")+
  scale_fill_discrete(name = "Feeding Group",
                      breaks = c("predator","scraper", "gather", 
                                 "filter","other"),
                      labels = c("Predators", "Scrapers", 
                                 "Collector - Gatherers", 
                                 "Collector - Filterers",
                                 "Other"))+
  theme_classic()+
  theme(
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5)
  )+
  coord_cartesian(ylim = c(0,1.1))

ggsave(file.path("1_Data","Output","FeedingGroup_by_Habitat_Plot.png"), device = "png",
       width = 15, height = 7, units = "in", dpi = 300)  


## ANOVAs =================================================================
feeding_anova_hab <- aov(cbind(perc_predator, perc_scraper,
                               perc_gather, perc_filter, perc_other) ~ area * period * treatment * habitat,
                         summary %>% filter(!is.na(perc_dom)))
summary(feeding_anova_hab)

# Tolerance Analysis ###################################################
tolerance <- summary |>
  dplyr::select(waterbody:habitat,avg_tol) |>
  mutate(treatment_area = paste0(area,"\n",
                                 #site_list,"\n",
                                 treatment)) |>
  mutate(year = stringr::str_sub(as.character(year),3,4)) %>%
  mutate(pd_yr = paste(period,year, sep = "-'"))
  
# group_by(area, period, treatment, habitat, feed) |>
# summarise(mean = mean(percent, na.rm = TRUE),
#           min = min(percent, na.rm = TRUE),
#           max = max(percent, na.rm = TRUE)) %>%
# filter(mean > 0)

ggpubr::ggboxplot(data = tolerance,
                  x = "period", y = "avg_tol",
                  fill = "treatment",
                  facet.by = c("habitat","treatment_area"))+
  stat_compare_means(comparisons = list(c("Before","After")),
                     label = "p.format",
                     size = 3)+
  #geom_hline(yintercept = 7, color = "blue", linewidth = 0.75)+
  labs(title = "Tolerance", x = "Treatment", 
       y = "Tolerance")+
  scale_fill_manual(breaks = trt_color$trt,
                    values = trt_color$color_cmt)+
  theme_classic()+
  theme(
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5)
  )+
  coord_cartesian(ylim = c(0,11.5))

ggsave(file.path("1_Data","Output","Tolerance_by_Habitat_Plot.png"), device = "png",
       width = 15, height = 7, units = "in", dpi = 300)  

ggpubr::ggboxplot(data = tolerance,
                  x = "pd_yr", y = "avg_tol",
                  fill = "treatment",
                  facet.by = c("habitat","treatment_area"))+
  stat_compare_means(comparisons = list(c("Before-'22","After-'23"),
                                        c("Before-'22","After-'24")),
                     label = "p.format",
                     size = 3)+
  #geom_hline(yintercept = 7, color = "blue", linewidth = 0.75)+
  labs(title = "Tolerance", x = "Treatment", 
       y = "Tolerance")+
  scale_fill_manual(breaks = trt_color$trt,
                    values = trt_color$color_cmt)+
  theme_classic()+
  theme(
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5)
  )+
  coord_cartesian(ylim = c(0,12.5))

ggsave(file.path("1_Data","Output","Tolerance_by_Habitat_Plot_yearsplit.png"), device = "png",
       width = 15, height = 7, units = "in", dpi = 300)  


## ANOVAs ==================================================================
tolerance_anova_hab <- aov(avg_tol ~ area * period * treatment * habitat,
                         tolerance)
summary(tolerance_anova_hab)

## IBI ################################################################


ibi <- dat |>  
  arrange(waterbody, location, year, area, period, site, treatment, habitat, desc(abundance)) |>
  group_by(waterbody, location, year, area, period, site, treatment, habitat) |>
  mutate(dom_taxa = taxon %in% (unique(taxon[!is.na(taxon)])[1:3]) &
           !is.na(taxon)) %>%
  summarise(all_rich = length(unique(taxon[!is.na(taxon)])),
            total_abund = case_when(
              !is.na(sum(abundance)) ~ sum(abundance),
              TRUE ~ 0),
            ephem_rich = sum(ifelse(order == "Ephemeroptera",1,0)),
            plec_rich = sum(ifelse(order == "Plecoptera",1,0)),
            tric_rich = sum(ifelse(order == "Tichoptera",1,0)),
            acari_rich = sum(ifelse(higher.classification == "Arachnida: Acari",1,0)),
            perc_midge = (sum(ifelse(grepl("Chironomidae",family),1,0))/all_rich)*100,
            perc_tol = (sum(ifelse(wy.hbi %in% c(7:10),1,0))/all_rich)*100,
            perc_shred = (sum(ifelse(feeding.group %in% c("SH"),abundance,0))/total_abund)*100,
            perc_dom = (sum(case_when(
              dom_taxa ~ abundance,
              all_rich == 1 ~ abundance,
              TRUE ~ 0
            ))/total_abund)*100,
            com_tol = sum((abundance[which(wy.hbi != 11)]/
                             sum(abundance[which(wy.hbi != 11)]))*
                            wy.hbi[which(wy.hbi != 11)]),
            pred_rich = sum(ifelse(feeding.group == "PR",1,0)),
            ept_abund = (sum(ifelse(order %in% c("Ephemeroptera","Tricoptera",
                                                 "Plecoptera"),abundance,0))/total_abund)*100
  ) %>%
  replace_na(list(
    "perc_shred" = 0,
    "com_tol" = 10
  )) %>%
  mutate(
    rich_code = case_when(
      all_rich <= 30 ~ 0,
      all_rich >= 50 ~ 10,
      TRUE ~ 10*(all_rich-30)/20
    ),
    eph_code = case_when(
      ephem_rich <= 3 ~ 0,
      ephem_rich >= 9 ~ 10,
      TRUE ~ 10*(ephem_rich-3)/6
    ),
    plec_code = case_when(
      plec_rich <= 1 ~ 0,
      plec_rich >= 6 ~ 10,
      TRUE ~ 10*(plec_rich-1)/5
    ),
    tric_code = case_when(
      tric_rich <= 2 ~ 0,
      tric_rich >= 8 ~ 10,
      TRUE ~ 10*(tric_rich-2)/6
    ),
    acari_code = case_when(
      acari_rich <= 1 ~ 0,
      acari_rich >= 6 ~ 10,
      TRUE ~ 10*(acari_rich-1)/5
    ),
    midge_code = case_when(
      perc_midge >= 43.4 ~ 0,
      perc_midge <= 26.4 ~ 10,
      TRUE ~ 10*(43.4 - perc_midge)/17
    ),
    tol_code = case_when(
      perc_tol >= 34.1 ~ 0,
      perc_tol <= 18.7 ~ 10,
      TRUE ~ 10*(34.1-perc_tol)/(15.4)
    ),
    shred_code = case_when(
      perc_shred == 0 ~ 0,
      perc_shred >= 2.7 ~ 10,
      TRUE ~ 10*perc_shred/(2.7)
    ),
    dom_code = case_when(
      perc_dom >= 65.9 ~ 0,
      perc_dom <= 42.9 ~ 10,
      TRUE ~ 10*(65.9-perc_dom)/(23)
    ),
    comm_tol_code = case_when(
      com_tol >= 5.79 ~ 0,
      com_tol <= 4.05 ~ 10,
      TRUE ~ 10*(5.79-com_tol)/(5.79-4.05)
    ),
    pred_code = case_when(
      pred_rich <= 7 ~ 0,
      pred_rich >= 16 ~ 10,
      TRUE ~ 10*(pred_rich - 7)/(16-7)
    ),
    ept_abund_code = case_when(
      ept_abund <= 17.5 ~ 0,
      ept_abund >= 59.1 ~ 10,
      TRUE ~ 10*(ept_abund - 17.5)/(59.1-17.5)
    ))|>
  mutate(treatment_area = paste0(area,"\n",
                                 #site_list,"\n",
                                 treatment)) |>
  mutate(year = stringr::str_sub(as.character(year),3,4)) %>%
  mutate(pd_yr = paste(period,year, sep = "-'")) %>%
  replace_na(list("rich_code" = 0,
             "eph_code" = 0,
             "plec_code" = 0,
             "tric_code" = 0,
             "acari_code" = 0,
             "midge_code" = 0,
             "tol_code" = 0,
             "shred_code" = 0,
             "dom_code" = 0,
             "comm_tol_code" = 0,
             "pred_code" = 0,
             "ept_abund_code" = 0)) %>%
  mutate(
         ibi_sum = rich_code + eph_code + plec_code + tric_code + acari_code + 
           midge_code + tol_code + shred_code + dom_code + comm_tol_code + 
           pred_code + ept_abund_code,
         mod_ibi_sum = 2*(rich_code + tol_code + dom_code + comm_tol_code + pred_code)
  )

ibi_control <- dat |>
  filter(period == "Before") |>
  arrange(waterbody, location, year, area, period, site, treatment, habitat, desc(abundance)) |>
  group_by(waterbody, location, year, area, period, site, treatment, habitat) |>
  mutate(dom_taxa = taxon %in% (unique(taxon)[1:3])) %>%
  summarise(all_rich = length(unique(taxon[!is.na(taxon)])),
            total_abund = case_when(
              !is.na(sum(abundance)) ~ sum(abundance),
              TRUE ~ 0),
            ephem_rich = sum(ifelse(order == "Ephemeroptera",1,0)),
            plec_rich = sum(ifelse(order == "Plecoptera",1,0)),
            tric_rich = sum(ifelse(order == "Tichoptera",1,0)),
            acari_rich = sum(ifelse(higher.classification == "Arachnida: Acari",1,0)),
            perc_midge = (sum(ifelse(grepl("Chironomidae",family),1,0))/all_rich)*100,
            perc_tol = (sum(ifelse(wy.hbi %in% c(7:10),1,0))/all_rich)*100,
            perc_shred = (sum(ifelse(feeding.group %in% c("SH"),abundance,0))/total_abund)*100,
            perc_dom = (sum(case_when(
              dom_taxa ~ abundance,
              all_rich == 1 ~ abundance,
              TRUE ~ 0
            ))/total_abund)*100,
            com_tol = sum((abundance[which(wy.hbi != 11)]/
                             sum(abundance[which(wy.hbi != 11)]))*
                            wy.hbi[which(wy.hbi != 11)]),
            pred_rich = sum(ifelse(feeding.group == "PR",1,0)),
            ept_abund = (sum(ifelse(order %in% c("Ephemeroptera","Tricoptera",
                                                 "Plecoptera"),abundance,0))/total_abund)*100
  ) %>%
  replace_na(list(
    "perc_shred" = 0,
    "com_tol" = 10
  )) %>%
  ungroup() %>%
  select(all_rich, ephem_rich:ept_abund)

quantile(ibi_control$all_rich, probs = c(0.1,0.9), type = 3)
quantile(ibi_control$perc_midge, probs = c(0.1,0.9), type = 3)
quantile(ibi_control$perc_tol, probs = c(0.1,0.9), type = 3)
quantile(ibi_control$perc_dom, probs = c(0.1,0.9), type = 3, na.rm = TRUE)
quantile(ibi_control$com_tol, probs = c(0.1,0.9), type = 3)
quantile(ibi_control$pred_rich, probs = c(0.1,0.9), type = 3)

mod_ibi <- dat |>  
  arrange(waterbody, location, year, area, period, site, treatment, habitat, desc(abundance)) |>
  group_by(waterbody, location, year, area, period, site, treatment, habitat) |>
  mutate(dom_taxa = taxon %in% (unique(taxon)[1:3])) %>%
  summarise(all_rich = length(unique(taxon[!is.na(taxon)])),
            total_abund = case_when(
              !is.na(sum(abundance)) ~ sum(abundance),
              TRUE ~ 0),
            ephem_rich = sum(ifelse(order == "Ephemeroptera",1,0)),
            plec_rich = sum(ifelse(order == "Plecoptera",1,0)),
            tric_rich = sum(ifelse(order == "Tichoptera",1,0)),
            acari_rich = sum(ifelse(higher.classification == "Arachnida: Acari",1,0)),
            perc_midge = (sum(ifelse(grepl("Chironomidae",family),1,0))/all_rich)*100,
            perc_tol = (sum(ifelse(wy.hbi %in% c(7:10),1,0))/all_rich)*100,
            perc_shred = (sum(ifelse(feeding.group %in% c("SH"),abundance,0))/total_abund)*100,
            perc_dom = (sum(case_when(
              dom_taxa ~ abundance,
              all_rich == 1 ~ abundance,
              TRUE ~ 0
            ))/total_abund)*100,
            com_tol = sum((abundance[which(wy.hbi != 11)]/
                             sum(abundance[which(wy.hbi != 11)]))*
                            wy.hbi[which(wy.hbi != 11)]),
            pred_rich = sum(ifelse(feeding.group == "PR",1,0)),
            ept_abund = (sum(ifelse(order %in% c("Ephemeroptera","Tricoptera",
                                                 "Plecoptera"),abundance,0))/total_abund)*100
  ) %>%
  replace_na(list(
    "perc_shred" = 0,
    "com_tol" = 10
  )) %>%
  mutate(
    rich_code = case_when(
      all_rich <= 3 ~ 0,
      all_rich >= 23 ~ 10,
      TRUE ~ 10*(all_rich-3)/(23-3)
    ),
    midge_code = case_when(
      perc_midge >= 67 ~ 0,
      perc_midge <= 29 ~ 10,
      TRUE ~ 10*(67 - perc_midge)/(67-29)
    ),
    tol_code = case_when(
      perc_tol >= 65 ~ 0,
      perc_tol <= 33 ~ 10,
      TRUE ~ 10*(65-perc_tol)/(65-33)
    ),
    dom_code = case_when(
      perc_dom >= 100 ~ 0,
      perc_dom <= 47 ~ 10,
      TRUE ~ 10*(100-perc_dom)/(100-47)
    ),
    comm_tol_code = case_when(
      com_tol >= 8.1 ~ 0,
      com_tol <= 6 ~ 10,
      TRUE ~ 10*(8.1-com_tol)/(8.1 - 6)
    ),
    pred_code = case_when(
      pred_rich <= 1 ~ 0,
      pred_rich >= 9 ~ 10,
      TRUE ~ 10*(pred_rich - 1)/(9-1)
    ),
    ibi_sum = rich_code + midge_code + tol_code + dom_code + comm_tol_code + 
      pred_code
  )|>
  mutate(treatment_area = paste0(area,"\n",
                                 #site_list,"\n",
                                 treatment)) |>
  mutate(year = stringr::str_sub(as.character(year),3,4)) %>%
  mutate(pd_yr = paste(period,year, sep = "-'"))

ggpubr::ggboxplot(data = mod_ibi,
                  x = "pd_yr", y = "ibi_sum",
                  fill = "treatment",
                  facet.by = c("habitat","treatment_area"))+
  stat_compare_means(comparisons = list(c("Before-'22","After-'23"),
                                        c("Before-'22","After-'24")),
                     label = "p.format",
                     size = 3)+
  #geom_hline(yintercept = 7, color = "blue", linewidth = 0.75)+
  labs(title = "Index of Biological Integrity", x = "Treatment", 
       y = "Index of Biological Integrity")+
  scale_fill_manual(breaks = trt_color$trt,
                    values = trt_color$color_cmt)+
  theme_classic()+
  theme(
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5)
  )

ggsave(file.path("1_Data","Output","IBI_by_Habitat_Plot_yrsplit.png"), device = "png",
       width = 15, height = 7, units = "in", dpi = 300)  

ggpubr::ggboxplot(data = mod_ibi,
                  x = "period", y = "ibi_sum",
                  fill = "treatment",
                  facet.by = c("habitat","treatment_area"))+
  stat_compare_means(comparisons = list(c("Before","After")),
                     label = "p.format",
                     size = 3)+
  #geom_hline(yintercept = 7, color = "blue", linewidth = 0.75)+
  labs(title = "Index of Biological Integrity", x = "Treatment", 
       y = "Index of Biological Integrity")+
  scale_fill_manual(breaks = trt_color$trt,
                    values = trt_color$color_cmt)+
  theme_classic()+
  theme(
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5)
  )

ggsave(file.path("1_Data","Output","IBI_by_Habitat_Plot_periodsplit.pdf"), device = "png",
       width = 15, height = 7, units = "in", dpi = 300)  

ibi_anova_hab <- aov(ibi_sum ~ area * period * treatment * habitat,
                     mod_ibi)
summary(ibi_anova_hab)

# EXPORTING ################################################################
library(officer)
library(flextable)

# Function to convert aov output to ANOVA table
aov_output_to_table <- function(aov_output) {
  # Initialize empty list to store ANOVA tables for each response variable
  anova_tables <- list()
  
  # Iterate over each response variable
  for (i in 1:length(aov_output)) {
    # Extract ANOVA summary for the current response variable
    aov_summary <- summary(aov_output)[i]
    
    # Extracting the response variable name from the summary
    response_name <- colnames(aov_output$coefficients)[i]
    
    # Extracting the ANOVA table from the summary
    anova_table <- aov_summary[[1]]
    anova_table$Source <- as.character(row.names(anova_table))
    
    # Append response variable name to the ANOVA table
    anova_table$Response <- response_name
    
    # Append ANOVA table to the list
    anova_tables[[i]] <- anova_table
  }
  
  # Combine all ANOVA tables into a single data frame
  combined_anova_table <- bind_rows(anova_tables)
  
  combined_anova_table <- combined_anova_table %>%
    as_tibble() %>%
    tibble::rownames_to_column()
  
  # Return the combined ANOVA table
  return(combined_anova_table)
}

feeding_group_table <- aov_output_to_table(feeding_anova_hab)
richness_table <- aov_output_to_table(richness_anova_hab)
dominance_table <- aov_output_to_table(dominant_anova_hab)
tolerance_table <- aov_output_to_table(tolerance_anova_hab)
ibi_table <- aov_output_to_table(ibi_anova_hab)

feeding_group_table <-feeding_group_table %>%
  dplyr::select(Response, Source,Df:`Pr(>F)`) %>%
  mutate(across(where(is.numeric), .fns = ~round(.x,3)))

richness_table <-richness_table %>%
  dplyr::select(Response, Source,Df:`Pr(>F)`) %>%
  mutate(across(where(is.numeric), .fns = ~round(.x,3)))

dominance_table <-dominance_table %>%
  dplyr::select(Source,Df:`Pr(>F)`) %>%
  mutate(across(where(is.numeric), .fns = ~round(.x,3)))

tolerance_table <-tolerance_table %>%
  dplyr::select(Source,Df:`Pr(>F)`) %>%
  mutate(across(where(is.numeric), .fns = ~round(.x,3)))

ibi_table <-ibi_table %>%
  dplyr::select(Source,Df:`Pr(>F)`) %>%
  mutate(across(where(is.numeric), .fns = ~round(.x,3)))


# Convert data frame to flextable
set_flextable_defaults(font.family = "Arial", font.size = 10, color = "black",
                       text.align = "center", padding.top = 0,
                       padding.bottom = 0, 
                       line_spacing = 1, table.layout = "autofit",
                       digits = 2, pct_digits = 1)
feeding_group_ft <- flextable(feeding_group_table)
richness_ft <- flextable(richness_table)
tolerance_ft <- flextable(tolerance_table)
dominance_ft <- flextable(dominance_table)
ibi_ft <- flextable(ibi_table)

# Save as Word document
doc <- read_docx()
doc %>%
  body_add_flextable(feeding_group_ft, align = "center",
                     ) %>%
  print(target = file.path("1_Data","Output",
                           "Feeding_Group_ANOVA_Table.docx"))

doc <- read_docx()
doc %>%
  body_add_flextable(dominance_ft, align = "center",
  ) %>%
  print(target = file.path("1_Data","Output",
                           "Dominance_ANOVA_Table.docx"))

doc <- read_docx()
doc %>%
  body_add_flextable(richness_ft) %>%
  print(target = file.path("1_Data","Output",
                           "Richness_ANOVA_Table.docx"))
doc <- read_docx()
doc %>%
  body_add_flextable(tolerance_ft) %>%
  print(target = file.path("1_Data","Output",
                           "tolerance_ANOVA_Table.docx"))

doc <- read_docx()
doc %>%
  body_add_flextable(ibi_ft, align = "center",
  ) %>%
  print(target = file.path("1_Data","Output",
                           "IBI_ANOVA_Table.docx"))

print("ANOVA table has been saved as 'ANOVA_Table.docx'")
