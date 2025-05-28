### BII data form natural history museum



library(tidyverse)
library(readr)
library(dplyr)

# Load the data
bii_data <- readRDS("data/long_data.rds")
bii_data <- read.csv("data/resource.csv")
# View the first few rows
head(bii_data)

#check scenarios.- different model used for each scen projections
unique(bii_data$scenario)

# Check the structure of the dataset
str(bii_data)
bii_data <- bii_data %>%
  mutate(area_code = gsub("-", ".", area_code))

library(jsonlite)


# Read area code JSON, countries only
area_mapping<- fromJSON("data/bii.json")[c(1:128,153:279),]
# Check structure of the JSON
str(area_mapping)
names(area_mapping)<-c("area_code","area_name","parent_code")

# add area_names
bii_data_named <- bii_data %>%
  left_join(area_mapping, by = "area_code")


# Choose a few example countries
selected_countries <- c("India", "Indonesia", "Australia")

bii_subset <- bii_data_named %>%
  filter(area_name %in% selected_countries, variable=="bii")

ggplot(bii_subset, aes(x = year, y = value, color = scenario)) +
  geom_point() +
  facet_wrap(~ area_name, scales = "free_y") +
  labs(title = "BII Time Series by Country",
       x = "Year", y = "Biodiversity Intactness Index (%)") +
  theme_minimal()

### data are in 5 year blocks

library(readr)

# Export full cleaned BII data with country names
write_csv(bii_data_named, "output/bii_timeseries_by_country.csv")