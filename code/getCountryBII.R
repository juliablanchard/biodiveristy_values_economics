######### 1. Get biodiveristy model data: using Biologicla Intactness Index data from natural history museum (PREDICTS)
# Data portal:https://www.nhm.ac.uk/our-science/services/data/biodiversity-intactness-index.html
# Data Citation: Helen Phillips; Adriana De Palma; Ricardo E Gonzalez; Sara Contu et al. (2021). The Biodiversity Intactness Index - country, region and global-level summaries for the year 1970 to 2050 under various scenarios [Data set]. Natural History Museum. https://doi.org/10.5519/he1eqmg1


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
names(area_mapping)<-c("area_code","country","parent_code")

# add area_names
bii_data_named <- bii_data %>%
  left_join(area_mapping, by = "area_code")


# Choose a few example countries
selected_countries <- c("India", "Indonesia", "Australia")

bii_subset <- bii_data_named %>%
  filter(country %in% selected_countries, variable=="bii")

ggplot(bii_subset, aes(x = year, y = value, color = scenario)) +
  geom_point() +
  facet_wrap(~ country, scales = "free_y") +
  labs(title = "BII Time Series by Country",
       x = "Year", y = "Biodiversity Intactness Index (%)") +
  theme_minimal()

### data are in 5 year blocks

library(readr)

# Export full cleaned BII data with country names
write_csv(bii_data_named, "output/bii_timeseries_by_country.csv")


########## 2. Combine biodiversity & economci indicators

install.packages("gapminder")
library(gapminder)
library(dplyr)
library(ggplot2)
#cite: https://www.gapminder.org/free-material/

# check it
ggplot(filter(gapminder,country=="Australia"), aes(x = year, y =  gdpPercap)) +
  geom_point() +
  geom_point(filter(bii_subset,country=="Australia"), aes(x = year, y =  value)) +
  facet_wrap(~ country, scales = "free_y") +
  labs(title = "GDP per capita  Time Series by Country",
       x = "Year", y = "GDPPP") +
  theme_minimal()

# get decadal averages by country

gapminder %>% summarise(gdpPercap,c(country, continent))

gapminder_decadal <- gapminder %>%
  mutate(decade = floor(year / 10) * 10) %>%
  group_by(decade,country, continent) %>%
  summarise(mean_GDPPP = mean(gdpPercap, na.rm = TRUE))%>%
  rename(year=decade) 


###### merge gapminder with BII data.
bii_data_named$country<-as.factor(bii_data_named$country)
country_list<-levels(gapminder$country) 

bii_gdp<- bii_data_named %>% filter(bii_data_named$country %in% country_list & variable=="bii") %>% 
  right_join(gapminder_decadal,by=c("country","year"))%>%
  filter(year >= 1970)%>%
  select(value,variable,year,country,continent,mean_GDPPP)%>%
  ungroup()

summary(bii_gdp) 


# snapshot to see how compares across countries ina given year- no clear pattern
ggplot(filter(bii_gdp,country=="Australia"), aes(y = value, x = year, color = continent)) +
  geom_point() +
  labs(title = "BII vs time Australia",
       x ="Year" , y = "Biodiversity Intactness Index (%)") +
  theme_minimal()

# change in BII vs change in GDP relative to 1970

bii_gdp <- bii_gdp %>%
  group_by(country) %>%
  arrange(year) %>%
  mutate(change_bii = 100 * (value - first(value)) / first(value)) %>%
  mutate(change_gdp = 100 * (mean_GDPPP - first(mean_GDPPP)) / first(mean_GDPPP)) %>%
  ungroup()


# snapshot to see how compares across countries ina given year- no clear pattern
ggplot(filter(bii_gdp, year == 1990), aes(y = change_gdp, x = change_bii, color = continent)) +
  geom_point() +  # data points
  geom_hline(yintercept = 0, color = "black", linetype = "solid") +  # horizontal line at y = 0
  geom_vline(xintercept = 0, color = "black", linetype = "solid") +  # vertical line at x = 0
  labs(
    title = "Change in BII vs Change in GDP",
    x = "Change BII (%)",
    y = "Change GDP (%)"
  ) +
  theme_minimal()

# Export merged data with country names
write_csv(bii_data_named, "output/bii_gdppp_timeseries_by_country.csv")



