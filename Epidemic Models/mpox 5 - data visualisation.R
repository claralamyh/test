library(tidyverse)
library(lubridate)

#WHO
who <- read.csv("mpox\\WHO mpox.csv")
Idata_uk <- who %>%
  filter(iso3 == "GBR") %>%
  select(month_lab, cases) %>%
  group_by(month_lab) %>%
  summarise(total_cases = sum(cases)) %>%
  mutate(month_lab = my(month_lab)) %>%
  arrange(month_lab)

yt_uk <- as.vector(Idata_uk$total_cases)
plot(yt_uk, ylab="Cases")


#UKHSA
ukhsa <- read.csv("mpox\\UKHSA mpox.csv")
plot(ukhsa[,2],ukhsa[,3], type="o", xlab="Day", ylab="Cases")


# not comparable, WHO only includes confirmed cases
# while UKHSA also includes highly probable cases
# but we can gain a brief idea on the disease development
# they both show a very similar trend

