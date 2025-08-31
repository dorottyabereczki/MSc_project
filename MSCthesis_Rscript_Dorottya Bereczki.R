# === Descriptive statistics of explanatory variables ====
# Loading the required package
library(psych)  # for describe()

# Selecting only the relevant variables
vars <- df_clean %>%
  select(DA_price, consumption, RES_generation, net_import,
         temperature, sunshine, wind, rain, TTF_price)

# Generating the descriptive statistics
describe(vars)[, c("mean", "sd", "min", "median", "max", "n")]




# ====== Panel Data Regression: Main model ======

# Step 0: Loading the libraries
# Loading the required libraries
library(dplyr)
library(plm)
library(lmtest)
library(sandwich)
library(broom)
library(car)

# Loading the data
load("~/MSC/df.RData")

# STEP 1: Datetime formatting
df$delivery_start <- as.POSIXct(df$delivery_start, format = "%d.%m.%Y %H:%M:%S")
df$delivery_end   <- as.POSIXct(df$delivery_end, format = "%d.%m.%Y %H:%M:%S")
# this changes the dates to the required format

# STEP 2: Cleaning and creating unique delivery_id because of the DST duplicates
df_clean <- df %>%
  filter(!is.na(DA_price), !is.na(consumption), !is.na(RES_generation),
         !is.na(net_import), !is.na(temperature), !is.na(sunshine),
         !is.na(wind), !is.na(rain)) %>%
  # this filters out the rows with missing data, so the model won't break
  
  mutate(
    day = as.factor(weekdays(delivery_start)),
    month = as.factor(format(delivery_start, "%m"))
  ) %>%
  # this creates the time dummies
  
  group_by(country, delivery_start) %>%
  mutate(duplicate_instance = row_number()) %>%
  ungroup() %>%
  # this identifies the duplicates by country and timestamp
  
  mutate(delivery_id = paste0(country, "_", format(delivery_start, "%Y-%m-%d %H:%M:%S"), "_", duplicate_instance)) %>%
  arrange(country, delivery_id)
# this creates the new unique delivery IDs

# STEP 3: Declaring the panel data structure
pdata <- pdata.frame(df_clean, index = c("country", "delivery_id"))

# STEP 4: Estimating the Fixed Effects model with weather and time dummies
fe_model <- plm(DA_price ~ consumption + RES_generation + net_import + temperature + sunshine + wind + rain + day + month,
                data = pdata, model = "within")
summary(fe_model)

# STEP 5: Hausman test for Fixed vs Random effects
re_model <- plm(DA_price ~ consumption + RES_generation + net_import + temperature + sunshine + wind + rain + day + month,
                data = pdata, model = "random", random.method = "walhus")
cat("\n--- Hausman Test: Fixed vs Random Effects ---\n")
hausman_result <- phtest(fe_model, re_model)
print(hausman_result)

# STEP 7: IV Estimation by instrumenting RES_generation with sunshine and wind
iv_model <- plm(
  DA_price ~ consumption + net_import + temperature + rain + RES_generation + day + month |
    consumption + net_import + temperature + rain + day + month + sunshine + wind,
  data = pdata,
  model = "within"
)
# on the left side of the "|" is the structural equation, and on the right side of it are the instruments

cat("\n--- IV Fixed Effects Model (2SLS) ---\n")
summary(iv_model)

# First stage IV regression to evaluate whether the instruments are strong
first_stage <- plm(
  RES_generation ~ consumption + net_import + temperature + rain + day + month + sunshine + wind,
  data = pdata,
  model = "within"
)

summary(first_stage)$fstatistic

# STEP 8: Checking for multicollinearity in the IV model 
cat("\n--- VIF Check ---\n")
ols_full <- lm(DA_price ~ consumption + RES_generation + net_import + temperature + rain + day + month + sunshine + wind, data = pdata)
print(vif(ols_full))

# STEP 9: Homoskedasticity test (Modified Wald test on IV model)
# ChatGPT4o was used as an aid for this step to make sure it is in line with Baum (2001) and the code is working, but it was adjusted to the model
cat("\n--- Modified Wald Test for Groupwise Heteroskedasticity (IV Model) ---\n")

# Extracting the residuals and number of units
iv_resid <- residuals(iv_model)
# the residuals from the 2SLS IV regression
unit_ids <- pdata$country
# the country each residual belongs to (e.g., DK, NO, SE)

# Creating the data frame with residuals and unit identifiers
resid_df <- data.frame(resid = iv_resid, unit = unit_ids)
# cleaning the structure, each row = one observation, with its residual and associated country

# Computing the σ̂²ᵢ for each group (estimating each country's residual variance sigma2_i)
sigma2_i <- resid_df %>%
  group_by(unit) %>%
  summarise(
    Ti = n(),
    sigma2_hat_i = mean(resid^2),
    .groups = "drop"
  )
# this is the estimated error variance for each country

# Computing the overall σ̂²p (average residual variance across all countries)
sigma2_bar <- mean(sigma2_i$sigma2_hat_i)
# this is the overall average variance to which each country is compared

# Computing the estimated variance of σ̂²ᵢ for each group 
Vi <- resid_df %>%
  group_by(unit) %>%
  summarise(
    Ti = n(),
    sigma2_hat_i = mean(resid^2),
    Vi = var(resid^2) / Ti,
    .groups = "drop"
  )
# this is the estimated variance of each estimated variances and it is needed for the Modified Wald test

# Merging into one table (the variance and the variance of variance tables to one)
wald_df <- merge(sigma2_i, Vi, by = c("unit", "Ti", "sigma2_hat_i"))
# it gives one dataframe

# Computing the test statistic
wald_df$numerator <- (wald_df$sigma2_hat_i - sigma2_bar)^2
wald_df$wald_component <- wald_df$numerator / wald_df$Vi

W <- sum(wald_df$wald_component)
df_wald <- length(unique(unit_ids))
p_value <- pchisq(W, df = df_wald, lower.tail = FALSE)
# this is the formula from Baum (2001)

# Printing the results
cat("Chi-square (", df_wald, ") = ", round(W, 2), ", p-value = ", signif(p_value, 4), "\n", sep = "")


# STEP 10: Wooldridge Test for Serial Correlation for the IV Model
#ChatGPT4o was used as an aid for this step to make sure the code is working and correctly follows Wooldridge (2001), but it was adjusted to the dataset
cat("\n--- Wooldridge Test for Serial Correlation in IV Model ---\n")

# Extracting the panel structure and residuals
panel_df_iv <- as.data.frame(index(iv_model))
# This gives the panel structure
panel_df_iv$resid <- resid(iv_model)
# This gives the residuals from the 2SLS (IV) model

# Adding observation ID and arranging by country and time
panel_df_iv <- panel_df_iv %>%
  mutate(obs_id = row_number()) %>%
  arrange(country, obs_id)
# This makes sure that all observations are sorted correclty within each panel unit and adds observation ID in case it's needed

# Creating the lagged residuals for each country
panel_df_iv$lag_resid <- ave(panel_df_iv$resid, panel_df_iv$country, FUN = function(x) dplyr::lag(x))
# this groups the residuals by country and creates the lagged version of each countries' residuals

# Removing the observations with missing lags
df_resid_iv <- panel_df_iv %>% filter(!is.na(lag_resid))

# Running the Wooldridge regression
wooldridge_iv_lm <- lm(resid ~ lag_resid, data = df_resid_iv)

# Applying the robust standard errors
wooldridge_iv_test <- coeftest(wooldridge_iv_lm, vcov = vcovHC(wooldridge_iv_lm, type = "HC1"))

# Printing the result
print(wooldridge_iv_test)

# STEP 11: Driscoll-Kraay SEs (IV Model)
cat("\n--- Driscoll-Kraay SEs (IV Model) ---\n")
coeftest(iv_model, vcovSCC(iv_model, type = "HC1"))




# ==== ROBUSTNESS TEST 1 ====

# Step 0: Loading the libraries and cleaning the base data again
library(dplyr)
library(plm)
library(lmtest)
library(sandwich)
library(zoo)
library(car)

df_clean_1 <- df %>%
  filter(!is.na(DA_price), !is.na(consumption), !is.na(RES_generation),
         !is.na(net_import), !is.na(temperature), !is.na(sunshine),
         !is.na(wind), !is.na(rain)) %>%
  mutate(
    day = as.factor(weekdays(delivery_start)),
    month = as.factor(format(delivery_start, "%m"))
  ) %>%
  group_by(country, delivery_start) %>%
  mutate(duplicate_instance = row_number()) %>%
  ungroup() %>%
  mutate(delivery_id = paste0(country, "_", format(delivery_start, "%Y-%m-%d %H:%M:%S"), "_", duplicate_instance)) %>%
  arrange(country, delivery_start)

# Step 1: Creating the interaction terms by country (omitting SE to avoid the dummy variable trap)
df1 <- df_clean_1 %>%
  mutate(
    dk = ifelse(country == "DK", 1, 0),
    no = ifelse(country == "NO", 1, 0),
    RES_dk = RES_generation * dk,
    RES_no = RES_generation * no
    # RES_se is omitted → Sweden is reference category
  )

# Step 2: Declaring the panel structure
pdata1 <- pdata.frame(df1, index = c("country", "delivery_id"))

# Step 3: Estimating the IV model with RES_dk and RES_no instrumented by sunshine and wind
iv_model_1 <- plm(
  DA_price ~ consumption + net_import + temperature + rain + RES_dk + RES_no + day + month |
    consumption + net_import + temperature + rain + day + month + dk:sunshine + dk:wind + no:sunshine + no:wind,
  data = pdata1,
  model = "within"
)

cat("\n--- IV Model Summary (RT1, RES × Country instrumented by weather × Country) ---\n")
print(summary(iv_model_1))

# First-stage for Denmark
first_stage_dk <- plm(
  RES_dk ~ consumption + net_import + temperature + rain + day + month + dk:sunshine + dk:wind + no:sunshine + no:wind,
  data = pdata1,
  model = "within"
)
summary(first_stage_dk)$fstatistic

# First-stage for Norway
first_stage_no <- plm(
  RES_no ~ consumption + net_import + temperature + rain + day + month + dk:sunshine + dk:wind + no:sunshine + no:wind,
  data = pdata1,
  model = "within"
)
summary(first_stage_no)$fstatistic

# Step 4: Checking for multicollinearity in the IV model 
cat("\n--- VIF Check ---\n")
print(vif(lm(DA_price ~ consumption + net_import + temperature + sunshine + wind + rain +
               RES_dk + RES_no + day + month, data = df1)))

# Step 5: Homoskedasticity test (Modified Wald test on IV model)
cat("\n--- Modified Wald Test for Groupwise Heteroskedasticity (IV Model, RT1) ---\n")

# Extracting the residuals and unit IDs
iv_resid_1 <- residuals(iv_model_1)
unit_ids_1 <- pdata1$country

# Creating the data frame with residuals and unit identifiers
resid_df_1 <- data.frame(resid = iv_resid_1, unit = unit_ids_1)

# Computing the σ̂²ᵢ for each group
sigma2_i_1 <- resid_df_1 %>%
  group_by(unit) %>%
  summarise(
    Ti = n(),
    sigma2_hat_i = mean(resid^2),
    .groups = "drop"
  )

# Computing the overall σ̂²
sigma2_bar_1 <- mean(sigma2_i_1$sigma2_hat_i)

# Computing the estimated variance of σ̂²ᵢ for each group
Vi_1 <- resid_df_1 %>%
  group_by(unit) %>%
  summarise(
    Ti = n(),
    sigma2_hat_i = mean(resid^2),
    Vi = var(resid^2) / Ti,
    .groups = "drop"
  )

# Merging together the two data frames
wald_df_1 <- merge(sigma2_i_1, Vi_1, by = c("unit", "Ti", "sigma2_hat_i"))

# Computing the Wald statistic
wald_df_1$numerator <- (wald_df_1$sigma2_hat_i - sigma2_bar_1)^2
wald_df_1$wald_component <- wald_df_1$numerator / wald_df_1$Vi

W_1 <- sum(wald_df_1$wald_component)
df_wald_1 <- length(unique(unit_ids_1))
p_value_1 <- pchisq(W_1, df = df_wald_1, lower.tail = FALSE)

# Printing the test result
cat("Chi-square (", df_wald_1, ") = ", round(W_1, 2), ", p-value = ", signif(p_value_1, 4), "\n", sep = "")


# Step 6: Wooldridge Test for Serial Correlation for the IV Model
cat("\n--- Wooldridge Test for Serial Correlation ---\n")

panel_df_1 <- as.data.frame(index(iv_model_1))
panel_df_1$resid <- resid(iv_model_1)

panel_df_1 <- panel_df_1 %>%
  mutate(obs_id = row_number()) %>%
  arrange(country, obs_id)

panel_df_1$lag_resid <- ave(panel_df_1$resid, panel_df_1$country, FUN = function(x) dplyr::lag(x))

df_resid_1 <- panel_df_1 %>%
  filter(!is.na(lag_resid))

wooldridge_lm_1 <- lm(resid ~ lag_resid, data = df_resid_1)
wooldridge_test_1 <- coeftest(wooldridge_lm_1, vcov = vcovHC(wooldridge_lm_1, type = "HC1"))
print(wooldridge_test_1)

# Step 7: Driscoll-Kraay SEs for IV model
cat("\n--- Driscoll-Kraay SE (IV Model) ---\n")
print(coeftest(iv_model_1, vcovSCC(iv_model_1, type = "HC1")))





# ==== ROBUSTNESS TEST 2 ====

# Step 0: Loading the libraries, cleaning the data again and adding time dummies
library(dplyr)
library(plm)
library(lmtest)
library(sandwich)
library(zoo)
library(car)

df_clean_2 <- df %>%
  filter(!is.na(DA_price), !is.na(consumption), !is.na(RES_generation),
         !is.na(net_import), !is.na(temperature), !is.na(sunshine),
         !is.na(wind), !is.na(rain), !is.na(TTF_price)) %>%
  mutate(
    day = as.factor(weekdays(delivery_start)),
    month = as.factor(format(delivery_start, "%m"))
  ) %>%
  group_by(country, delivery_start) %>%
  mutate(duplicate_instance = row_number()) %>%
  ungroup() %>%
  mutate(delivery_id = paste0(country, "_", format(delivery_start, "%Y-%m-%d %H:%M:%S"), "_", duplicate_instance)) %>%
  arrange(country, delivery_start)

# Step 1: Declaring the panel structure
pdata2 <- pdata.frame(df_clean_2, index = c("country", "delivery_id"))

# Step 2: Estimating IV model with RES_generation instrumented by sunshine and wind
iv_model_2 <- plm(
  DA_price ~ consumption + net_import + temperature + rain + TTF_price + RES_generation + day + month |
    consumption + net_import + temperature + rain + TTF_price + day + month + sunshine + wind,
  data = pdata2,
  model = "within"
)

cat("\n--- IV Model Summary (RT2, RES instrumented by sunshine & wind) ---\n")
print(summary(iv_model_2))

# Step 3: Checking for multicollinearity in the IV model 
cat("\n--- VIF Check ---\n")
print(vif(lm(DA_price ~ consumption + RES_generation + net_import + temperature + sunshine +
               wind + rain + TTF_price + day + month, data = df_clean_2)))

# Step 4: Homoskedasticity test (Modified Wald test on IV model)
cat("\n--- Modified Wald Test for Groupwise Heteroskedasticity (IV Model, RT2) ---\n")

# Extracting the residuals and unit IDs
iv_resid_2 <- residuals(iv_model_2)
unit_ids_2 <- pdata2$country

# Creating the data frame with residuals and unit identifiers
resid_df_2 <- data.frame(resid = iv_resid_2, unit = unit_ids_2)

# Computing the σ̂²ᵢ for each group
sigma2_i_2 <- resid_df_2 %>%
  group_by(unit) %>%
  summarise(
    Ti = n(),
    sigma2_hat_i = mean(resid^2),
    .groups = "drop"
  )

# Computing the overall σ̂²
sigma2_bar_2 <- mean(sigma2_i_2$sigma2_hat_i)

# Computing the estimated variance of σ̂²ᵢ for each group
Vi_2 <- resid_df_2 %>%
  group_by(unit) %>%
  summarise(
    Ti = n(),
    sigma2_hat_i = mean(resid^2),
    Vi = var(resid^2) / Ti,
    .groups = "drop"
  )

# Merging the the two data frames
wald_df_2 <- merge(sigma2_i_2, Vi_2, by = c("unit", "Ti", "sigma2_hat_i"))

# Computing the Wald statistic
wald_df_2$numerator <- (wald_df_2$sigma2_hat_i - sigma2_bar_2)^2
wald_df_2$wald_component <- wald_df_2$numerator / wald_df_2$Vi

W_2 <- sum(wald_df_2$wald_component)
df_wald_2 <- length(unique(unit_ids_2))
p_value_2 <- pchisq(W_2, df = df_wald_2, lower.tail = FALSE)

# Printing the test result
cat("Chi-square (", df_wald_2, ") = ", round(W_2, 2), ", p-value = ", signif(p_value_2, 4), "\n", sep = "")


# Step 5: Wooldridge Test for Serial Correlation for the IV Model
cat("\n--- Wooldridge Test for Serial Correlation ---\n")

panel_df_2 <- as.data.frame(index(iv_model_2))
panel_df_2$resid <- resid(iv_model_2)

panel_df_2 <- panel_df_2 %>%
  mutate(obs_id = row_number()) %>%
  arrange(country, obs_id)

panel_df_2$lag_resid <- ave(panel_df_2$resid, panel_df_2$country, FUN = function(x) dplyr::lag(x))

df_resid_2 <- panel_df_2 %>%
  filter(!is.na(lag_resid))

wooldridge_lm_2 <- lm(resid ~ lag_resid, data = df_resid_2)
wooldridge_test_2 <- coeftest(wooldridge_lm_2, vcov = vcovHC(wooldridge_lm_2, type = "HC1"))
print(wooldridge_test_2)

# Step 6: Driscoll-Kraay SE (IV Model)
cat("\n--- Driscoll-Kraay SE (IV Model) ---\n")
print(coeftest(iv_model_2, vcovSCC(iv_model_2, type = "HC1")))





# ==== ROBUSTNESS TEST 3 ====

# Step 0: Loading the libraries, cleaning the data and adding time dummies again
library(dplyr)
library(plm)
library(lmtest)
library(sandwich)
library(zoo)
library(car)

df_clean_3 <- df %>%
  filter(!is.na(DA_price), !is.na(consumption), !is.na(RES_generation),
         !is.na(net_import), !is.na(temperature), !is.na(sunshine),
         !is.na(wind), !is.na(rain), !is.na(TTF_price)) %>%
  mutate(
    day = as.factor(weekdays(delivery_start)),
    month = as.factor(format(delivery_start, "%m"))
  ) %>%
  group_by(country, delivery_start) %>%
  mutate(duplicate_instance = row_number()) %>%
  ungroup() %>%
  mutate(delivery_id = paste0(country, "_", format(delivery_start, "%Y-%m-%d %H:%M:%S"), "_", duplicate_instance)) %>%
  arrange(country, delivery_start)

# Step 1: Creating the interaction terms between TTF_price and country (omitting SE to avoid the dummy variable trap)
df3 <- df_clean_3 %>%
  mutate(
    dk = ifelse(country == "DK", 1, 0),
    no = ifelse(country == "NO", 1, 0),
    TTF_dk = TTF_price * dk,
    TTF_no = TTF_price * no
  )

# Step 2: Declaring the panel structure
pdata3 <- pdata.frame(df3, index = c("country", "delivery_id"))

# Step 3: Estimating the IV model with RES_generation instrumented by sunshine and wind
iv_model_3 <- plm(
  DA_price ~ consumption + net_import + temperature + rain + RES_generation + TTF_dk + TTF_no + day + month |
    consumption + net_import + temperature + rain + TTF_dk + TTF_no + day + month + sunshine + wind,
  data = pdata3,
  model = "within"
)

cat("\n--- IV Model Summary (RT3, RES instrumented by sunshine & wind) ---\n")
print(summary(iv_model_3))

# Step 4: Checking for multicollinearity in the IV model 
cat("\n--- VIF Check ---\n")
print(vif(lm(DA_price ~ consumption + RES_generation + net_import + temperature + sunshine +
               wind + rain + TTF_dk + TTF_no + day + month, data = df3)))

# Step 5: Homoskedasticity test (Modified Wald test on IV model)
cat("\n--- Modified Wald Test for Groupwise Heteroskedasticity (IV Model, RT3) ---\n")

# Extracting the residuals and unit IDs
iv_resid_3 <- residuals(iv_model_3)
unit_ids_3 <- pdata3$country

# Creating the data frame with residuals and unit identifiers
resid_df_3 <- data.frame(resid = iv_resid_3, unit = unit_ids_3)

# Computting the σ̂²ᵢ for each group
sigma2_i_3 <- resid_df_3 %>%
  group_by(unit) %>%
  summarise(
    Ti = n(),
    sigma2_hat_i = mean(resid^2),
    .groups = "drop"
  )

# Computing the overall σ̂²
sigma2_bar_3 <- mean(sigma2_i_3$sigma2_hat_i)

# Computing the estimated variance of σ̂²ᵢ for each group
Vi_3 <- resid_df_3 %>%
  group_by(unit) %>%
  summarise(
    Ti = n(),
    sigma2_hat_i = mean(resid^2),
    Vi = var(resid^2) / Ti,
    .groups = "drop"
  )

# Merge=ing the two data frames
wald_df_3 <- merge(sigma2_i_3, Vi_3, by = c("unit", "Ti", "sigma2_hat_i"))

# Computing the Wald statistic
wald_df_3$numerator <- (wald_df_3$sigma2_hat_i - sigma2_bar_3)^2
wald_df_3$wald_component <- wald_df_3$numerator / wald_df_3$Vi

W_3 <- sum(wald_df_3$wald_component)
df_wald_3 <- length(unique(unit_ids_3))
p_value_3 <- pchisq(W_3, df = df_wald_3, lower.tail = FALSE)

# Printing the test result
cat("Chi-square (", df_wald_3, ") = ", round(W_3, 2), ", p-value = ", signif(p_value_3, 4), "\n", sep = "")


# Step 6: Wooldridge Test for Serial Correlation for the IV Model
cat("\n--- Wooldridge Test for Serial Correlation  ---\n")

panel_df_3 <- as.data.frame(index(iv_model_3))
panel_df_3$resid <- resid(iv_model_3)

panel_df_3 <- panel_df_3 %>%
  mutate(obs_id = row_number()) %>%
  arrange(country, obs_id)

panel_df_3$lag_resid <- ave(panel_df_3$resid, panel_df_3$country, FUN = function(x) dplyr::lag(x))

df_resid_3 <- panel_df_3 %>%
  filter(!is.na(lag_resid))

wooldridge_lm_3 <- lm(resid ~ lag_resid, data = df_resid_3)
wooldridge_test_3 <- coeftest(wooldridge_lm_3, vcov = vcovHC(wooldridge_lm_3, type = "HC1"))
print(wooldridge_test_3)

# Step 7: Driscoll-Kraay SEs for IV model
cat("\n--- Driscoll-Kraay SE (IV Model) ---\n")
print(coeftest(iv_model_3, vcovSCC(iv_model_3, type = "HC1")))





# ==== ROBUSTNESS TEST 4 ====

# Step 0: Loading the libraries 
library(dplyr)
library(plm)
library(lmtest)
library(sandwich)
library(zoo)
library(car)

# Step 1: Cleaning the data and adding time dummies
df_clean_6 <- df %>%
  filter(!is.na(DA_price), !is.na(consumption), !is.na(RES_generation),
         !is.na(net_import), !is.na(temperature), !is.na(sunshine),
         !is.na(wind), !is.na(rain)) %>%
  mutate(
    day = relevel(factor(weekdays(delivery_start)), ref = "Friday"),
    month = relevel(factor(format(delivery_start, "%m")), ref = "01")
  ) %>%
  group_by(country, delivery_start) %>%
  mutate(duplicate_instance = row_number()) %>%
  ungroup() %>%
  mutate(delivery_id = paste0(country, "_", format(delivery_start, "%Y-%m-%d %H:%M:%S"), "_", duplicate_instance)) %>%
  arrange(country, delivery_start)

# Step 2: Creating the dummy for high RES production days (top 5%)
quantile_res <- quantile(df_clean_6$RES_generation, 0.95, na.rm = TRUE)
df5 <- df_clean_6 %>%
  mutate(high_res_dummy = ifelse(RES_generation >= quantile_res, 1, 0))

# Step 3: Declaring the panel structure
pdata5 <- pdata.frame(df5, index = c("country", "delivery_id"))

# Step 4: Estimating the IV model with RES_generation instrumented by sunshine and wind
iv_model_5 <- plm(
  DA_price ~ consumption + net_import + temperature + rain + RES_generation + high_res_dummy + day + month |
    consumption + net_import + temperature + rain + high_res_dummy + day + month + sunshine + wind,
  data = pdata5,
  model = "within"
)

cat("\n--- IV Model Summary (RT4, RES instrumented by sunshine & wind) ---\n")
print(summary(iv_model_5))

# Step 5: Checking for multicollinearity in the IV model 
cat("\n--- Step 5: VIF Check ---\n")
print(vif(lm(DA_price ~ consumption + RES_generation + net_import + temperature + sunshine +
               wind + rain + high_res_dummy + day + month, data = df5)))

# Step 6: Homoskedasticity test (Modified Wald test on IV model)
cat("\n--- Step 6: Modified Wald Test for Groupwise Heteroskedasticity (IV Model, RT4) ---\n")

# Extracting the residuals and unit IDs
iv_resid_5 <- residuals(iv_model_5)
unit_ids_5 <- pdata5$country

# Creating the data frame with residuals and unit identifiers
resid_df_5 <- data.frame(resid = iv_resid_5, unit = unit_ids_5)

# Computing the σ̂²ᵢ for each group
sigma2_i_5 <- resid_df_5 %>%
  group_by(unit) %>%
  summarise(
    Ti = n(),
    sigma2_hat_i = mean(resid^2),
    .groups = "drop"
  )

# Computing the overall σ̂²
sigma2_bar_5 <- mean(sigma2_i_5$sigma2_hat_i)

# Computing the estimated variance of σ̂²ᵢ for each group
Vi_5 <- resid_df_5 %>%
  group_by(unit) %>%
  summarise(
    Ti = n(),
    sigma2_hat_i = mean(resid^2),
    Vi = var(resid^2) / Ti,
    .groups = "drop"
  )

# Merging the the two data frames
wald_df_5 <- merge(sigma2_i_5, Vi_5, by = c("unit", "Ti", "sigma2_hat_i"))

# Computing the Wald statistic
wald_df_5$numerator <- (wald_df_5$sigma2_hat_i - sigma2_bar_5)^2
wald_df_5$wald_component <- wald_df_5$numerator / wald_df_5$Vi

W_5 <- sum(wald_df_5$wald_component)
df_wald_5 <- length(unique(unit_ids_5))
p_value_5 <- pchisq(W_5, df = df_wald_5, lower.tail = FALSE)

# Printing the test result
cat("Chi-square (", df_wald_5, ") = ", round(W_5, 2), ", p-value = ", signif(p_value_5, 4), "\n", sep = "")


# Step 7: Wooldridge Test for Serial Correlation for the IV Model
cat("\n--- Step 7: Wooldridge Test for Serial Correlation ---\n")

panel_df_5 <- as.data.frame(index(iv_model_5))
panel_df_5$resid <- resid(iv_model_5)

panel_df_5 <- panel_df_5 %>%
  mutate(obs_id = row_number()) %>%
  arrange(country, obs_id)

panel_df_5$lag_resid <- ave(panel_df_5$resid, panel_df_5$country, FUN = function(x) dplyr::lag(x))

df_resid_5 <- panel_df_5 %>%
  filter(!is.na(lag_resid))

wooldridge_lm_5 <- lm(resid ~ lag_resid, data = df_resid_5)
wooldridge_test_5 <- coeftest(wooldridge_lm_5, vcov = vcovHC(wooldridge_lm_5, type = "HC1"))
print(wooldridge_test_5)

# Step 8: Driscoll-Kraay robust standard errors for IV model
cat("\n--- Step 8: Driscoll-Kraay SE (IV Model) ---\n")
print(coeftest(iv_model_5, vcovSCC(iv_model_5, type = "HC1")))





# === ROBUSTNESS TEST 5: High RES dummy × Country ===

# Step 0: Loading the needed libraries
library(dplyr)
library(plm)
library(lmtest)
library(sandwich)
library(zoo)
library(car)

# Step 1: Cleaning the data and adding time dummies
df_clean_6 <- df %>%
  filter(!is.na(DA_price), !is.na(consumption), !is.na(RES_generation),
         !is.na(net_import), !is.na(temperature), !is.na(sunshine),
         !is.na(wind), !is.na(rain)) %>%
  mutate(
    day = relevel(factor(weekdays(delivery_start)), ref = "Friday"),
    month = relevel(factor(format(delivery_start, "%m")), ref = "01")
  ) %>%
  group_by(country, delivery_start) %>%
  mutate(duplicate_instance = row_number()) %>%
  ungroup() %>%
  mutate(
    delivery_id = paste0(country, "_", format(delivery_start, "%Y-%m-%d %H:%M:%S"), "_", duplicate_instance)
  ) %>%
  arrange(country, delivery_start)

# Step 2: Creating the high RES dummy for the countries (top 5%) and interaction terms
quantiles_by_country <- df_clean_6 %>%
  group_by(country) %>%
  summarise(res_95 = quantile(RES_generation, 0.95, na.rm = TRUE)) %>%
  ungroup()

df6 <- df_clean_6 %>%
  left_join(quantiles_by_country, by = "country") %>%
  # this merges the quantile data back to the main dataframe with matching the countries
  mutate(
    dk = ifelse(country == "DK", 1, 0),
    no = ifelse(country == "NO", 1, 0),
    se = ifelse(country == "SE", 1, 0),
    
    high_res_dk_dummy = ifelse(country == "DK" & RES_generation >= res_95, 1, 0),
    high_res_no_dummy = ifelse(country == "NO" & RES_generation >= res_95, 1, 0),
    high_res_se_dummy = ifelse(country == "SE" & RES_generation >= res_95, 1, 0),
    
    high_res_dk = high_res_dk_dummy * dk,
    high_res_no = high_res_no_dummy * no,
    high_res_se = high_res_se_dummy * se
  )

# Step 3: Declaring the panel structure
pdata6 <- pdata.frame(df6, index = c("country", "delivery_id"))

# Step 4: Estimating the IV model with RES_generation instrumented by sunshine and wind
iv_model_6 <- plm(
  DA_price ~ consumption + net_import + temperature + rain + RES_generation +
    high_res_dk + high_res_no + day + month |
    consumption + net_import + temperature + rain +
    high_res_dk + high_res_no + day + month + sunshine + wind,
  data = pdata6,
  model = "within"
)

cat("\n--- IV Model Summary (RT6, RES instrumented by sunshine & wind) ---\n")
print(summary(iv_model_6))

# Step 5: Checking for multicollinearity in the IV model 
cat("\n--- Step 5: VIF Check ---\n")
print(vif(lm(DA_price ~ consumption + RES_generation + net_import + temperature + sunshine +
               wind + rain + day + month + high_res_dk + high_res_no, data = df6)))

# Step 6: Homoskedasticity test (Modified Wald test on IV model)
cat("\n--- Step 6: Modified Wald Test for Groupwise Heteroskedasticity (IV Model, RT6) ---\n")

# Extracting the residuals and unit IDs
iv_resid_6 <- residuals(iv_model_6)
unit_ids_6 <- pdata6$country

# Creating the data frame with residuals and unit identifiers
resid_df_6 <- data.frame(resid = iv_resid_6, unit = unit_ids_6)

# Computing the σ̂²ᵢ for each group
sigma2_i_6 <- resid_df_6 %>%
  group_by(unit) %>%
  summarise(
    Ti = n(),
    sigma2_hat_i = mean(resid^2),
    .groups = "drop"
  )

# Computing the overall σ̂²
sigma2_bar_6 <- mean(sigma2_i_6$sigma2_hat_i)

# Computing the estimated variance of σ̂²ᵢ for each group
Vi_6 <- resid_df_6 %>%
  group_by(unit) %>%
  summarise(
    Ti = n(),
    sigma2_hat_i = mean(resid^2),
    Vi = var(resid^2) / Ti,
    .groups = "drop"
  )

# Merging the the two data frames
wald_df_6 <- merge(sigma2_i_6, Vi_6, by = c("unit", "Ti", "sigma2_hat_i"))

# Computing the Wald statistic
wald_df_6$numerator <- (wald_df_6$sigma2_hat_i - sigma2_bar_6)^2
wald_df_6$wald_component <- wald_df_6$numerator / wald_df_6$Vi

W_6 <- sum(wald_df_6$wald_component)
df_wald_6 <- length(unique(unit_ids_6))
p_value_6 <- pchisq(W_6, df = df_wald_6, lower.tail = FALSE)

# Printing the test result
cat("Chi-square (", df_wald_6, ") = ", round(W_6, 2), ", p-value = ", signif(p_value_6, 4), "\n", sep = "")


# Step 7: Wooldridge Test for Serial Correlation for the IV Model
cat("\n--- Step 7: Wooldridge Test for Serial Correlation  ---\n")

panel_df_6 <- as.data.frame(index(iv_model_6))
panel_df_6$resid <- resid(iv_model_6)

panel_df_6 <- panel_df_6 %>%
  mutate(obs_id = row_number()) %>%
  arrange(country, obs_id)

panel_df_6$lag_resid <- ave(panel_df_6$resid, panel_df_6$country, FUN = function(x) dplyr::lag(x))

df_resid_6 <- panel_df_6 %>%
  filter(!is.na(lag_resid))

wooldridge_lm_6 <- lm(resid ~ lag_resid, data = df_resid_6)
wooldridge_test_6 <- coeftest(wooldridge_lm_6, vcov = vcovHC(wooldridge_lm_6, type = "HC1"))
print(wooldridge_test_6)

# Step 8: Driscoll-Kraay Robust SEs
cat("\n--- Step 8: Driscoll-Kraay SE (RT6 IV) ---\n")
print(coeftest(iv_model_6, vcovSCC(iv_model_6, type = "HC1")))





# ==== Random Effects Model ====

#Step 0:Loading the libraries and cleaning the data
# Loading the required libraries
library(dplyr)
library(plm)
library(lmtest)
library(sandwich)
library(car)
library(multiwayvcov)

# Loading the cleaned data
load("~/MSC/df.RData")

# Step 1: Making sure that the datetime formatting is correct and adding the time dummies
df$delivery_start <- as.POSIXct(df$delivery_start, format = "%d.%m.%Y %H:%M:%S")
df$delivery_end   <- as.POSIXct(df$delivery_end, format = "%d.%m.%Y %H:%M:%S")

# Adding the day and month dummies 
df$day <- as.factor(weekdays(df$delivery_start))
df$month <- as.factor(format(df$delivery_start, "%m"))

# Step 2: Cleaning and creating delivery IDs, which is needed due to the DST duplications
df_clean <- df %>%
  filter(!is.na(DA_price), !is.na(consumption), !is.na(RES_generation),
         !is.na(net_import), !is.na(temperature), !is.na(sunshine)) %>%
  group_by(country, delivery_start) %>%
  mutate(duplicate_instance = row_number()) %>%
  ungroup() %>%
  mutate(delivery_id = paste0(format(delivery_start, "%Y-%m-%d %H:%M:%S"), "_", duplicate_instance)) %>%
  arrange(country, delivery_id)

# Step 3: Declaring the panel data structure
pdata <- pdata.frame(df_clean, index = c("country", "delivery_id"))

# STEP 4: Estimating the Fixed Effects model
fe_model <- plm(DA_price ~ consumption + RES_generation + net_import + temperature +
                  sunshine + wind + rain + day + month,
                data = pdata, model = "within")
summary(fe_model)

# STEP 5: Estimating the Random Effects model
re_model <- plm(DA_price ~ consumption + RES_generation + net_import + temperature +
                  sunshine + wind + rain + day + month,
                data = pdata, model = "random", random.method = "walhus")
summary(re_model)

# STEP 6: Hausman Test
cat("\n--- Hausman Test: Fixed vs Random Effects ---\n")
hausman_test <- phtest(fe_model, re_model)
print(hausman_test
      
# STEP 7: Durbin-Watson test for serial correlation
cat("\n--- Durbin-Watson Test for Serial Correlation (Panel Data) ---\n")
dw_test <- pdwtest(re_model)
print(dw_test)
      
# STEP 8: Multicollinearity - VIF
cat("\n--- VIF check ---\n")
ols_model_re <- lm(DA_price ~ consumption + RES_generation + net_import + temperature +
                           sunshine + wind + rain + day + month, data = pdata)
print(vif(ols_model_re))
      
# STEP 9: Breusch-Pagan test for heteroskedasticity 
# ChatGPT4o was used as an aid for this as there is no R function, but it was manually adjusted for the model and dataset
      
# Extracting the residuals and fitted values
df_clean$resid_re <- resid(re_model)
df_clean$fitted_re <- fitted(re_model)

# Squaring the residuals
df_clean$resid_sq <- df_clean$resid_re^2
      
# Auxiliary regression: squared residuals ~ regressors (or fitted values)
aux_model <- lm(resid_sq ~ consumption + RES_generation + net_import + temperature +
                        sunshine + wind + rain + day + month,
                      data = df_clean)
      
# Computing the test statistic
R2_aux <- summary(aux_model)$r.squared
n_obs <- nrow(df_clean)
df_aux <- length(coef(aux_model)) - 1  # exclude intercept
LM_stat <- n_obs * R2_aux
pval_bp <- 1 - pchisq(LM_stat, df_aux)
      
# Output
cat("LM statistic =", round(LM_stat, 3), "\n")
cat("Degrees of freedom =", df_aux, "\n")
cat("p-value =", round(pval_bp, 4), "\n")
if (pval_bp < 0.05) {
cat("Result: Heteroskedasticity likely present (reject H0)\n")
} else {
cat("Result: No strong evidence of heteroskedasticity (fail to reject H0)\n")
}
      
# STEP 10: Clustered Standard Errors
cat("\n--- Clustered Standard Errors ---\n")
lm_model <- lm(DA_price ~ consumption + RES_generation + net_import + temperature +
                       sunshine + wind + rain + day + month, data = pdata)
clustered_vcov <- cluster.vcov(lm_model, pdata$country)
coeftest(lm_model, clustered_vcov)
      