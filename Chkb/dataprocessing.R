#################################################################
# - Research Project:
#   Interconnectedness, Performance, Risk, and Response 
#   to Shocks: Evidence from the Banking Network in Chile.
# - Paper:
#   Chanci, Kumbhakar and Bobadilla - 
#   Interconnectedness and Banking Performance
# - File:
#   dataprocessing.R
#   This code is for processing the raw (CSV) data and 
#   creating the dataframes: df.Rda and all_W_df.Rda 
#   Raw data to use (located at Chkb/Data/Raw):
#   - CMF: 
#     raw data C18_Muestra.CSV, 
#     Empl_muestra.CSV, MB2_Muestra.CSV 
#     The staff at the CMF would need to change 
#     or update the content in those files with the full sample.
#   - Central Bank of Chile:
#     IPC_06072024162200.xlsx with the CPI,
#     which was downloaded from the Central Bank's website.
# - By: 
#   Luis Chanci, 
#   Contact: luischanci@santotomas.cl
#   2025
####################################################

# Inital setup (checking for required packages)

required_pkg.R <- c("tidyverse", "readxl", "igraph")
new_packages   <- required_pkg.R[!(required_pkg.R %in% installed.packages()[,"Package"])]

if(length(new_packages)) install.packages(new_packages)
lapply(required_pkg.R, library, character.only = TRUE)


####################################
# Part I. 
# - Balance-Sheet info (form MB2) as a Panel Dataframe 
# - Datasets: 
#   * It first creates the a (raw) dataset MB2
#   * Then, it creates the panel dataframe df
#     which contains the main variables to use in the estimation.

#####
# I.1. MB2 organized by columns

# The location of the folder wd is set in the main.jl file

MB2_Muestra  <- read.csv(paste0(wd,"/Data/Raw/MB2_Muestra.csv"))
MB2_Muestra  <- MB2_Muestra[order(MB2_Muestra$INS_COD, MB2_Muestra$PERIODO, MB2_Muestra$CODIGO),]
       
Empl_muestra <- read.csv(paste0(wd,"/Data/Raw/Empl_muestra.csv"))
Empl_muestra <- Empl_muestra|> 
                  rename(PERIODO = I06_PERIODO)

MB2_Muestra  <- left_join(MB2_Muestra, Empl_muestra, by=c('INS_COD'='INS_COD', 'PERIODO'='PERIODO')) # Adding Labor

# Archivo C04: Activo ponderado por riesgo consolidado (APC), Patrimonio efectivo consolidado (PEC), Indicador reglamentario de solvencia (IRS)
APR_Muestra  <- read.csv(paste0(wd,"/Data/Raw/APR_Muestra.csv"))
APR_Muestra  <- APR_Muestra|> 
                  rename(PERIODO = C04_PERIODO)

MB2_Muestra  <- left_join(MB2_Muestra, APR_Muestra, by=c('INS_COD'='INS_COD', 'PERIODO'='PERIODO')) # Adding Labor

# Main codes to use and Panel setup: 
cod <- sort(unique(MB2_Muestra$CODIGO))

MB2 <- MB2_Muestra|>
          filter(CODIGO %in% cod[1])|>
            mutate(MONTO1000 = MONTO_TOTAL/1000)|>              # values in thousands of CLP
              select(INS_COD, PERIODO, E_Perm, E_Ext, starts_with("C04"), MONTO1000)|>
                rename(!!cod[1] := MONTO1000)

for (ii in 2:length(cod)){ 
  MB2 <-left_join(
         MB2,
          MB2_Muestra|>
           filter(CODIGO %in% cod[ii])|>
            mutate(MONTO1000 = MONTO_TOTAL/1000)|>          
             select(INS_COD, PERIODO, MONTO1000) |>
              rename(!!cod[ii] := MONTO1000), 
        by=c('INS_COD'='INS_COD', 'PERIODO'='PERIODO'))
}

# Inflation-adjusted values (Real values):
IPC <- read_excel(paste0(wd,"/Data/Raw/IPC_06072024162200.xlsx"), range = "A3:D195")|> 
        rename(IPC = starts_with('1.IPC'))|> 
          filter(PERIOD>=201604 & PERIOD<=201712)|>
            select(PERIOD, IPC)|>
              mutate(IPC        = IPC/100, 
                     time_trend = row_number())|>
                rename(PERIODO = PERIOD)

MB2 <- left_join(MB2, IPC, by=c('PERIODO'='PERIODO')) 


MB2 <- MB2 |>
        mutate(across(all_of(cod), ~.x / IPC, .names = "{.col}_r"))

MB2 <- MB2|>
        mutate( date = as.Date(paste0(PERIODO,'01'), format='%Y%m%d') )


# Final Adjustments in the MB2:
MB2 <- MB2|> 
        filter(!is.na(E_Perm) &!is.na(E_Ext)) # Keep INS with full info in the labor variable (later can explore using the total value)


# Notes:
# The following step should be used with caution, as it eliminates institutions with few observations
# within the provided time frame (2016-2017). If the sample is updated, skip this step. 
# Alternatively, you can filter the data later to retain only
# the 20 institutions used in the final sample.

tokeep  <-   c("OR","PP","PS","PU","PW","QN","RO","RR","SN","SP","SR","SV","TP",
               "TR","TT","TV","UP","UQ","SOU","ST")
MB2     <- MB2|> filter(INS_COD%in%tokeep)

# If the previous steps ran correctly, the result is a panel dataset
# containing balance sheet information in real thousands of Chilean Pesos (CLP).


#####
# I.2. Creating the panel (dataframe) df with the variables to be used in the estimation

# Note: This section defines variables using two approaches based on available accounting codes.
# The active commands use the currently available codes. 

# -----------------------------
# Inputs X and Prices P
# -----------------------------

# I.2.1. Input: Labor -----------------------------------------------

  X_labor <- rowSums(MB2[,c(
             "E_Perm", "E_Ext")])                 # Use both permanent and external headcount (available)

  P_labor <- MB2[, "C_4600000_r"]                 # price of labor: Total personnel expenses (a single-line proxy)

# I.2.2. Input: Physical Capital ------------------------------------

  X_capital  <- MB2[, "C_1600000_r"]                          # Net capital 
  P_capital  <- rowSums(MB2[,c("C_4610101_r","C_4620001_r")]) # maintenance/taxes + depreciation/amortization

# I.2.3. Inputs: Deposits  --------------------------

  X_deposit <- rowSums(MB2[,c(
               "C_2100000_r",         # demand/sight deposits (aggregate)
               "C_2200000_r")],       # time/savings deposits (aggregate)
                na.rm = TRUE)
  P_deposit <- rowSums(MB2[,c(
               "C_4150101_r",         # interest on current accounts
               "C_4150102_r",         # interest on other sight obligations
               "C_4150301_r","C_4150302_r","C_4150309_r")])  # interest on time/savings


# I.2.4. Definition of Qusi-fixed input: Equity
  X_equity  <- MB2[, "C_5001000_r"]   
  

# I.2.5.Definition of Input: Risk
  risk      <- MB2[, "C04_IRS"]/100 # A proxy for risk: Capital Adequacy Ratio, a widely used inverse proxy for bank risk
                                             

# -----------------------------
# Outputs Y
# -----------------------------

# I.2.7.Definition of Output: Y1
  Y_cloan  <- MB2[, 'C_1302000_r'] # Commercial (industrial) loans

# I.2.8.Definition of Output: Y2
  Y_rsloan <- MB2[, 'C_1304000_r'] # Real State Loans or Mortgages

# I.2.9.Definition of Output: Y3
  Y_hhloan <- MB2[, 'C_1305000_r'] # Consumer or Household Loans

# I.2.10.Definition of Output: Y4
  Y_secur  <- rowSums(MB2[,c('C_1350000_r', 'C_1360000_r')], na.rm = TRUE) # Securities or any other income generating activities


# -----------------------------
# Dependent Variable: Total Variable Cost (TVC)
# -----------------------------

# --- Preferred (using available codes: matches the constructed input prices; Exclude depreciation/maintenance provisions):
  TVC <- abs(rowSums(MB2[, c(
           "C_4600000_r",                               # labor cost
           "C_4150101_r","C_4150102_r",                 # sight deposit interest expenses
           "C_4150301_r","C_4150302_r","C_4150309_r",   # time/savings interest
           "C_4150400_r","C_4150500_r","C_4150600_r")], # borrowings/bonds/other funds
               na.rm = TRUE))

# -----------------------------
# Assemble df
# -----------------------------

df <- data.frame(
      INS_COD = MB2$INS_COD,  PERIODO = MB2$PERIODO,  time_trend = MB2$time_trend,
      Cost      = TVC,       X_equity  = X_equity, risk     = risk,
      Y_cloan   = Y_cloan,   Y_rsloan  = Y_rsloan, Y_hhloan = Y_hhloan, Y_secur  = Y_secur,
      X_labor   = X_labor,   P_labor   = P_labor,
      X_capital = X_capital, P_capital = P_capital,
      X_deposit = X_deposit, P_deposit = P_deposit)

inp <- c("labor", "capital", "deposit")
df  <- df |>
        mutate(
          purrr::map_dfc(inp, ~
          tibble(
            "{.x}" := {
              X <- get(paste0("X_", .x))
              P <- get(paste0("P_", .x))
              ifelse(!is.na(X) & X > 0, abs(P) / X, NA_real_)
            }
          ))|>
          set_names(paste0("p_", inp))
        )|>
        group_by(INS_COD)|>                        # 'Imputing' the p50 value to few observations with missing input prices
        mutate(across(all_of(paste0("p_", inp)), ~ ifelse(is.na(.x), median(.x, na.rm = TRUE), .x)))|>
        ungroup()                                  

# Extra: 
# For the (complete) model using determinants of mean efficiency,
# - Ownership: Add dummy variables for any bank with a large portion of loans with
#   State Guarantees (e.g., 1302.6.01 , 1302.6.02) or funding from the Public 
#   Sector (e.g., 2500.1) or with many CuentaRUT (e.g., 2100.2.04)
# - Segment: A dummy variable for those with overwhelming concentration in 
#   consumer assets (e.g., 1305 , 1305.4) and lack or  minimal involvement in 
#   corporate lending (e.g., 1302) or unique liability account (e.g., 2500.2.01)
# - Foreign: Add dummy variables for those banks with a high probability of being foreign
#   liabilities to foreign banks (e.g., 2302), loans to foreign entities (e.g., 1270.2 , 1302.1.02),
#   large net foreign exchange accounts (e.g., abs(4350)).

# However, I could only set the following (which may need to be checked later with the CMF):

  df <- df|>
        mutate(d_ownership = ifelse(INS_COD == "PS", 1, 0),
               d_segment   = ifelse(INS_COD %in% c("TT", "TR"), 1, 0))


# Note: If the previous steps executed correctly, 'df' is a panel dataset
# containing the core balance sheet information. All monetary values are expressed
# in real thousands of Chilean Pesos (CLP) and are ready for estimation.
# The dataframe is now saved for use in the estimation algorithm.

save(df, file = paste0(wd,"/Data/Processed/df.Rda")) 



####################################
# Part II. 
# Interconnectedness (Adjacency Matrix) 
# - This part is for organizing the 
#   interbank transactions (the C18)
#   creating time-varying adjacency matrices.


# Load data
C18_Muestra <- read.csv(paste0(wd, "/Data/Raw/C18_Muestra.csv"))

# Initial cleaning and total obligations construction
C18_Muestra <- C18_Muestra|> 
  rename(LENDER = C18_COD_IFI_ACREED, 
         PERIODO= C18_PERIODO)|>
  mutate(across(all_of(
    c("C18_CTAS_CTES", 
      "C18_OT_OBLIG_VISTA", 
      "C18_CONTR_REC_PTMOS",
      "C18_DEP_CAP", 
      "C18_CONTR_DER_FIN", 
      "C18_OBLIG_BCOS")), as.numeric),
    TOT_OBLIG = (C18_CTAS_CTES + C18_OT_OBLIG_VISTA + 
                   C18_CONTR_REC_PTMOS + C18_DEP_CAP + 
                   C18_CONTR_DER_FIN + C18_OBLIG_BCOS),
    UNSEC_OBL = TOT_OBLIG - C18_MTO_GAR )

# Checking (again) to filter for main banks
tokeep <- c("OR","PP","PS","PU","PW","QN","RO","RR","SN","SP","SR","SV","TP",
            "TR","TT","TV","UP","UQ","SOU","ST")

C18_Muestra <- C18_Muestra |> 
                filter(INS_COD %in% tokeep, LENDER %in% tokeep)


# Function fn_W_df: Construct Adjacency Matrix (W)
# Note: A dataframe-centric approach is used to carefully manage
#       and align borrower and lender IDs. This is crucial for
#       constructing a correctly specified matrix for each time period.

fn_W_df <- function(
    data, 
    maturity_value = 2,                     # Options (maturity): 1,2,3
    instrument_var = "C18_CONTR_DER_FIN") { # C18_CONTR_DER_FIN as default
  
  df.W <- data |>
    filter(C18_PZO_RES_VENC == maturity_value) |>
    select(INS_COD, LENDER, PERIODO, all_of(instrument_var)) |>
    drop_na(all_of(instrument_var))
  
  unique_periods <- sort(unique(df.W$PERIODO))
  all_nodes      <- tokeep     # all_nodes <- sort(unique(c(df.W$INS_COD, df.W$LENDER)))
  
  create_adjacency_matrix <- function(tt) {
    temp  <- df.W|> filter(PERIODO == tt)
    edges <- temp|> transmute(from  = INS_COD,                   # W: B->L
                              to    = LENDER,
                              weight= .data[[instrument_var]])
    g     <- graph_from_data_frame(d        = edges, 
                                   directed = TRUE, 
                                   vertices = all_nodes)
    A    <- as_adjacency_matrix(g, 
                                attr = "weight", 
                                type = "both", 
                                names= TRUE)
    W <- matrix(0, 
                nrow     = length(all_nodes), 
                ncol     = length(all_nodes),
                dimnames = list(all_nodes, all_nodes))
    W[rownames(A), colnames(A)] <- as.matrix(A)
    
    W_df        <- as.data.frame(W)
    W_df$PERIOD <- tt
    W_df$ID     <- rownames(W_df)
    W_df        <- W_df|>select(PERIOD, ID, all_of(all_nodes))
    return(W_df)
  }
  
  all_W_df <- do.call(rbind, lapply(unique_periods, create_adjacency_matrix))
  return(all_W_df)
}


#####
# II.1. C18: First version for the adjacency matrix W_t
#            Variation across terms to maturity

all_W_df_Tot_1 <- fn_W_df(C18_Muestra, 
                          maturity_value = 1, 
                          instrument_var = "TOT_OBLIG")
all_W_df_Tot_2 <- fn_W_df(C18_Muestra, 
                          maturity_value = 2, 
                          instrument_var = "TOT_OBLIG")
all_W_df_Tot_3 <- fn_W_df(C18_Muestra, 
                          maturity_value = 3, 
                          instrument_var = "TOT_OBLIG")

all_W_df_Tot_sum  <- rbind(all_W_df_Tot_1, all_W_df_Tot_2, all_W_df_Tot_3)
network_node_cols <- names(all_W_df_Tot_sum)[!(names(all_W_df_Tot_sum) %in% c("PERIOD", "ID"))]
all_W_df_Tot      <- all_W_df_Tot_sum%>%
  group_by(PERIOD, ID) %>%
  summarise(
    across(all_of(network_node_cols), sum, .names = "{.col}"),
    .groups = "drop"
  )%>%
  ungroup() 

save(all_W_df_Tot_1, file = paste0(wd,"/Data/Processed/all_W_df_Tot_1.Rda"))
save(all_W_df_Tot_2, file = paste0(wd,"/Data/Processed/all_W_df_Tot_2.Rda"))
save(all_W_df_Tot_3, file = paste0(wd,"/Data/Processed/all_W_df_Tot_3.Rda"))
save(all_W_df_Tot,   file = paste0(wd,"/Data/Processed/all_W_df_Tot.Rda"))

#####
# II.2. C18: second version for the adjacency matrix
#            Variation across maturity for (only) financial instruments

all_W_df_FIN_1 <- fn_W_df(C18_Muestra, 
                          maturity_value = 1, 
                          instrument_var = "C18_CONTR_DER_FIN")

all_W_df_FIN_2 <- fn_W_df(C18_Muestra, 
                          maturity_value = 2, 
                          instrument_var = "C18_CONTR_DER_FIN")

all_W_df_FIN_3 <- fn_W_df(C18_Muestra, 
                          maturity_value = 3, 
                          instrument_var = "C18_CONTR_DER_FIN")

all_W_df_FIN_sum  <- rbind(all_W_df_FIN_1, all_W_df_FIN_2, all_W_df_FIN_3)
network_node_cols <- names(all_W_df_FIN_sum)[!(names(all_W_df_FIN_sum) %in% c("PERIOD", "ID"))]
all_W_df_FIN      <- all_W_df_FIN_sum%>%
  group_by(PERIOD, ID) %>%
  summarise(
    across(all_of(network_node_cols), sum, .names = "{.col}"),
    .groups = "drop"
  )%>%
  ungroup() 

save(all_W_df_FIN_1, file = paste0(wd,"/Data/Processed/all_W_df_FIN_1.Rda"))
save(all_W_df_FIN_2, file = paste0(wd,"/Data/Processed/all_W_df_FIN_2.Rda"))
save(all_W_df_FIN_3, file = paste0(wd,"/Data/Processed/all_W_df_FIN_3.Rda"))
save(all_W_df_FIN,   file = paste0(wd,"/Data/Processed/all_W_df_FIN.Rda"))


#####
# II.3. C18: Third version for the adjacency matrix
#            Variation across maturity for (only) instruments related to traditional funding, such as term deposits

# Alternative settings using "Traditional Funding Network" 
# and C18_PZO_RES_VENC == 1,2,3,all:
all_W_df_CAP_1 <- fn_W_df(C18_Muestra, 
                          maturity_value = 1, 
                          instrument_var = "C18_DEP_CAP")
all_W_df_CAP_2 <- fn_W_df(C18_Muestra, 
                          maturity_value = 2, 
                          instrument_var = "C18_DEP_CAP")
all_W_df_CAP_3 <- fn_W_df(C18_Muestra, 
                          maturity_value = 3, 
                          instrument_var = "C18_DEP_CAP")
all_W_df_CAP_sum  <- rbind(all_W_df_CAP_1, all_W_df_CAP_2, all_W_df_CAP_3)
network_node_cols <- names(all_W_df_CAP_sum)[!(names(all_W_df_CAP_sum) %in% c("PERIOD", "ID"))]
all_W_df_CAP   <- all_W_df_CAP_sum%>%
  group_by(PERIOD, ID) %>%
  summarise(
    across(all_of(network_node_cols), sum, .names = "{.col}"),
    .groups = "drop"
  )%>%
  ungroup()

save(all_W_df_CAP, file = paste0(wd,"/Data/Processed/all_W_df_CAP.Rda"))


#####
# II.4. C18: Fourth version for the adjacency matrix
#            Variation across maturity for (only) total obligations subtracting the value of collateral

# Alternative settings using "Unsecured Exposure Network" 
# and C18_PZO_RES_VENC == 1,2,3,all:
all_W_df_UNSEC_1 <- fn_W_df(C18_Muestra, 
                            maturity_value = 1, 
                            instrument_var = "UNSEC_OBL")
all_W_df_UNSEC_2 <- fn_W_df(C18_Muestra, 
                            maturity_value = 2, 
                            instrument_var = "UNSEC_OBL")
all_W_df_UNSEC_3 <- fn_W_df(C18_Muestra, 
                            maturity_value = 3, 
                            instrument_var = "UNSEC_OBL")
all_W_df_UNSEC_sum<- rbind(all_W_df_UNSEC_1, all_W_df_UNSEC_2, all_W_df_UNSEC_3)
network_node_cols <- names(all_W_df_UNSEC_sum)[!(names(all_W_df_UNSEC_sum) %in% c("PERIOD", "ID"))]
all_W_df_UNSEC    <- all_W_df_UNSEC_sum%>%
  group_by(PERIOD, ID) %>%
  summarise(
    across(all_of(network_node_cols), sum, .names = "{.col}"),
    .groups = "drop"
  )%>%
  ungroup()

save(all_W_df_UNSEC, file = paste0(wd,"/Data/Processed/all_W_df_UNSEC.Rda"))

####################################
# End of dataprocessing.R
# LChanci
####################################