# --- 0. 載入必要套件 ---
# 請確保已安裝以下套件: install.packages(c("shiny", "shinyjs", "metafor", "ggplot2", "dplyr", "DT"))
library(shiny)
library(shinyjs) # For enable/disable
library(ggplot2)
library(stats)   # 使用 pnorm, qnorm, qt 等基本統計函數
library(dplyr)   # 用於處理結果表格的資料
library(DT)      # 用於顯示互動式表格

set.seed(123)    # 設定隨機種子以確保預測檢查中的模擬結果可重現

# --- 定義 *預設* 先驗分配 (固定參數，放在 UI/Server 之外) ---
# Note: These priors are defined on the log(effect size) scale (e.g., log(OR/RR/HR) or lnRR).
# Users can now modify these in the UI.
default_priors <- list(
  skeptic = list(name = "懷疑型 (Skeptic)", mean = 0, sd = 0.205),       # 95% CI on RR scale approx 0.67-1.5
  neutral = list(name = "中性型 (Neutral)", mean = 0, sd = 0.355),       # 95% CI on RR scale approx 0.5-2.0
  optimistic = list(name = "樂觀型 (Optimistic)", mean = log(0.8), sd = 0.215), # Assumes benefit (RR=0.8). 95% CI on RR scale approx 0.54-1.19
  pessimistic = list(name = "悲觀型 (Pessimistic)", mean = log(1.2), sd = 0.18)  # Assumes harm (RR=1.2). 95% CI on RR scale approx 0.84-1.71. P(RR<1) = pnorm(0, log(1.2), 0.18) ≈ 0.156
)

# --- 輔助函數 (放在 UI/Server 之外) ---

# 貝氏更新函數 (常態-常態模型)
bayesian_update_normal <- function(theta_hat, se_theta, prior_mean, prior_sd) {
  # Ensure prior_sd is positive and finite, se_theta is positive and finite
  if (!is.finite(prior_sd) || prior_sd <= 0) {
      stop("Prior standard deviation must be positive and finite.")
  }
  if (!is.finite(se_theta) || se_theta <= 0) {
      stop("Likelihood standard error (SE) must be positive and finite.")
  }
  if (!is.finite(theta_hat) || !is.finite(prior_mean)) {
      stop("Likelihood estimate (theta_hat) and prior mean must be finite.")
  }

  like_precision <- 1 / se_theta^2
  prior_precision <- 1 / prior_sd^2
  post_precision <- like_precision + prior_precision
  # Check for potential division by zero if precisions are extremely small
  if (post_precision <= 1e-12) { # Use a small threshold
      stop("Posterior precision is too close to zero, likely due to very large variances/SE.")
  }
  post_variance <- 1 / post_precision
  post_sd <- sqrt(post_variance)
  post_mean <- (like_precision * theta_hat + prior_precision * prior_mean) / post_precision
  return(list(mean = post_mean, sd = post_sd))
}

# 計算結果函數 (HPDI, 機率) - Adapted to handle different scales
calculate_results <- function(post_mean, post_sd, poi_orig, is_ratio_or_lnrr) {
  # poi_orig is the POI on the original scale (e.g., 0.9 for RR/OR/HR, or a specific MD value)
  # is_ratio_or_lnrr is TRUE if the analysis scale is log(OR/RR/HR) or lnRR

  if (is_ratio_or_lnrr) {
    # Analysis is on log scale
    null_value <- 0 # Null is log(1) = 0
    if (poi_orig <= 0) {
        warning("POI for ratio measures should be > 0. Calculation P(Effect < POI) might be invalid.")
        log_poi <- -Inf
    } else {
        log_poi <- log(poi_orig)
    }
    calc_poi <- log_poi
  } else {
    stop("MD calculation path should no longer be reached.") # Obsolete path
  }

  # Calculate results on the analysis scale (log)
  alpha <- 0.05
  hpdi_lower <- qnorm(alpha / 2, mean = post_mean, sd = post_sd)
  hpdi_upper <- qnorm(1 - alpha / 2, mean = post_mean, sd = post_sd)
  hpdi_calc_scale <- c(hpdi_lower, hpdi_upper)

  prob_lt_null <- pnorm(null_value, mean = post_mean, sd = post_sd) # P(log(Effect) < 0) <=> P(Effect < 1)
  prob_lt_poi <- pnorm(calc_poi, mean = post_mean, sd = post_sd)    # P(log(Effect) < log(POI)) <=> P(Effect < POI)

  # Calculate posterior median and 95% CI for the *original* scale
  if (is_ratio_or_lnrr) {
      post_median_orig <- exp(post_mean) # Median on log scale transforms to median on original scale
      post_ci_orig <- exp(hpdi_calc_scale) # Transform HPDI bounds
  } else { # MD - THIS BLOCK IS NOW OBSOLETE
      stop("MD calculation path should no longer be reached.")
  }

  return(list(hpdi_calc_scale = hpdi_calc_scale, # HPDI on calculation scale (log)
              prob_lt_null = prob_lt_null,
              prob_lt_poi = prob_lt_poi,
              poi_orig_used = poi_orig, # Store the original scale POI
              post_median_orig = post_median_orig, # Original scale median
              post_ci_orig = post_ci_orig)) # Original scale 95% CI/HPDI
}

# 先驗預測檢查函數
prior_predictive_check <- function(theta_hat, se_theta, prior_mean, prior_sd, n_sim = 10000) {
  if (is.na(prior_mean) || is.na(prior_sd) || is.na(se_theta)) return(NA) # Handle NA inputs
  prior_theta_sim <- rnorm(n_sim, mean = prior_mean, sd = prior_sd)
  prior_pred_sim <- rnorm(n_sim, mean = prior_theta_sim, sd = se_theta) # Simulate data given prior thetas
  return(prior_pred_sim)
}

# 後驗預測檢查函數
posterior_predictive_check <- function(theta_hat, se_theta, post_mean, post_sd, n_sim = 10000) {
  if (is.na(post_mean) || is.na(post_sd) || is.na(se_theta)) return(NA) # Handle NA inputs
  post_theta_sim <- rnorm(n_sim, mean = post_mean, sd = post_sd) # Sample thetas from posterior
  post_pred_rep <- rnorm(n_sim, mean = post_theta_sim, sd = se_theta) # Simulate replicated data given posterior thetas
  return(post_pred_rep)
}

# Function to calculate and format prior 95% CI
format_prior_ci <- function(mean_log, sd_log) {
  if (!is.numeric(mean_log) || !is.numeric(sd_log) || sd_log <= 0) {
    return("無效輸入")
  }
  ci_log <- qnorm(c(0.025, 0.975), mean = mean_log, sd = sd_log)
  ci_orig <- exp(ci_log)
  sprintf("%.2f 到 %.2f", ci_orig[1], ci_orig[2])
}

# --- Shiny UI 定義 (使用者介面) ---
ui <- fluidPage(
  useShinyjs(), # Initialize shinyjs
  titlePanel("貝氏再分析工具 (Bayesian Re-analysis Tool)"),

  sidebarLayout(
    sidebarPanel(
      width = 4, # Increase sidebar width slightly
      # --- NEW: Select Input Method ---
      radioButtons("input_method", "步驟 1: 選擇數據輸入方式",
                   choices = c("單一效果量 (OR, RR, HR)" = "single_effect",
                               "兩組連續數據 (計算 lnRR)" = "two_groups"),
                   selected = "single_effect"),
      hr(),

      # --- Conditional Panel for Single Effect Input ---
      conditionalPanel(
        condition = "input.input_method == 'single_effect'",
        h5("步驟 1a: 輸入單一效果量數據 (原始尺度)"),
        helpText("請提供效果量點估計值及其 95% 信賴區間。效果量應 > 0。"),
        selectInput("effect_type", "效果量類型:",
                    choices = c("勝算比 (Odds Ratio)" = "OR",
                                "風險比 (Risk Ratio)" = "RR",
                                "風險函數比 (Hazard Ratio)" = "HR"),
                    selected = "OR"),
        numericInput("effect_mean_orig", "觀察效果估計值:", value = 0.75, step = 0.01, min=0.001),
        numericInput("ci_lower_orig", "95% CI 下界:", value = 0.55, step = 0.01, min=0.001),
        numericInput("ci_upper_orig", "95% CI 上界:", value = 1.02, step = 0.01, min=0.001)
      ),

      # --- Conditional Panel for Two Group Input ---
      conditionalPanel(
        condition = "input.input_method == 'two_groups'",
        h5("步驟 1b: 輸入兩組連續數據 (原始尺度)"),
        helpText("請提供兩組的平均值 (Mean)、95% 信賴區間 (CI) 及樣本數 (N)。平均值應 > 0 以計算對數反應比 (lnRR)。"),
        h6("第一組 (例如：治療組)"),
        fluidRow(
          column(6, numericInput("mean1", "平均值 1:", value = 130, step = 1, min=0.001)),
          column(6, numericInput("n1", "樣本數 1:", value = 50, min = 2, step = 1))
        ),
        fluidRow(
          column(6, numericInput("lower_ci1", "95% CI 下界 1:", value = 125, step = 1)),
          column(6, numericInput("upper_ci1", "95% CI 上界 1:", value = 135, step = 1))
        ),
        hr(),
        h6("第二組 (例如：對照組)"),
        fluidRow(
          column(6, numericInput("mean2", "平均值 2:", value = 140, step = 1, min=0.001)),
          column(6, numericInput("n2", "樣本數 2:", value = 50, min = 2, step = 1))
        ),
        fluidRow(
          column(6, numericInput("lower_ci2", "95% CI 下界 2:", value = 133, step = 1)),
          column(6, numericInput("upper_ci2", "95% CI 上界 2:", value = 147, step = 1))
        )
      ),
      hr(),
      h4("步驟 2: 設定分析與先驗參數"),
      numericInput("poi", "感興趣點 (POI) (原始尺度):", value = 0.9, step = 0.05, min=0.001),
      helpText("輸入效果閾值（對 OR/RR/HR/lnRR，此值 > 0）。將計算 P(真實效果 < POI)。"),
      hr(),

      # --- NEW: User Modifiable Prior Settings ---
      h5("先驗分配參數 (對數尺度)"),
      helpText("修改下方先驗分配的平均值 (Mean) 與標準差 (SD)，或使用預設值。"),
      checkboxInput("use_defaults", "使用預設先驗參數", value = TRUE),
      br(),

      # Skeptic Prior Inputs
      strong("懷疑型 (Skeptic):"),
      fluidRow(
        column(6, numericInput("prior_mean_skeptic", "Mean (log)", value = default_priors$skeptic$mean, step = 0.01)),
        column(6, numericInput("prior_sd_skeptic", "SD (log)", value = default_priors$skeptic$sd, min = 0.001, step = 0.01))
      ),
      textOutput("prior_ci_skeptic"),
      br(),

      # Neutral Prior Inputs
      strong("中性型 (Neutral):"),
      fluidRow(
        column(6, numericInput("prior_mean_neutral", "Mean (log)", value = default_priors$neutral$mean, step = 0.01)),
        column(6, numericInput("prior_sd_neutral", "SD (log)", value = default_priors$neutral$sd, min = 0.001, step = 0.01))
      ),
      textOutput("prior_ci_neutral"),
      br(),

      # Optimistic Prior Inputs
      strong("樂觀型 (Optimistic):"),
      fluidRow(
        column(6, numericInput("prior_mean_optimistic", "Mean (log)", value = default_priors$optimistic$mean, step = 0.01)),
        column(6, numericInput("prior_sd_optimistic", "SD (log)", value = default_priors$optimistic$sd, min = 0.001, step = 0.01))
      ),
      textOutput("prior_ci_optimistic"),
      br(),

      # Pessimistic Prior Inputs
      strong("悲觀型 (Pessimistic):"),
      fluidRow(
        column(6, numericInput("prior_mean_pessimistic", "Mean (log)", value = default_priors$pessimistic$mean, step = 0.01)),
        column(6, numericInput("prior_sd_pessimistic", "SD (log)", value = default_priors$pessimistic$sd, min = 0.001, step = 0.01))
      ),
      textOutput("prior_ci_pessimistic"),
      hr(),

      # --- End User Modifiable Priors ---

      h4("步驟 3: 選擇先驗進行預測檢查"),
      selectInput("plot_prior_select", "選擇用於預測檢查圖的先驗:",
                  choices = c("懷疑型" = "skeptic",
                              "中性型" = "neutral",
                              "樂觀型" = "optimistic",
                              "悲觀型" = "pessimistic"),
                  selected = "neutral"),
      helpText("先驗分配代表分析前對效果量的信念（在對數尺度上定義）。詳細說明請見「背景 & 先驗」分頁。"),
      hr(),
      tags$h5("參考文獻"),
      tags$a(href = "https://www.atsjournals.org/doi/suppl/10.1164/rccm.202006-2381CP?role=tab",
             target = "_blank",
             "Using Bayesian Methods to Augment the Interpretation of Critical Care Trials. An Overview of Theory and Example Reanalysis of the ART Trial"
      ),
      tags$br(),
      tags$br(),
      tags$a(href = "https://onlinelibrary.wiley.com/doi/full/10.1002/sim.10114",
             target = "_blank",
             "Standardization and other approaches to meta-analyze differences in means"
      ),
      hr()
    ),

    mainPanel(
      width = 8,
      tabsetPanel(
        tabPanel("分析結果 (Analysis Results)",
          h4("分析結果"),
          p("以下顯示根據您輸入的數據計算出的效果量與標準誤，以及結合數據與不同先驗分配後得到的貝氏分析結果。"),
          h5("計算得到的效果量與標準誤 (SE)"),
          helpText("效果量及其 SE 在對數尺度上計算（log(OR/RR/HR) 或 lnRR）。"),
          verbatimTextOutput("effect_se_output"),
          hr(),
          h5("貝氏再分析結果摘要表"),
          DTOutput("results_table"),
          helpText(HTML("<b>表格欄位說明:</b><ul>
                          <li><b>先驗 (Prior):</b> 使用的先驗分配名稱。</li>
                          <li><b>先驗平均/標準差 (對數尺度):</b> 您設定（或預設）的先驗參數。</li>
                          <li><b>後驗平均/標準差 (計算尺度):</b> 效果量估計及其不確定性（對數尺度）。</li>
                          <li><b>95% 置信區間 HPDI (計算尺度) 下界/上界:</b> 後驗區間（對數尺度）。</li>
                          <li><b>後驗中位數 (原始):</b> 後驗分配轉換回原始尺度的中位數估計。</li>
                          <li><b>95% CI (原始) 下界/上界:</b> 後驗分配轉換回原始尺度的 95% 可信區間。</li>
                          <li><b>P(效果 < Null):</b> 效果小於無效假設（RR/OR/HR/RR < 1）的後驗機率。</li>
                          <li><b>P(效果 < POI):</b> 真實效果量小於您設定的感興趣點 (POI) 的後驗機率。</li></ul>
                          <b>結果判讀:</b> 比較不同先驗下的後驗結果。觀察後驗估計是否偏移、區間是否縮小、以及相關機率是否支持您的研究假設。")),
          hr(),
          h5("預測檢查 (Predictive Checks)"),
          helpText("預測檢查用於評估模型是否能合理地產生出我們實際觀察到的數據（在計算尺度上）。請在左側選擇要繪圖的先驗。"),
          plotOutput("predictive_plot"),
          verbatimTextOutput("plot_interpretation")
        ),
        tabPanel("背景 & 先驗 (Background & Priors)",
          h4("貝氏定理 (Bayes' Theorem)"),
          p("貝氏定理提供了一種基於新證據更新信念（先驗機率）的方法，以產生更新的信念（後驗機率）。其核心概念是："),
          withMathJax("$$ P(\\theta | Data) = \\frac{P(Data | \\theta) \\times P(\\theta)}{P(Data)} $$"),
          p("其中："),
          tags$ul(
            tags$li(withMathJax("\\( P(\\theta | Data) \\) 是後驗機率：給定觀察到的數據 (Data) 後，參數 \\( \\theta \\) （例如，真實的效果量）的機率。")),
            tags$li(withMathJax("\\( P(Data | \\theta) \\) 是似然率 (Likelihood)：假設參數 \\( \\theta \\) 為某特定值時，觀察到目前數據的可能性。")),
            tags$li(withMathJax("\\( P(\\theta) \\) 是先驗機率 (Prior Probability)：在觀察數據之前，我們對參數 \\( \\theta \\) 的信念。")),
            tags$li(withMathJax("\\( P(Data) \\) 是證據 (Evidence)或邊際似然率：無論參數是多少，觀察到數據的總體機率，作為歸一化常數。"))
          ),
          p(withMathJax("簡單來說， 後驗 \\( \\propto \\) 似然率 \\( \\times \\) 先驗。")),
          h4("高斯先驗與似然率的解析解"),
          p(withMathJax("在此應用程式中，我們假設效果量（在對數尺度上）的似然率 \\( P(Data | \\theta) \\) 遵循常態分佈，其平均值為觀察到的效果量估計值 (\\( \\hat{\\theta} \\))，標準差為該估計值的標準誤 (SE)。我們同時使用常態分佈作為先驗 \\( P(\\theta) \\) 。")),
          tags$ul(
            tags$li(strong("共軛先驗 (Conjugate Priors):"), "當先驗分佈和似然函數的形式使得計算出的後驗分佈屬於同一分佈家族時，該先驗稱為共軛先驗。常態分佈是（當變異數已知或估計時）常態似然率平均參數的共軛先驗."),
            withMathJax(
              tags$ul(
                tags$li(strong("常態-常態模型與解析解:"), "由於使用了共軛的常態先驗和常態似然率，後驗分佈 \\( P(\\theta | Data) \\) 也會是常態分佈。其參數（後驗平均值 \\( \\mu_{post} \\) 和後驗標準差 \\( \\sigma_{post} \\)）可以直接透過以下公式計算，無需複雜模擬："),
                tags$li("$$ \\text{似然率精密度: } \\tau_{like} = \\frac{1}{SE^2} $$"),
                tags$li("$$ \\text{先驗精密度: } \\tau_{prior} = \\frac{1}{\\sigma_{prior}^2} $$"),
                tags$li("$$ \\text{後驗精密度: } \\tau_{post} = \\tau_{like} + \\tau_{prior} $$"),
                tags$li("$$ \\text{後驗平均值: } \\mu_{post} = \\frac{\\tau_{like} \\times \\hat{\\theta} + \\tau_{prior} \\times \\mu_{prior}}{\\tau_{post}} $$"),
                tags$li("$$ \\text{後驗標準差: } \\sigma_{post} = \\sqrt{\\frac{1}{\\tau_{post}}} $$"),
                tags$li(strong("非共軛先驗:"), "如果選擇的先驗分佈不是似然率的共軛先驗（例如，使用 t 分佈先驗與常態似然率），則通常無法得到簡單的後驗分佈公式。在這種情況下，需要使用數值方法，如馬可夫鏈蒙地卡羅 (Markov Chain Monte Carlo, MCMC) 模擬，來從後驗分佈中抽樣並估計其特性。")
              )
            )
          ),
          hr(),
          h4("先驗分配詳細說明 (Prior Explanations)"),
          p("這些先驗分配是在效果量的 **對數尺度** 上定義的（例如 log(OR)，log(RR)，log(HR)，lnRR）。您可以修改側邊欄中的參數。下方說明基於預設值。"),
          p(strong("大部分的 RR, OR, HR, 反應比介於 0.5-2："), "「悲觀」不適用於臨床試驗，因為只有預期「樂觀」的效果你才會想要進行試驗。「樂觀」和「悲觀」的解釋取決於效果量指標的含義。在此應用中，預設假設效果量 < 1（對數尺度 < 0）代表有益效果（例如，降低風險）。如果您的研究中，效果量 > 1 代表有益效果（例如，增加成功率的反應比），那麼原本的「樂觀」先驗（均值 < 0）實際上變成了「悲觀」，而「悲觀」先驗（均值 > 0）則變成了「樂觀」。請根據您的具體情況調整參數和解釋。當先驗是常態分布、只有一個參數、無資訊/平坦先驗（均勻分布、gaussian(0, 1000)）、似然率的變異數已知、統計模型正確時，單尾 p 值幾乎等於後驗機率。當統計測試（t 試驗、線性回歸的 Wald 試驗）是對稱分布時，單尾 p 值等於雙尾 p 值/2。非對稱分布的統計測試：卡方分布、F 分布（變異數分析）、邏輯回歸的 Wald 試驗、無母數分析、似然比試驗（模型比較）、log-rank（生存分析）等。"),
          br(),
          h5("1. 懷疑型 (Skeptic Prior): 預設 N(mean = 0, sd = 0.205)"),
          uiOutput("desc_skeptic"),
          h5("2. 中性型 (Neutral Prior): 預設 N(mean = 0, sd = 0.355)"),
          uiOutput("desc_neutral"),
          h5("3. 樂觀型 (Optimistic Prior): 預設 N(mean = log(0.8) ≈ -0.223, sd = 0.215)"),
          uiOutput("desc_optimistic"),
          p(strong("選擇建議："), "樂觀先驗的平均效果量可以基於先前的研究、臨床意義的閾值，或甚至基於原始隨機對照試驗 (RCT) 計算樣本數時所假設的效果大小來設定."),
          h5("4. 悲觀型 (Pessimistic Prior): 預設 N(mean = log(1.2) ≈ 0.182, sd = 0.18)"),
          uiOutput("desc_pessimistic")
        )
      )
    )
  )
)

# --- Shiny Server 定義 (伺服器邏輯) ---
server <- function(input, output, session) {
  # --- Reactive values for priors ---
  current_priors <- reactiveValues(
    skeptic = list(name = "懷疑型 (Skeptic)", mean = default_priors$skeptic$mean, sd = default_priors$skeptic$sd),
    neutral = list(name = "中性型 (Neutral)", mean = default_priors$neutral$mean, sd = default_priors$neutral$sd),
    optimistic = list(name = "樂觀型 (Optimistic)", mean = default_priors$optimistic$mean, sd = default_priors$optimistic$sd),
    pessimistic = list(name = "悲觀型 (Pessimistic)", mean = default_priors$pessimistic$mean, sd = default_priors$pessimistic$sd)
  )

  # --- Observer to update priors based on UI inputs ---
  observe({
    if (input$use_defaults) {
      # Reset to defaults and disable inputs
      current_priors$skeptic$mean <- default_priors$skeptic$mean
      current_priors$skeptic$sd <- default_priors$skeptic$sd
      current_priors$neutral$mean <- default_priors$neutral$mean
      current_priors$neutral$sd <- default_priors$neutral$sd
      current_priors$optimistic$mean <- default_priors$optimistic$mean
      current_priors$optimistic$sd <- default_priors$optimistic$sd
      current_priors$pessimistic$mean <- default_priors$pessimistic$mean
      current_priors$pessimistic$sd <- default_priors$pessimistic$sd

      updateNumericInput(session, "prior_mean_skeptic", value = default_priors$skeptic$mean)
      updateNumericInput(session, "prior_sd_skeptic", value = default_priors$skeptic$sd)
      updateNumericInput(session, "prior_mean_neutral", value = default_priors$neutral$mean)
      updateNumericInput(session, "prior_sd_neutral", value = default_priors$neutral$sd)
      updateNumericInput(session, "prior_mean_optimistic", value = default_priors$optimistic$mean)
      updateNumericInput(session, "prior_sd_optimistic", value = default_priors$optimistic$sd)
      updateNumericInput(session, "prior_mean_pessimistic", value = default_priors$pessimistic$mean)
      updateNumericInput(session, "prior_sd_pessimistic", value = default_priors$pessimistic$sd)

      disable("prior_mean_skeptic"); disable("prior_sd_skeptic")
      disable("prior_mean_neutral"); disable("prior_sd_neutral")
      disable("prior_mean_optimistic"); disable("prior_sd_optimistic")
      disable("prior_mean_pessimistic"); disable("prior_sd_pessimistic")
    } else {
      # Enable inputs and use their current values
      enable("prior_mean_skeptic"); enable("prior_sd_skeptic")
      enable("prior_mean_neutral"); enable("prior_sd_neutral")
      enable("prior_mean_optimistic"); enable("prior_sd_optimistic")
      enable("prior_mean_pessimistic"); enable("prior_sd_pessimistic")

      # Validate and update reactive values from UI inputs
      req(input$prior_mean_skeptic, input$prior_sd_skeptic,
          input$prior_mean_neutral, input$prior_sd_neutral,
          input$prior_mean_optimistic, input$prior_sd_optimistic,
          input$prior_mean_pessimistic, input$prior_sd_pessimistic)

      # Add validation for positive SDs
      validate(
        need(input$prior_sd_skeptic > 0, "懷疑型 SD 必須 > 0"),
        need(input$prior_sd_neutral > 0, "中性型 SD 必須 > 0"),
        need(input$prior_sd_optimistic > 0, "樂觀型 SD 必須 > 0"),
        need(input$prior_sd_pessimistic > 0, "悲觀型 SD 必須 > 0")
      )

      current_priors$skeptic$mean <- input$prior_mean_skeptic
      current_priors$skeptic$sd <- input$prior_sd_skeptic
      current_priors$neutral$mean <- input$prior_mean_neutral
      current_priors$neutral$sd <- input$prior_sd_neutral
      current_priors$optimistic$mean <- input$prior_mean_optimistic
      current_priors$optimistic$sd <- input$prior_sd_optimistic
      current_priors$pessimistic$mean <- input$prior_mean_pessimistic
      current_priors$pessimistic$sd <- input$prior_sd_pessimistic
    }
  })

  # --- Render dynamic prior CIs in sidebar ---
  output$prior_ci_skeptic <- renderText({ paste("=> 原始尺度 95% CI 約: 1 (", format_prior_ci(current_priors$skeptic$mean, current_priors$skeptic$sd)) })
  output$prior_ci_neutral <- renderText({ paste("=> 原始尺度 95% CI 約: 1 (", format_prior_ci(current_priors$neutral$mean, current_priors$neutral$sd)) })
  output$prior_ci_optimistic <- renderText({ paste("=> 原始尺度 95% CI 約: 0.8 (", format_prior_ci(current_priors$optimistic$mean, current_priors$optimistic$sd)) })
  output$prior_ci_pessimistic <- renderText({ paste("=> 原始尺度 95% CI 約: 1.2 (", format_prior_ci(current_priors$pessimistic$mean, current_priors$pessimistic$sd)) })

  # --- Render dynamic prior descriptions in Background Tab ---
  output$desc_skeptic <- renderUI({
    m <- current_priors$skeptic$mean
    s <- current_priors$skeptic$sd
    ci_str <- format_prior_ci(m, s)
    tags$ul(
      tags$li(strong("基本原理:"), sprintf("此先驗將信念集中在平均值 = %.3f（對數尺度），對應於原始尺度效果 = %.2f。", m, exp(m))),
      tags$li(strong("確定性:"), sprintf("標準差 SD = %.3f 代表信念的確定性程度。數值越小，信念越集中。", s)),
      tags$li(strong("解釋:"), sprintf("此設定意味著先驗假設真實效果（原始尺度）有 95%% 的機率位於 %s。預設值 (mean=0, sd=0.205) 代表相對強烈的「無效果」信念。", ci_str))
    )
  })
  output$desc_neutral <- renderUI({
    m <- current_priors$neutral$mean
    s <- current_priors$neutral$sd
    ci_str <- format_prior_ci(m, s)
    tags$ul(
      tags$li(strong("基本原理:"), sprintf("此先驗將信念集中在平均值 = %.3f（對數尺度），對應於原始尺度效果 = %.2f。", m, exp(m))),
      tags$li(strong("確定性:"), sprintf("標準差 SD = %.3f。預設值 (sd=0.355) 比懷疑型先驗更大，反映了較低的確定性或較弱的資訊。", s)),
      tags$li(strong("解釋:"), sprintf("此設定意味著先驗假設真實效果（原始尺度）有 95%% 的機率位於 %s。預設值 (mean=0, sd=0.355) 通常被認為是「弱資訊性」先驗。", ci_str))
    )
  })
  output$desc_optimistic <- renderUI({
    m <- current_priors$optimistic$mean
    s <- current_priors$optimistic$sd
    ci_str <- format_prior_ci(m, s)
    prob_lt_1 <- pnorm(0, mean=m, sd=s) # P(log(Effect) < 0) <=> P(Effect < 1)
    tags$ul(
      tags$li(strong("基本原理:"), sprintf("此先驗將信念集中在平均值 = %.3f（對數尺度），對應於原始尺度效果 = %.2f。", m, exp(m))),
      tags$li(strong("確定性:"), sprintf("標準差 SD = %.3f。", s)),
      tags$li(strong("解釋:"), sprintf("此設定代表了一種先驗信念，即干預可能是有效的（預設 mean=log(0.8) 假設效果降低20%%）。95%% 的先驗機率範圍大約是 %s。對於此先驗，效果 < 1 的先驗機率是 %.2f。", ci_str, prob_lt_1))
    )
  })
  output$desc_pessimistic <- renderUI({
    m <- current_priors$pessimistic$mean
    s <- current_priors$pessimistic$sd
    ci_str <- format_prior_ci(m, s)
    prob_lt_1 <- pnorm(0, mean=m, sd=s) # P(log(Effect) < 0) <=> P(Effect < 1)
    tags$ul(
      tags$li(strong("基本原理:"), sprintf("此先驗將信念集中在平均值 = %.3f（對數尺度），對應於原始尺度效果 = %.2f。", m, exp(m))),
      tags$li(strong("確定性:"), sprintf("標準差 SD = %.3f。", s)),
      tags$li(strong("解釋:"), sprintf("此設定代表了一種先驗信念，即干預可能是有害的（預設 mean=log(1.2) 假設效果增加20%%）。95%% 的先驗機率範圍大約是 %s。對於此先驗，效果 < 1 的先驗機率是 %.2f (預設值下約為 0.156)。", ci_str, prob_lt_1))
    )
  })

  # --- 反應性計算核心 ---
  analysis_output <- reactive({
    # --- Read common inputs ---
    input_method <- input$input_method
    point_of_interest <- input$poi # POI on original scale
    req(point_of_interest)

    # --- Read current prior settings ---
    priors_to_use <- reactiveValuesToList(current_priors) # Get the current list of priors

    # --- Process inputs based on selected method ---
    if (input_method == "single_effect") {
      # --- Validate and process single effect inputs ---
      req(input$effect_type, input$effect_mean_orig, input$ci_lower_orig, input$ci_upper_orig)
      effect_type <- input$effect_type

      validate(
        need(input$effect_mean_orig > 0, "觀察效果估計值 (原始尺度) 必須 > 0。"),
        need(input$ci_lower_orig > 0, "信賴區間下界 (原始尺度) 必須 > 0。"),
        need(input$ci_upper_orig > 0, "信賴區間上界 (原始尺度) 必須 > 0。"),
        need(input$ci_upper_orig > input$ci_lower_orig, "信賴區間上界 (原始尺度) 必須大於下界。"),
        need(point_of_interest > 0, "感興趣點 (POI) 必須 > 0。")
      )

      # Log transform inputs
      theta_hat <- log(input$effect_mean_orig)
      ci_low <- log(input$ci_lower_orig)
      ci_upp <- log(input$ci_upper_orig)
      # Calculate SE from log-transformed CI
      se_theta <- (ci_upp - ci_low) / (2 * qnorm(0.975))
      scale_label <- "對數尺度"
      calc_scale_is_log <- TRUE
      null_value <- 0 # Null on log scale is 0
      null_label <- "1" # Null on original ratio scale is 1
      # Map effect type to full label
      effect_type_map <- c("OR" = "勝算比 (Odds Ratio)", "RR" = "風險比 (Risk Ratio)", "HR" = "風險函數比 (Hazard Ratio)")
      effect_type_label_full <- effect_type_map[effect_type]

    } else if (input_method == "two_groups") {
      # --- Validate and process two-group inputs ---
      req(input$mean1, input$lower_ci1, input$upper_ci1, input$n1,
          input$mean2, input$lower_ci2, input$upper_ci2, input$n2)
      validate(
        need(input$n1 >= 2, "第一組樣本數 (N 1) 必須至少為 2。"),
        need(input$n2 >= 2, "第二組樣本數 (N 2) 必須至少為 2。"),
        need(input$mean1 > 0, "第一組平均值必須 > 0 以計算 lnRR。"),
        need(input$mean2 > 0, "第二組平均值必須 > 0 以計算 lnRR。"),
        need(input$upper_ci1 > input$lower_ci1, "第一組信賴區間上界必須大於下界。"),
        need(input$upper_ci2 > input$lower_ci2, "第二組信賴區間上界必須大於下界。"),
        need(point_of_interest > 0, "感興趣點 (POI) 對於反應比必須 > 0。")
      )

      m1 <- input$mean1; lci1 <- input$lower_ci1; uci1 <- input$upper_ci1; n1 <- input$n1
      m2 <- input$mean2; lci2 <- input$lower_ci2; uci2 <- input$upper_ci2; n2 <- input$n2

      # Calculate SD from Mean and 95% CI using t-distribution
      se1 <- (uci1 - lci1) / (2 * qt(0.975, df = n1 - 1)); sd1 <- se1 * sqrt(n1)
      se2 <- (uci2 - lci2) / (2 * qt(0.975, df = n2 - 1)); sd2 <- se2 * sqrt(n2)
      validate(need(sd1 > 1e-9, "從第一組CI計算出的SD無效或為零。"), need(sd2 > 1e-9, "從第二組CI計算出的SD無效或為零。"))

      # Calculate lnRR and its SE
      theta_hat <- log(m1) - log(m2) # lnRR
      var_lnRR <- (sd1^2 / (n1 * m1^2)) + (sd2^2 / (n2 * m2^2))
      se_theta <- sqrt(var_lnRR)      # SE(lnRR)

      scale_label <- "對數反應比 (lnRR)"
      calc_scale_is_log <- TRUE
      null_value <- 0 # Null on log scale is 0
      null_label <- "1" # Null for RR is 1
      effect_type_label_full <- "反應比 (Response Ratio)"

    } else {
      stop("未知的輸入方式")
    }

    # --- Common validation for calculated values ---
    validate(
      need(is.finite(theta_hat) && is.finite(se_theta), "計算出的效果量或 SE 不是有限數值。請檢查輸入。"),
      need(se_theta > 1e-9, "計算出的標準誤 (SE) 過小或無效。")
    )

    # --- 針對每種先驗分配進行貝氏分析 ---
    results_list <- list()
    # Use names from the reactive priors list
    for (prior_key in names(priors_to_use)) {
      prior_spec <- priors_to_use[[prior_key]]
      prior_mean <- prior_spec$mean # From reactive value
      prior_sd <- prior_spec$sd     # From reactive value

      # Input validation for prior SD happens in the observer now,
      # but double-check here just in case.
      if(!is.finite(prior_sd) || prior_sd <= 0) {
        warning(paste("Invalid Prior SD used for", prior_spec$name, "- skipping this prior."))
        results_list[[prior_key]] <- NULL # Skip this prior
        next # Move to next prior
      }

      # 1. Bayesian update (uses theta_hat, se_theta on calculation scale)
      posterior <- bayesian_update_normal(theta_hat, se_theta, prior_mean, prior_sd)
      post_mean <- posterior$mean # Posterior mean on calculation scale (log)
      post_sd <- posterior$sd     # Posterior SD on calculation scale (log)

      # 2. Calculate posterior results (HPDI, probabilities, original scale median/CI)
      analysis_results <- calculate_results(post_mean, post_sd, poi_orig = point_of_interest, is_ratio_or_lnrr = TRUE) # Always TRUE now

      # 3. Predictive checks (uses calculation scale inputs/results)
      # Pass current prior mean/sd to predictive check functions
      prior_pred_samples <- prior_predictive_check(theta_hat, se_theta, prior_mean, prior_sd)
      post_pred_samples <- posterior_predictive_check(theta_hat, se_theta, post_mean, post_sd)

      # 4. Store results
      results_list[[prior_key]] <- list(
        prior_name_nice = prior_spec$name,
        prior_mean = prior_mean, # Log scale (user defined or default)
        prior_sd = prior_sd,     # Log scale (user defined or default)
        posterior_mean_calc = post_mean,  # Calculation scale (log)
        posterior_sd_calc = post_sd,      # Calculation scale (log)
        hpdi_calc = analysis_results$hpdi_calc_scale, # Calculation scale (log)
        prob_lt_null = analysis_results$prob_lt_null, # P(Effect < 1)
        prob_lt_poi = analysis_results$prob_lt_poi, # P(Effect < POI)
        post_median_orig = analysis_results$post_median_orig, # Original scale
        post_ci_orig = analysis_results$post_ci_orig,         # Original scale
        poi_used = analysis_results$poi_orig_used,            # Original scale POI
        prior_pred_samples = prior_pred_samples,              # Calculation scale (log)
        post_pred_samples = post_pred_samples,                # Calculation scale (log)
        observed_theta_hat = theta_hat                        # Calculation scale (log)
      )
    } # End prior loop

    # Filter out any NULL results from skipped priors
    results_list <- results_list[!sapply(results_list, is.null)]

    # --- Return all necessary values ---
    return(list(
      theta_hat = theta_hat,        # Calculation scale effect estimate (log)
      se_theta = se_theta,          # Calculation scale SE (log)
      results = results_list,       # List of results for each valid prior
      priors_used = priors_to_use,  # Include the actual priors used
      poi = point_of_interest,      # Original scale POI
      effect_type_label = effect_type_label_full, # Full label for display
      scale_label = scale_label,    # Label for calculation scale (e.g., "對數尺度", "lnRR")
      null_label = null_label,      # Label for null value on original scale ("1")
      input_method = input_method   # Store input method used
    ))
  }) # End reactive analysis_output

  # --- 輸出渲染區塊 ---

  # Display calculated effect size and SE (in Results tab)
  output$effect_se_output <- renderPrint({
    analysis_data <- analysis_output()
    req(analysis_data)
    theta <- analysis_data$theta_hat
    se <- analysis_data$se_theta
    scale_lab <- analysis_data$scale_label
    effect_lab <- analysis_data$effect_type_label

    effect_desc <- switch(analysis_data$input_method,
                          "single_effect" = effect_lab, # Use OR, RR, HR label
                          "two_groups" = scale_lab, # Use "對數反應比 (lnRR)"
                          "效果量") # Default

    cat(sprintf("計算得到的效果量 (%s): %.4f\n", effect_desc, theta))
    cat(sprintf("標準誤 (SE) (%s): %.4f\n", scale_lab, se))
    cat(sprintf("隱含的概似精密度 (Precision = 1/SE²): %.4f", 1/se^2))
  })

  # Display results summary table (in Results tab)
  output$results_table <- renderDT({
    analysis_data <- analysis_output()
    req(analysis_data)
    res <- analysis_data$results
    # Check if results list is empty (e.g., all priors had invalid SD)
    if (length(res) == 0) {
      return(datatable(data.frame(訊息 = "沒有有效的先驗分配可用於分析。請檢查側邊欄中的先驗參數。"), rownames = FALSE, options = list(searching = FALSE, lengthChange = FALSE, paging = FALSE, info = FALSE)))
    }

    poi_val <- analysis_data$poi
    scale_lab <- analysis_data$scale_label # Calculation scale (log)
    null_lab <- analysis_data$null_label   # Original scale null ("1")

    prior_scale_lab <- "對數尺度" # Priors are always log scale

    df <- bind_rows(lapply(res, function(r) {
      data.frame(
        Prior = r$prior_name_nice,
        Prior_Mean = r$prior_mean,
        Prior_SD = r$prior_sd,
        Posterior_Mean_Calc = r$posterior_mean_calc, # Calc scale (log)
        Posterior_SD_Calc = r$posterior_sd_calc,     # Calc scale (log)
        HPDI_Calc_Lower = r$hpdi_calc[1],            # Calc scale (log)
        HPDI_Calc_Upper = r$hpdi_calc[2],            # Calc scale (log)
        Post_Median_Orig = r$post_median_orig,       # Original scale
        CI_Orig_Lower = r$post_ci_orig[1],           # Original scale
        CI_Orig_Upper = r$post_ci_orig[2],           # Original scale
        Prob_Effect_lt_Null = r$prob_lt_null,
        Prob_Effect_lt_POI = r$prob_lt_poi
      )
    }))

    colnames(df) <- c("先驗",
                     paste("先驗平均 (", prior_scale_lab, ")", sep=""),
                     paste("先驗標準差 (", prior_scale_lab, ")", sep=""),
                     paste("後驗平均 (", scale_lab, ")", sep=""),
                     paste("後驗標準差 (", scale_lab, ")", sep=""),
                     paste("95% HPDI 下界 (", scale_lab, ")", sep=""),
                     paste("95% HPDI 上界 (", scale_lab, ")", sep=""),
                     "後驗中位數 (原始)",
                     "95% CI 下界 (原始)",
                     "95% CI 上界 (原始)",
                     sprintf("P(效果 < %s)", null_lab), # P(Effect < 1)
                     sprintf("P(效果 < %.3f)", poi_val)) # P(Effect < POI)

    # Identify numeric columns first
    numeric_cols <- names(df)[sapply(df, is.numeric)]

    datatable(df,
              rownames = FALSE,
              extensions = 'Buttons',
              options = list(
                scrollX = TRUE,
                searching = FALSE,
                lengthChange = FALSE,
                paging = FALSE,
                info = FALSE,
                dom = 'Bfrtip',
                buttons = list(
                  list(extend = 'copy', text = '複製表格'),
                  list(extend = 'csv', text = '下載 CSV', filename = paste0("bayesian_reanalysis_results_", Sys.Date())),
                  list(extend = 'excel', text = '下載 Excel', filename = paste0("bayesian_reanalysis_results_", Sys.Date()))
                ),
                language = list(url = '//cdn.datatables.net/plug-ins/1.10.25/i18n/Chinese-Taiwan.json')
              ),
              caption = htmltools::tags$caption(style = 'caption-side: top; text-align: left; margin-bottom: 10px;',
                                               '註: 表格中部分估計值在計算尺度（對數尺度），部分在原始尺度，請注意欄位標示。HPDI = 最高後驗密度區間 (Highest Posterior Density Interval).')
    ) %>%
      formatRound(columns = numeric_cols, digits = 3)
  })

  # Display predictive check plot (in Results tab)
  output$predictive_plot <- renderPlot({
    analysis_data <- analysis_output()
    req(analysis_data)
    selected_prior_key <- input$plot_prior_select
    req(selected_prior_key)

    # Check if results exist for the selected prior
    if (!(selected_prior_key %in% names(analysis_data$results))) {
      plot.new()
      title("無法生成圖表：所選先驗的結果無效。")
      return()
    }

    plot_data <- analysis_data$results[[selected_prior_key]]
    # Check if predictive samples were generated (they might be NA if prior was invalid)
    if(is.null(plot_data) || any(is.na(plot_data$prior_pred_samples)) || any(is.na(plot_data$post_pred_samples))) {
      plot.new()
      title("無法生成圖表：預測檢查樣本無效。")
      return()
    }

    selected_prior_spec <- analysis_data$priors_used[[selected_prior_key]]
    observed_value <- analysis_data$theta_hat
    effect_type_label <- analysis_data$effect_type_label
    scale_lab <- analysis_data$scale_label

    df_prior_pred <- data.frame(value = plot_data$prior_pred_samples, type = "先驗預測 (Prior Predictive)")
    df_post_pred <- data.frame(value = plot_data$post_pred_samples, type = "後驗預測 (Posterior Predictive)")
    df_combined <- rbind(df_prior_pred, df_post_pred)

    gg <- ggplot(df_combined, aes(x = value, fill = type)) +
      geom_density(alpha = 0.6) +
      geom_vline(aes(xintercept = observed_value), color = "red", linetype = "dashed", linewidth = 1) +
      scale_fill_manual(name = "預測分佈類型", values = c("先驗預測 (Prior Predictive)" = "skyblue", "後驗預測 (Posterior Predictive)" = "lightcoral")) +
      labs(
        title = sprintf("預測檢查圖：%s (使用 %s 先驗)",
                        effect_type_label, selected_prior_spec$name),
        x = sprintf("模擬的效果量估計值 (%s)", scale_lab),
        y = "密度"
      ) +
      annotate("text", x = observed_value, y = 0, label = sprintf(" 觀察到的\n 估計值\n(%.3f)", observed_value), vjust = 1.1, hjust = -0.1, color = "red", size = 3.5) +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")

    print(gg)
  })

  # Display predictive check interpretation (in Results tab)
  output$plot_interpretation <- renderPrint({
    analysis_data <- analysis_output()
    req(analysis_data)
    scale_lab <- analysis_data$scale_label
    observed_value <- analysis_data$theta_hat

    cat("圖表判讀指南:\n")
    cat(sprintf("- 圖中曲線顯示了根據模型預期會觀察到的效果量估計值的分佈 (%s)。\n", scale_lab))
    cat(sprintf("  - 藍色 (先驗預測): 如果僅根據所選的『先驗分配』為真，重複實驗可能得到的結果分佈 (%s)。\n", scale_lab))
    cat(sprintf("  - 紅色 (後驗預測): 如果結合了『先驗與數據』得到的『後驗分配』為真，重複實驗可能得到的結果分佈 (%s)。\n", scale_lab))
    cat(sprintf("- 紅色虛線是從您輸入數據計算/轉換出的『觀察到的效果量估計值』(%.3f on %s)。\n\n", observed_value, scale_lab))
    cat("判讀重點:\n")
    cat(sprintf("1. 先驗預測檢查 (藍色曲線 vs 紅線): 您的觀察值 (%.3f) 是否落在先驗預期範圍內？若觀察值遠離藍色曲線的中心區域，表示先驗假設與實際觀察到的數據可能不太一致。\n", observed_value))
    cat(sprintf("2. 後驗預測檢查 (紅色曲線 vs 紅線): 您的觀察值 (%.3f) 是否落在後驗預期範圍內？理想情況下，觀察值（紅線）應大致落在後驗預測分佈（紅色曲線）的中心區域。如果觀察值落在紅色曲線的尾部，可能表示模型擬合不佳或數據中存在模型未捕捉到的異常情況。\n", observed_value))
  })
} # End server function

# --- 執行 Shiny App ---
shinyApp(ui = ui, server = server)
