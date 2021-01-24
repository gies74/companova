run <- function(file) {
  data <- p_readFileAsTable(file)

  vars <- p_auxiliaryVars(data)

  p_reportStatistics(vars)

}



p_readFileAsTable <- function(file) {
  rawData <- read.csv2(file, sep = ' ', header=F)

  # input validation
  if (length(rawData) != 8) {
    stop(paste("Each line of the input file should contain exactly 8 space-delimited columns. ",length(rawData),"were encountered"))
  }


  len <- nrow(rawData)
  nItems <- ceiling(sqrt(len))
  if (len < 6 | nItems * (nItems - 1) != len) {
    stop(paste("The number of lines (", len, ") in the input file should match nItems * ( nItems - 1 ). For best approximation nItems=", nItems, "this currently evaluates to", nItems * (nItems - 1)))
  }


  lineNum <- 1
  for (i in 1:(nItems - 1)) {
    for (j in (i+1):nItems) {
      if (rawData$V1[lineNum] != i) {
        stop(paste("At column 1 of line", lineNum, "item id", i, "was expected, but", rawData$V1[lineNum],"was encountered."))
      }
      if(rawData$V1[lineNum + 1] != j) {
        stop(paste("At column 1 of line", lineNum+1, "item id", j, "was expected, but", rawData$V1[lineNum+1],"was encountered."))
      }
      lineNum <- lineNum + 2
    }
  }

  nJudges <- sum(rawData[1, 2:8])
  for (i in 2:len) {
    if (sum(rawData[i, 2:8]) != nJudges) {
      stop(paste("Total number of listeners expected to be",nJudges,", like sum of columns 2...8 at line 1. However, at line",i,sum(rawData[i,2:8]),"were encountered."))
    }
  }

  return(list(rawData, nItems, nJudges))

}

p_auxiliaryVars <- function(data) {
  table <- data[[1]]
  for (i in 1:nrow(table)) {
    table$V9[i] <- -3 * table$V2[i] - 2 * table$V3[i] - 1 * table$V4[i] + 1 * table$V6[i] + 2 * table$V7[i] + 3 * table$V8[i]
    table$V10[i] <- table$V9[i] / nrow(table)
    if (i %% 2 == 0) {
      table$V11[i - 1] <- +.5 * (table$V10[i - 1] - table$V10[i])
      table$V11[i]     <- -.5 * (table$V10[i - 1] - table$V10[i])
    }
  }

  nItems <- data[[2]]
  a <- list(nItems)
  for (i in 1:nItems) {
    a[i] <- 0.0
    for (j in 1:nrow(table)) {
      if (table$V1[j] == i) {
        a[i] <- a[[i]] + table$V11[j]
      }
    }
    a[i] <- a[[i]] / nItems
  }

  p <- list()
  for (i in 0:(nrow(table)/2 - 1)) {
    p[i+1] <- table$V11[i*2 + 1]
  }

  st <- 3^2 * (sum(table$V2) + sum(table$V8)) + 2^2 * (sum(table$V3) + sum(table$V7)) + sum(table$V4) + sum(table$V6)

  sm <- 0
  for (i in 1:nrow(table)) {
    sm <- sm + table$V10[i]^2
  }
  sm <- sm * nrow(table)

  sa <- 0
  for (i in 1:nItems) {
    sa <- sa + a[[i]] ^ 2
  }
  nListeners <- data[[3]]
  sa <- 2 * nListeners * nItems * sa

  sp <- 0
  for (i in 1:(nrow(table)/2)) {
    sp <- sp + p[[i]] ^2
  }
  sp <- 2 * nrow(table) * sp

  return (list(table, nItems, nListeners, a, p, st, sm, sa, sp))
}

p_reportStatistics <- function(vars) {

  SA <- vars[[8]]
  SM <- vars[[7]]
  SP <- vars[[9]]
  SD <- SM - SP
  SG <- SP - SA
  M  <- vars[[2]]*(vars[[2]] - 1)/2

  ST <- vars[[6]]
  SE <- ST - SM

  DFA <- vars[[2]] - 1
  DFG <- M - vars[[2]]+1
  DFP <- M
  DFD <- M
  DFM <- 2 * M
  DFE <- 2 * M * (vars[[3]] - 1)
  DFT <- 2 * M * nrow(vars[[1]])
  QSA <- SA / DFA
  QSG <- SG / DFG
  QSD <- SD / DFD
  QSE <- SE / DFE

  FM <- QSA / QSE
  FG <- QSG / QSE
  FD <- QSD / QSE

  Y  <- sqrt(QSE / (2 * M * nrow(vars[[1]])))


  FORMAT <- "%-30s %10s %6s %8s %12s\n"

  cat(sprintf("%50s\n\n", "ANALYSIS OF VARIANCE"))

  cat(sprintf(FORMAT, "SOURCE", "SS", "DF", "MS", "F-RATIO"))
  cat(sprintf(FORMAT, "MAIN EFFECTS", sprintf("%.2f", SA), sprintf("%d", DFA), sprintf("%.2f", QSA), sprintf("%.2f", FM)))
  cat(sprintf(FORMAT, "DEVIAT. FROM SUBTRACTIVITY", sprintf("%.2f", SG), sprintf("%d", DFG), sprintf("%.2f", QSG), sprintf("%.2f", FG)))
  cat(sprintf(FORMAT, "AVERAGE PREFERENCES", sprintf("%.2f", SP), sprintf("%d", DFP), "", ""))
  cat(sprintf(FORMAT, "ORDER EFFECTS", sprintf("%.2f", SD), sprintf("%d", DFD), sprintf("%.2f", QSD), sprintf("%.2f", FD)))
  cat(sprintf(FORMAT, "MEANS", sprintf("%.2f", SM), sprintf("%d", DFM), "", ""))
  cat(sprintf(FORMAT, "ERROR", sprintf("%.2f", SE), sprintf("%d", DFE), sprintf("%.2f", QSE), ""))
  cat(sprintf(FORMAT, "TOTAL", sprintf("%.2f", ST), sprintf("%d", DFT), "", ""))
  cat(sprintf("\n\n---------------------------------\n\n\n"))

  for (i in 1:vars[[2]]) {
    cat(sprintf("%s = %8s\n", sprintf("A%d", i), sprintf("%.5f", vars[[4]][i])))
  }

  cat(sprintf("\n\n"))
  cat(sprintf("YARDSTICK Y = Q * %.4f\n\n", Y))

  cat(sprintf("%3s", ""))
  for (i in 1:vars[[2]]) {
    cat(sprintf("%7s", sprintf("A%d", i)))
  }
  cat("\n")

  for (i in 1:vars[[2]]) {
    cat(sprintf("%3s", sprintf("A%d", i)))
    for (j in 1:vars[[2]]) {
      if (i != j) {
        mem <- vars[[4]][[i]] - vars[[4]][[j]]
        cat(sprintf("%7s", sprintf("%.3f", abs(mem))))
      } else {
        cat(sprintf("%7s", "-"))
      }
    }
    cat("\n")
  }
}
