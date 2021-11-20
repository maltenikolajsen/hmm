#' Plotting Hidden Markov Models.
#'
#'
#' @param model Class `hmm` model.
#' @param legend_position Position of legend in plots.
#' @param xlab Text for x-axis. Default 'Time'.
#' @param ylab Text for y-axis. Default 'Observed emissions'.
#' @param ...
#'
#' @return Yields four plots stored in a list:
#'
#' |              |                                                                   |
#' |--------------|-------------------------------------------------------------------|
#' | `p1`         | Plot of observed emissions by most likely hidden state (Global)   |
#' | `p2`         | Plot of observed emissions by most likely hidden state (Local)    |
#' | `p3`         | Plot of probs. used in Viterbi algorithm (Global)                 |
#' | `p4`         | Plot of state probabilities (Local)                               |
#'
#' @examples
#' ## Continuation of Earthquake example
#' \dontshow{example(hmm)}
#' plot(hmm.EQ)
#'
#' @export
#'
plot.hmm <- function(model,
                     legend_position = "topright",
                     xlab = "Time",
                     ylab = "Observed emissions",
                     cols = 1:model$m,
                     ...){
  m <- model$m; n <- model$n

  #Plot global decoding.
  if(!is.null(model$viterbi_s)){
    t <- 1:n
    p1 <- {plot(t,
                model$x, main = "Global decoding",
                xlab = xlab,
                ylab = ylab,
                ...)
      for(i in 1:m){
        sub_x <- model$x[model$viterbi_s == i]
        sub_t <- t[model$viterbi_s == i]
        if(hasArg('type') && list(...)$type %in% c('h', 'l')){
          lines(sub_t, sub_x, col = cols[i], ...)
        }
        else{
          points(sub_t, sub_x, col = cols[i], pch = 16)
        }
      }
      legend(legend_position, legend=1:m, pch=rep(16, m), col = cols, title="State")
    }
  }

  #Plot local decoding.
  if(!is.null(model$posterior_s)){
    t <- 1:n
    p2 <- {plot(t,
                model$x, main = "Posterior decoding",
                xlab = xlab,
                ylab = ylab,
                ...)
      for(i in 1:m){
        sub_x <- model$x[model$posterior_s == i]
        sub_t <- t[model$posterior_s == i]
        if(hasArg('type') && list(...)$type %in% c('h', 'l')){
          lines(sub_t, sub_x, col = cols[i], ...)
        }
        else{
          points(sub_t, sub_x, col = cols[i], pch = 16)
        }
      }
      legend(legend_position, legend=1:m, pch=rep(16, m), col = cols, title="State")
    }
  }

  #Stacked bar plot for global decoding probabilities.
  if(!is.null(model$viterbi_p)){
    tmp <- as.table(model$viterbi_p)
    colnames(tmp) <- 1:n
    rownames(tmp) <- 1:m
    p3 <- {
      barplot(tmp, col = cols, main = "Global probabilities", xlab = xlab, ylab = "Probability of observing emission")
      legend("topright", legend=1:m, pch=rep(16, m), col = cols, title="State")
    }
  }

  #Stacked bar plot for local decoding probabilities.
  if(!is.null(model$posterior_p)){
    tmp <- as.table(model$posterior_p)
    colnames(tmp) <- 1:n
    rownames(tmp) <- 1:m
    p4 <- {
      barplot(tmp,
              col = cols,
              main = "Posterior probabilities",
              xlab = xlab,
              ylab = "Probability of observing emission")
      legend("topright", legend=1:m, pch=rep(16, m), col = cols, title="State")
    }
  }

  invisible()
}
