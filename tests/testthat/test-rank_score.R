
rank_score_riqite <- function( n, method.list = list(name = "Wilcoxon") ){
  
  if ( method.list$name == "DIM" ) {
    stop( "Can't calculate rank scores for DIM" )
  }
  
  if(method.list$name == "Wilcoxon"){
    score = c(1:n)
    score = score/max(score)
    return(score)
  }
  
  if(method.list$name == "Stephenson"){
    score = choose( c(1:n) - 1, method.list$s - 1 )
    score = score/max(score)
    return(score)
  }
  
}

# non-scaled one and non-standardized one should match with RIQITE's result

test_that("rank score", {
  
  n = sample(10^3, 1)
  r = sample(sample(10^3, 1))
  
  method.list.riqite.wil = list(name = "Wilcoxon")
  method.list.riqite.ste = list(name = "Stephenson",
                                s = r)
  
  method.list.new.wil = list(name = "Wilcoxon",
                             scale = FALSE)
  method.list.new.ste = list(name = "Stephenson",
                             scale = FALSE,
                             s = r)
  
  expect_equal( rank_score_riqite(n, method.list.riqite.wil), rank_score(n, method.list.new.wil)) 
  expect_equal( rank_score_riqite(n, method.list.riqite.ste), rank_score(n, method.list.new.ste)) 
})
