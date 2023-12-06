two.list.overlap <- function(one,two) { 
  #input is in form of a ordered gene list (character vector). Returns genelists of various intersections!
  #READ f[which(!f %in% g)]  AS "the values of f that are not in g"
  #based on VennDiagram::calculate.overlap
  n12 <- intersect(one, two)
  a1 = n12
  a2 = one[which(!one %in% a1)]
  a3 = two[which(!two %in% a1)]
  overlap <- list('uniqueone&two' = a1,'uniqueone' = a2, 'uniquetwo' = a3)
}

bi.venn <- function(x, colors = NULL, labs, title=NULL) { #the output of three.list.overlap can vary as the fit changes slightly with each iteration, may need to run multiple times to get accurate plot.
  require(eulerr)
  set.seed(1)
  y <- c(one = length(x$uniqueone), two = length(x$uniquetwo), "one&two" = length(x$`uniqueone&two`))
  fit.circ <- euler(y)
  fit.ellp <- euler(y, shape = "ellipse")
  
  if (missing(labs)) {
    if (fit.circ$stress > fit.ellp$stress ) {
      plot(fit.ellp, quantities = TRUE, fills = colors, main = title)
    } else {
      plot(fit.circ, quantities = TRUE, fills = colors, main = title)
    }
  } else {
    if (fit.circ$stress > fit.ellp$stress ) {
      plot(fit.ellp, quantities = TRUE, labels = labs, fills = colors, main = title)
    } else {
      plot(fit.circ, quantities = TRUE, labels = labs, fills = colors, main = title)
    }
  }
}

three.list.overlap <- function(one,two,three) { 
  #input is in form of a ordered gene list (character vector). Returns genelists of various intersections!
  #READ f[which(!f %in% g)]  AS "the values of f that are not in g"
  #based on VennDiagram::calculate.overlap
  n12 <- intersect(one, two)
  n13 <- intersect(one, three)
  n23 <- intersect(two, three)
  n123 <- intersect(n12, three)
  a1 = n123
  a2 = n12[which(!n12 %in% a1)]
  a3 = n13[which(!n13 %in% a1)]
  a4 = n23[which(!n23 %in% a1)]
  a5 = one[which(!one %in% c(a1, a2, a3))]
  a6 = two[which(!two %in% c(a1, a2, a4))]
  a7 = three[which(!three %in% c(a1, a3, a4))]
  overlap <- list('uniqueone&two&three' = a1, 'uniqueone&two' = a2, 'uniqueone&three' = a3, 
                  'uniquetwo&three' = a4, 'uniqueone' = a5, 'uniquetwo' = a6, 'uniquethree' = a7)
}

tri.venn <- function(x, colors = NULL, labs, title=NULL) { #the output of three.list.overlap can vary as the fit changes slightly with each iteration, may need to run multiple times to get accurate plot.
  require(eulerr)
  set.seed(1)
  y <- c(one = length(x$uniqueone), two = length(x$uniquetwo), three = length(x$uniquethree),
         "one&two" = length(x$`uniqueone&two`), "one&three" = length(x$`uniqueone&three`),
         "two&three" = length(x$`uniquetwo&three`), "one&two&three" = length(x$`uniqueone&two&three`))
  fit.circ <- euler(y)
  fit.ellp <- euler(y, shape = "ellipse")
  
  if (missing(labs)) {
    if (fit.circ$stress > fit.ellp$stress ) {
      plot(fit.ellp, quantities = TRUE, fills = colors, main = title)
    } else {
      plot(fit.circ, quantities = TRUE, fills = colors, main = title)
    }
  } else {
    if (fit.circ$stress > fit.ellp$stress ) {
      plot(fit.ellp, quantities = TRUE, labels = labs, fills = colors, main = title)
    } else {
      plot(fit.circ, quantities = TRUE, labels = labs, fills = colors, main = title)
    }
  }
}

four.list.overlap <- function(one,two,three,four) { 
  #input is in form of a ordered gene list (character vector). Returns genelists of various intersections!
  #READ f[which(!f %in% g)]  AS "the values of f that are not in g"
  #based on VennDiagram::calculate.overlap
  n12 <- intersect(one, two)
  n13 <- intersect(one, three)
  n14 <- intersect(one, four)
  n23 <- intersect(two, three)
  n24 <- intersect(two, four)
  n34 <- intersect(three, four)
  n123 <- intersect(n12, three)
  n124 <- intersect(n12, four)
  n134 <- intersect(n13, four)
  n234 <- intersect(n23, four)
  n1234 <- intersect(n123, four)
  a6 = n1234
  a12 = n123[which(!n123 %in% a6)]
  a11 = n124[which(!n124 %in% a6)]
  a5 = n134[which(!n134 %in% a6)]
  a7 = n234[which(!n234 %in% a6)]
  a15 = n12[which(!n12 %in% c(a6, a11, a12))]
  a4 = n13[which(!n13 %in% c(a6, a5, a12))]
  a10 = n14[which(!n14 %in% c(a6, a5, a11))]
  a13 = n23[which(!n23 %in% c(a6, a7, a12))]
  a8 = n24[which(!n24 %in% c(a6, a7, a11))]
  a2 = n34[which(!n34 %in% c(a6, a5, a7))]
  a9 = one[which(!one %in% c(a4, a5, a6, a10, a11, a12, a15))]
  a14 = two[which(!two %in% c(a6, a7, a8, a11, a12, a13, a15))]
  a1 = three[which(!three %in% c(a2, a4, a5, a6, a7, a12, a13))]
  a3 = four[which(!four %in% c(a2, a5, a6, a7, a8, a10, a11))]
  overlap <- list('one&two&three&four' = a6, 'uniqueone&two&three' = a12, 'uniqueone&two&four' = a11, 'uniqueone&three&four' = a5, 
                  'uniquetwo&three&four' = a7, 'uniqueone&two' = a15, 'uniqueone&three' = a4, 'uniqueone&four' = a10, 'uniquetwo&three' = a13, 
                  'uniquetwo&four' = a8, 'uniquethree&four' = a2, 'uniqueone' = a9, 'uniquetwo' = a14, 'uniquethree' = a1, 'uniquefour' = a3)
}

quad.venn <- function(x, colors = NULL, labs, title=NULL) { #the output of four.list.overlap can vary as the fit changes slightly with each iteration, may need to run multiple times to get accurate plot.
  require(eulerr)
  set.seed(1)
  y <- c(one = length(x$uniqueone), two = length(x$uniquetwo), three = length(x$uniquethree), four = length(x$uniquefour),
         "one&two" = length(x$`uniqueone&two`), "one&three" = length(x$`uniqueone&three`), "one&four" = length(x$`uniqueone&four`),
         "two&three" = length(x$`uniquetwo&three`), "two&four" = length(x$`uniquetwo&four`), "three&four" = length(x$`uniquethree&four`),
         "one&two&three" = length(x$`uniqueone&two&three`), "one&two&four" = length(x$`uniqueone&two&four`), 
         "one&three&four" = length(x$`uniqueone&three&four`), "two&three&four" = length(x$`uniquetwo&three&four`), "one&two&three&four" = length(x$`one&two&three&four`))
  fit.circ <- euler(y)
  fit.ellp <- euler(y, shape = "ellipse")
  
  if (missing(labs)) {
    if (fit.circ$stress > fit.ellp$stress ) {
      plot(fit.ellp, quantities = TRUE, fills = colors, main = title)
    } else {
      plot(fit.circ, quantities = TRUE, fills = colors, main = title)
    }
  } else {
    if (fit.circ$stress > fit.ellp$stress ) {
      plot(fit.ellp, quantities = TRUE, labels = labs, fills = colors, main = title)
    } else {
      plot(fit.circ, quantities = TRUE, labels = labs, fills = colors, main = title)
    }
  }
}

five.list.overlap <- function(one,two,three,four,five) { 
  #input is in form of a ordered gene list (character vector). Returns genelists of various intersections!
  #READ f[which(!f %in% g)]  AS "the values of f that are not in g"
  #based on VennDiagram::calculate.overlap
  n12 <- intersect(one, two)
  n13 <- intersect(one, three)
  n14 <- intersect(one, four)
  n23 <- intersect(two, three)
  n24 <- intersect(two, four)
  n34 <- intersect(three, four)
  n15 <- intersect(one, five)
  n25 <- intersect(two, five)
  n35 <- intersect(three, five)
  n45 <- intersect(four, five)
  n123 <- intersect(n12, three)
  n124 <- intersect(n12, four)
  n134 <- intersect(n13, four)
  n234 <- intersect(n23, four)
  n125 <- intersect(n12, five)
  n135 <- intersect(n13, five)
  n145 <- intersect(n14, five)
  n235 <- intersect(n23, five)
  n245 <- intersect(n24, five)
  n345 <- intersect(n34, five)
  n1234 <- intersect(n123, four)
  n1235 <- intersect(n123, five)
  n1245 <- intersect(n124, five)
  n1345 <- intersect(n134, five)
  n2345 <- intersect(n234, five)
  n12345 <- intersect(n1234, five)
  a1 = n12345
  a2 = n1234[which(!n1234 %in% a1)]
  a3 = n1235[which(!n1235 %in% a1)]
  a4 = n1245[which(!n1245 %in% a1)]
  a5 = n1345[which(!n1345 %in% a1)]
  a6 = n2345[which(!n2345 %in% a1)]
  a7 = n123[which(!n123 %in% c(a1,a2,a3))]
  a8 = n124[which(!n124 %in% c(a1,a2,a4))]
  a9 = n125[which(!n125 %in% c(a1,a3,a4))]
  a10 = n134[which(!n134 %in% c(a1,a2,a5))]
  a11 = n135[which(!n135 %in% c(a1,a3,a5))]
  a27 = n145[which(!n145 %in% c(a1,a4,a5))]
  a12 = n234[which(!n234 %in% c(a1,a2,a6))]
  a13 = n235[which(!n235 %in% c(a1,a3,a6))]
  a14 = n245[which(!n245 %in% c(a1,a4,a6))]
  a15 = n345[which(!n345 %in% c(a1,a5,a6))]
  a16 = n12[which(!n12 %in% c(a1, a2, a3, a4, a7, a8, a9))]
  a17 = n13[which(!n13 %in% c(a1, a3, a3, a5, a7, a10, a11))]
  a18 = n14[which(!n14 %in% c(a1, a2, a4, a5, a8, a10, a27))]
  a19 = n15[which(!n15 %in% c(a1, a3, a4, a5, a9, a11, a27))]
  a20 = n23[which(!n23 %in% c(a1, a2, a3, a6, a7, a12, a13))]
  a21 = n24[which(!n24 %in% c(a1, a2, a4, a6, a8, a12, a14))]
  a22 = n25[which(!n25 %in% c(a1, a3, a4, a6, a9, a13, a14))]
  a23 = n34[which(!n34 %in% c(a1, a2, a5, a6, a10, a12, a15))]
  a24 = n35[which(!n35 %in% c(a1, a3, a5, a6, a11, a13, a15))]
  a25 = n45[which(!n45 %in% c(a1, a4, a5, a6, a14, a15, a27))]
  a26 = one[which(!one %in% c(a1, a2, a3, a4, a5, a7, a8, a9, a10, a11, a16, a17, a18, a19, a27))]
  a28 = two[which(!two %in% c(a1, a2, a3, a4, a6, a7, a8, a9, a12, a13, a14, a16, a20, a21, a22))]
  a29 = three[which(!three %in% c(a1, a2, a3, a5, a6, a7, a10, a11, a12, a13, a15, a17, a20, a23, a24))]
  a30 = four[which(!four %in% c(a1, a2, a4, a5, a6, a8, a10, a12, a14, a15, a18, a21, a23, a25, a27))]
  a31 = five[which(!five %in% c(a1, a3, a4, a5, a6, a9, a11, a27, a13, a14, a15, a19, a22, a24, a25))]
  overlap1 <- list('one&two&three&four&five' = a1, 'uniqueone&two&three&four' = a2, 'uniqueone&two&three&five' = a3, 'uniqueone&two&four&five'= a4, 
                   'uniqueone&three&four&five' = a5, 'uniquetwo&three&four&five' = a6, 'uniqueone&two&three' = a7, 'uniqueone&two&four' = a8, 'uniqueone&two&five' = a9, 
                   'uniqueone&three&four' = a10, 'uniqueone&three&five' = a11, 'uniquetwo&three&four' = a12, 'uniquetwo&three&five' = a13, 'uniquetwo&four&five' = a14, 
                   'uniquethree&four&five' = a15, 'uniqueone&two' = a16, 'uniqueone&three' = a17, 'uniqueone&four' = a18, 'uniqueone&five' = a19, 'uniquetwo&three' = a20,
                   'uniquetwo&four' = a21, 'uniquetwo&five' = a22, 'uniquethree&four' = a23, 'uniquethree&five' = a24, 'uniquefour&five' = a25,
                   'uniqueone' = a26,'uniquetwo' = a28, 'uniquethree' = a29, 'uniquefour' = a30, 'uniquefive' = a31
  )
}

quint.venn <- function(x, colors=NULL, labs, title=NULL) { #the output of five.list.overlap can vary as the fit changes slightly with each iteration, may need to run multiple times to get accurate plot.
  require(eulerr)
  set.seed(1)
  y <- c(one = length(x$uniqueone), two = length(x$uniquetwo), three = length(x$uniquethree), four = length(x$uniquefour), five = length(x$uniquefive),
         "one&two" = length(x$`uniqueone&two`), "one&three" = length(x$`uniqueone&three`), "one&four" = length(x$`uniqueone&four`), "two&three" = length(x$`uniquetwo&three`), 
         "two&four" = length(x$`uniquetwo&four`), "three&four" = length(x$`uniquethree&four`), "one&five" = length(x$`uniqueone&five`), "two&five" = length(x$`uniquetwo&five`), 
         "three&five" = length(x$`uniquethree&five`), "four&five" = length(x$`uniquefour&five`), "one&two&three" = length(x$`uniqueone&two&three`), 
         "one&two&four" = length(x$`uniqueone&two&four`), "one&two&five" = length(x$`uniqueone&two&five`), "one&three&four" = length(x$`uniqueone&three&four`), 
         "one&three&five" = length(x$`uniqueone&three&five`), "two&three&four" = length(x$`uniquetwo&three&four`), "two&three&five" = length(x$`uniquetwo&three&five`), 
         "two&four&five" = length(x$`uniquetwo&four&five`), "three&four&five" = length(x$`uniquethree&four&five`),'one&two&three&four&five' = length(x$`one&two&three&four&five`), 
         'one&two&three&five' = length(x$`uniqueone&two&three&five`), 'one&two&four&five'= length(x$`uniqueone&two&four&five`), 'one&three&four&five' = length(x$`uniqueone&three&four&five`), 
         "two&three&four&five" = length(x$`uniquetwo&three&four&five`),"one&two&three&four" = length(x$`uniqueone&two&three&four`))
  fit.circ <- euler(y)
  fit.ellp <- euler(y, shape = "ellipse")
  
  if (missing(labs)) {
    if (fit.circ$stress > fit.ellp$stress ) {
      plot(fit.ellp, quantities = TRUE, fills = colors, main = title)
    } else {
      plot(fit.circ, quantities = TRUE, fills = colors, main - title)
    }
  } else {
    if (fit.circ$stress > fit.ellp$stress ) {
      plot(fit.ellp, quantities = TRUE, labels = labs, fills = colors, main = title)
    } else {
      plot(fit.circ, quantities = TRUE, labels = labs, fills = colors, main = title)
    }
  }
}