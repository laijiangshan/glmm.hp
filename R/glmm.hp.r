#' Hierarchical Partitioning of Marginal R2 for Generalized Mixed-Effect Models

#' @param  mod  Fitted lme4 or nlme model objects.

#' @details This function conducts hierarchical partitioning to calculate the individual contributions of each predictor towards marginal R2 for Generalized Mixed-effect Model. The marginal R2 is the output of r.squaredGLMM in MuMIn package.

#' @return \item{Total.Marginal.R2}{The marginal R2 (fixed effect) for the full model.}
#' @return \item{Hier.part}{A matrix containing individual effects and percentage of individual effects towards total marginal R2 for each predictor.}

#' @author {Jiangshan Lai} \email{lai@ibcas.ac.cn}
#' @author {Kim Nimon} \email{kim.nimon@gmail.com}


#' @references
#' \itemize{
#' \item Lai J.,Zou Y., Zhang J.,Peres-Neto P.(2022) Generalizing hierarchical and variation partitioning in multiple regression and canonical analyses using the rdacca.hp R package.Methods in Ecology and Evolution,13(4):782-788<DOI:10.1111/2041-210X.13800>
#' \item Chevan, A. & Sutherland, M. (1991). Hierarchical partitioning. American Statistician, 45, 90-96. doi:10.1080/00031305.1991.10475776
#' \item Nimon, K., Oswald, F.L. & Roberts, J.K. (2013). Yhat: Interpreting regression effects. R package version 2.0.0.
#' \item Nakagawa, S., & Schielzeth, H. (2013). A general and simple method for obtaining R2 from generalized linear mixed-effects models. Methods in Ecology and Evolution, 4(2), 133-142.
#' \item Nakagawa, S., Johnson, P. C., & Schielzeth, H. (2017). The coefficient of determination R2 and intra-class correlation coefficient from generalized linear mixed-effects models revisited and expanded. Journal of the Royal Society Interface, 14(134), 20170213.
#' }


#'@export
#'@examples
#'library(MuMIn)
#'library(lme4)
#'mod1 <- lmer(Sepal.Length ~ Petal.Length + Petal.Width +(1 | Species), data = iris)
#'r.squaredGLMM(mod1)
#'glmm.hp(mod1)
#'plot(glmm.hp(mod1))

glmm.hp <- function(mod)
{
  # initial checks
  if (!inherits(mod, c("merMod","lme"))) stop("glmm.hp only supports lme or merMod objects at the moment")
  if(inherits(mod, "merMod"))
  {# interaction checks
  if("*"%in%strsplit(as.character(mod@call$formula)[3],"")[[1]])stop("glmm.hp does not supports interaction terms at the moment")
  varname <- strsplit(strsplit(as.character(mod@call$formula)[3],"(",fixed=T)[[1]][1]," ")[[1]]
  ivname <- varname[seq(1,length(varname),2)]
  }
  
  if(inherits(mod, "lme"))
  {# interaction checks
  if("*"%in%strsplit(as.character(mod$call$fixed)[3],"")[[1]])stop("glmm.hp does not supports interaction terms at the moment")
  ivname <- strsplit(as.character(mod$call$fixed)[3]," + ",fixed=T)[[1]]
  }
  
  iv.name <- ivname
  nvar <- length(iv.name)
  if (nvar < 2)
    stop("Analysis not conducted. Insufficient number of predictors.")

  totalN <- 2^nvar - 1
  binarymx <- matrix(0, nvar, totalN)
  for (i in 1:totalN) {
    binarymx <- creatbin(i, binarymx)
  }

outr2  <- r.squaredGLMM(mod)
r2type  <-  row.names(outr2)
nr2type   <-  length(r2type)
if(nr2type==0)
{nr2type <- 1
r2type <- 'hierarchical.partitioning'
}
#ifelse(class(mod)=="merMod",dat <- eval(mod@call$data),dat <- eval(mod$call$data))
if(inherits(mod, "merMod"))
{dat <- na.omit(eval(mod@call$data))
if(!inherits(dat, "data.frame")){stop("Please change the name of data object in the original (g)lmm analysis then try again.")}
to_del <- paste(paste("-", iv.name, sep= ""), collapse = " ")
# reduced formula
 modnull<- stats::update(stats::formula(mod), paste(". ~ . ", to_del, sep=""))
 mod_null <-  stats::update(object = mod, formula. = modnull, data = dat)
 }

if(inherits(mod, "lme"))
{dat <- eval(mod$call$data)
mod_null <- stats::update(object = mod,data=dat,fixed=~1)
}

outputList  <- list()
outputList[[1]] <- outr2
for (k in 1:nr2type)
{
  commonM <- matrix(nrow = totalN, ncol = 3)
  for (i in 1:totalN) {
    tmp.name <- iv.name[as.logical(binarymx[, i])]
   if(inherits(mod, "merMod"))
   {to_add <- paste(paste("+", tmp.name, sep= ""), collapse = " ")
    modname <- stats::update(stats::formula(mod_null), paste(". ~ . ", to_add, sep=""))
    modnew  <- stats::update(object = mod_null, formula. = modname, data = dat) 
	}
	
	if(inherits(mod, "lme"))
	{to_add <- paste("~",paste(tmp.name,collapse = " + "),sep=" ")
	  modnew  <- stats::update(object = mod_null, data = dat,fixed=to_add) 
	}
	
	commonM[i, 2]  <- MuMIn::r.squaredGLMM(modnew)[k,1]
  }

  commonlist <- vector("list", totalN)

  seqID <- vector()
  for (i in 1:nvar) {
    seqID[i] = 2^(i-1)
  }


  for (i in 1:totalN) {
    bit <- binarymx[1, i]
    if (bit == 1)
      ivname <- c(0, -seqID[1])
    else ivname <- seqID[1]
    for (j in 2:nvar) {
      bit <- binarymx[j, i]
      if (bit == 1) {
        alist <- ivname
        blist <- genList(ivname, -seqID[j])
        ivname <- c(alist, blist)
      }
      else ivname <- genList(ivname, seqID[j])
    }
    ivname <- ivname * -1
    commonlist[[i]] <- ivname
  }

  for (i in 1:totalN) {
    r2list <- unlist(commonlist[i])
    numlist  <-  length(r2list)
    ccsum  <-  0
    for (j in 1:numlist) {
      indexs  <-  r2list[[j]]
      indexu  <-  abs(indexs)
      if (indexu != 0) {
        ccvalue  <-  commonM[indexu, 2]
        if (indexs < 0)
          ccvalue  <-  ccvalue * -1
        ccsum  <-  ccsum + ccvalue
      }
    }
    commonM[i, 3]  <-  ccsum
  }

  orderList <- vector("list", totalN)
  index  <-  0
  for (i in 1:nvar) {
    for (j in 1:totalN) {
      nbits  <-  sum(binarymx[, j])
      if (nbits == i) {
        index  <-  index + 1
        commonM[index, 1] <- j
      }
    }
  }

  outputcommonM <- matrix(nrow = totalN + 1, ncol = 2)
  totalRSquare <- sum(commonM[, 3])
  for (i in 1:totalN) {
    outputcommonM[i, 1] <- round(commonM[commonM[i,
                                                 1], 3], digits = 4)
    outputcommonM[i, 2] <- round((commonM[commonM[i,
                                                  1], 3]/totalRSquare) * 100, digits = 2)
  }
  outputcommonM[totalN + 1, 1] <- round(totalRSquare,
                                        digits = 4)
  outputcommonM[totalN + 1, 2] <- round(100, digits = 4)
  rowNames  <-  NULL
  for (i in 1:totalN) {
    ii  <-  commonM[i, 1]
    nbits  <-  sum(binarymx[, ii])
    cbits  <-  0
    if (nbits == 1)
      rowName  <-  "Unique to "
    else rowName  <-  "Common to "
    for (j in 1:nvar) {
      if (binarymx[j, ii] == 1) {
        if (nbits == 1)
          rowName  <-  paste(rowName, iv.name[j], sep = "")
        else {
          cbits  <-  cbits + 1
          if (cbits == nbits) {
            rowName  <-  paste(rowName, "and ", sep = "")
            rowName  <-  paste(rowName, iv.name[j], sep = "")
          }
          else {
            rowName  <-  paste(rowName, iv.name[j], sep = "")
            rowName  <-  paste(rowName, ", ", sep = "")
          }
        }
      }
    }
    rowNames  <-  c(rowNames, rowName)
  }
  rowNames  <-  c(rowNames, "Total")
  rowNames <- format.default(rowNames, justify = "left")
  colNames <- format.default(c("Fractions", " % Total"),
                             justify = "right")
  dimnames(outputcommonM) <- list(rowNames, colNames)

 # VariableImportance <- matrix(nrow = nvar, ncol = 4)
 VariableImportance <- matrix(nrow = nvar, ncol = 2)
  for (i in 1:nvar) {
	#VariableImportance[i, 3] <-  round(sum(binarymx[i, ] * (commonM[,3]/apply(binarymx,2,sum))), digits = 4)
	VariableImportance[i, 1] <-  round(sum(binarymx[i, ] * (commonM[,3]/apply(binarymx,2,sum))), digits = 4)
  }
  
  #VariableImportance[,1] <- outputcommonM[1:nvar,1]
 # VariableImportance[,2] <- VariableImportance[,3]-VariableImportance[,1]
  
  #total=round(sum(VariableImportance[,3]),digits = 4)
  total=round(sum(VariableImportance[,1]),digits = 4)
  #VariableImportance[, 4] <- round(100*VariableImportance[, 3]/total,2)
VariableImportance[, 2] <- round(100*VariableImportance[, 1]/total,2)
  dimnames(VariableImportance) <- list(iv.name, c("Individual","I.perc(%)"))

#if(commonality)
#{outputList <- list(Total.Marginal.R2=total,commonality = outputcommonM, Hier.part = VariableImportance)}
#else
#{outputList<-list(Total.Marginal.R2=total,Hier.part= VariableImportance)}
#outputList[[j]]<- list(Total.Marginal.R2=total,Hier.part= VariableImportance)
outputList[[k+1]]<-VariableImportance
}
names(outputList) <- c("r.squaredGLMM",r2type)
class(outputList) <- "glmmhp" # Class definition
outputList
}

