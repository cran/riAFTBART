RF2List <-
  function(rf){
    treeList <- NULL
    treeList$ntree <- rf$ntree
    treeList$list <- vector("list",rf$ntree)
    for(i in 1:treeList$ntree){
      treeList$list[[i]] <- RRF::getTree(rf,i,labelVar=FALSE)
    }
    return(treeList)
  }

extractRules <-
  function(treeList,X,ntree=100,maxdepth=6,random=FALSE,digits=NULL){
    if(is.numeric(digits)) digits <- as.integer(abs(digits))

    levelX = list()
    for(iX in 1:ncol(X))
      levelX <- c(levelX,list(levels(X[,iX])))
    # X <- NULL; target <- NULL
    ntree=min(treeList$ntree,ntree)
    allRulesList = list()
    for(iTree in 1:ntree){
      if(random==TRUE){max_length = sample(1:maxdepth,1,replace=FALSE)}else{
        max_length = maxdepth}
      rule = list(); count = 0; rowIx = 1;
      # tree = RRF::getTree(rf,iTree,labelVar=FALSE)
      tree <- treeList$list[[iTree]]
      if(nrow(tree)<=1) next # skip if there is no split
      ruleSet = vector("list", length(which(tree[,"status"]==-1)))
      res = treeVisit(tree,rowIx = rowIx,count,ruleSet,rule,levelX,length=0,max_length=max_length,digits=digits)
      allRulesList = c(allRulesList, res$ruleSet)
    }
    allRulesList <- allRulesList[!unlist(lapply(allRulesList, is.null))]
    message(length(allRulesList), " rules (length<=", max_length, ") were extracted from the first ", ntree, " trees.")
    rulesExec <- ruleList2Exec(X,allRulesList)
    return(rulesExec)
  }

getRuleMetric <-
  function(ruleExec, X, target){
    #typeX = getTypeX(X)
    #ruleExec <- unique(t(sapply(allRulesList,RuleList2Exec,typeX=typeX)))
    #colnames(ruleExec) <- c("len","condition")
    ruleMetric <- t(sapply(ruleExec[,"condition",drop=FALSE],measureRule,X,target))
    rownames(ruleMetric) = NULL;
    # ruleMetric <- cbind( ruleExec[,1] ,  ruleMetric )
    colnames(ruleMetric) <- c("len","freq","err","condition","pred")
    dIx <- which(ruleMetric[,"len"]=="-1")
    if(length(dIx)>0){
      ruleMetric <- ruleMetric[-dIx,]
      message(paste(length(dIx), " paths are ignored."))
    }
    return(ruleMetric)
    #qIx = order((1- as.numeric(ruleMetric[,"err"])),
    #            as.numeric(ruleMetric[,"freq"]),
    #            -as.numeric(ruleMetric[,"len"]),
    #            decreasing=TRUE)
    #return(ruleMetric[qIx,])
  }

pruneRule <-
  function(rules,X,target, maxDecay = 0.05, typeDecay = 2){
    newRuleMetric <- NULL
    for(i in 1:nrow(rules)){
      newRuleMetric <- rbind(newRuleMetric, pruneSingleRule(rules[i,],X,target, maxDecay, typeDecay))
    }
    return(newRuleMetric)
  }

selectRuleRRF <-
  function(ruleMetric,X,target){
    ruleI = sapply(ruleMetric[,"condition"],rule2Table,X,target)
    coefReg <- 0.95 - 0.01*as.numeric(ruleMetric[,"len"])/max(as.numeric(ruleMetric[,"len"]))
    rf <- RRF::RRF(ruleI,as.factor(target), flagReg = 1, coefReg=coefReg, mtry = (ncol(ruleI)*1/2) , ntree=50, maxnodes= 10,replace=FALSE)
    imp <- rf$importance/max(rf$importance)
    feaSet <- which(imp > 0.01)
    ruleSetPrunedRRF <- cbind(ruleMetric[feaSet,,drop=FALSE],impRRF=imp[feaSet])
    ix = order(as.numeric(ruleSetPrunedRRF[,"impRRF"]),
               - as.numeric(ruleSetPrunedRRF[,"err"]),
               - as.numeric(ruleSetPrunedRRF[,"len"]),
               decreasing=TRUE)
    ruleSelect <- ruleSetPrunedRRF[ix,,drop=FALSE]
    return(ruleSelect)
  }

treeVisit <-
  function(tree,rowIx,count,ruleSet,rule,levelX,length,max_length,digits=NULL)
  {
    if( tree[rowIx,"status"] == -1 | length == max_length ){
      count = count + 1
      ruleSet[[count]] = rule
      return(list(ruleSet = ruleSet, count=count))
    }
    xIx <- tree[rowIx,"split var"]
    xValue <- tree[rowIx,"split point"]
    if(is.integer(digits)) xValue <- round(tree[rowIx,"split point"], digits)

    if(is.null(levelX[[xIx]])){
      lValue <- paste("X[,",xIx, "]<=",xValue,sep="")
      rValue <- paste("X[,",xIx, "]>",xValue,sep="")
    }else{
      xValue<- which(as.integer(intToBits(as.integer(xValue)))>0)
      lValue <- levelX[[xIx]][xValue]
      rValue <- setdiff(levelX[[xIx]],lValue)
      #   lValue <- paste("X[,",xIx, "]%in% '",lValue,"'",sep="")
      #   rValue <- paste("X[,",xIx, "]%in% '",rValue,"'",sep="")
    }
    xValue <- NULL
    ruleleft <- rule
    if(length(ruleleft)==0)
    {
      ruleleft[[as.character(xIx)]] <- lValue
    }else{
      if(as.character(xIx) %in% ls(ruleleft)) {
        if(!is.null(levelX[[xIx]])){
          lValue <- intersect(ruleleft[[as.character(xIx)]],lValue)
          ruleleft[[as.character(xIx)]] <- lValue
        }else{
          ruleleft[[as.character(xIx)]] <- paste(ruleleft[[as.character(xIx)]], "&", lValue)
        }
      }else{
        ruleleft[[as.character(xIx)]] <- lValue
      }
    }

    #thisItem = paste("X[,",xIx, "] %in% ", nxValue, sep="")
    ruleright <- rule
    if(length(ruleright)==0)
    {
      ruleright[[as.character(xIx)]] <- rValue
    }else{
      if(as.character(xIx) %in% ls(ruleright)) {
        if(!is.null(levelX[[xIx]])){
          rValue <- intersect(ruleright[[as.character(xIx)]],rValue)
          ruleright[[as.character(xIx)]] <- rValue
        }else{
          ruleright[[as.character(xIx)]] <- paste(ruleright[[as.character(xIx)]], "&", rValue)
        }
      }else{
        ruleright[[as.character(xIx)]] <- rValue
      }
    }

    thisList = treeVisit(tree, tree[rowIx,"left daughter"],count,ruleSet,ruleleft,levelX,length+1,max_length,digits)
    ruleSet = thisList$ruleSet; count = thisList$count

    thisList = treeVisit(tree, tree[rowIx,"right daughter"],count,ruleSet,ruleright,levelX,length+1,max_length,digits)
    ruleSet = thisList$ruleSet; count = thisList$count

    return(list(ruleSet = ruleSet, count=count))
  }

ruleList2Exec <-
  function(X,allRulesList){
    typeX = getTypeX(X)
    ruleExec <- unique(t(sapply(allRulesList,singleRuleList2Exec,typeX=typeX)))
    ruleExec <- t(ruleExec)
    colnames(ruleExec) <- "condition"
    return(ruleExec)
  }

getTypeX <-
  function(X){
    typeX = rep(0,ncol(X))
    for(i in 1:ncol(X)){ #numeric: 1; categorical: 2s
      if(is.numeric(X[,i])){ typeX[i] = 1 }else{
        typeX[i] = 2
      }
    }
    return(typeX)
  }

singleRuleList2Exec <-
  function(ruleList,typeX){ #numeric: 1; categorical: 2s
    #ruleExec <- "which("
    ruleExec <- ""
    vars <- ls(ruleList)
    #ruleL <- length(unique(vars))
    vars <- vars[order(as.numeric(vars))]
    for(i in 1:length(vars)){
      if(typeX[as.numeric(vars[i])]==2){
        values <- paste("c(",paste(  paste("'",ruleList[[vars[i]]],"'",sep="")    ,collapse=","),")",sep="")
        tmp = paste("X[,",vars[i], "] %in% ", values, sep="")
      }else{
        tmp = ruleList[[vars[i]]]
      }
      if(i==1)ruleExec <- paste(ruleExec, tmp,sep="")
      if(i>1)ruleExec <- paste(ruleExec, " & ", tmp, sep="")
    }
    #ruleExec <- paste(ruleExec,")",sep="")
    return(c(ruleExec))
  }

measureRule <-
  function(ruleExec,X,target,pred=NULL,regMethod="mean"){
    len <- length(unlist(strsplit(ruleExec, split=" & ")))
    origRule <- ruleExec
    ruleExec <- paste("which(", ruleExec, ")")
    ixMatch <- eval(parse(text=ruleExec))
    if(length(ixMatch)==0){
      v <- c("-1","-1", "-1", "", "")
      names(v) <- c("len","freq","err","condition","pred")
      return(v)
    }
    ys <- target[ixMatch]
    freq <- round(length(ys)/nrow(X),digits=3)

    if(is.numeric(target))
    {
      if(regMethod == "median"){
        ysMost = stats::median(ys)
      }else{
        ysMost <- mean(ys)
      }
      err <- sum((ysMost - ys)^2)/length(ys)
    }else{
      if(length(pred)>0){ #if pred of the rule is provided
        ysMost = as.character(pred)
      }else{
        ysMost <- names(which.max(  table(ys))) # get back the first max
      }
      ly <- sum(as.character(ys)==ysMost)
      conf <- round(ly/length(ys),digits=3)
      err <- 1 - conf
    }
    rule <- origRule

    v <- c(len, freq, err, rule, ysMost)
    names(v) <- c("len","freq","err","condition","pred")
    return(v)
  }

pruneSingleRule <-
  function(rule, X, target, maxDecay, typeDecay){
    # typeDecay = 1: relative error increase; otherwise: absolute error increase

    #A <- gregexpr("X\\[,[0-9]+\\]", s)
    newRuleMetric <- measureRule(rule["condition"],X,target)
    errOrig <- as.numeric(newRuleMetric["err"])
    ruleV <- unlist(strsplit(rule["condition"],split= " & "))
    pred <- rule["pred"]
    # newRule <- NULL
    if(length(ruleV)==1) return(newRuleMetric)
    for(i in length(ruleV):1){
      restRule <- ruleV[-i]
      restRule <- paste(restRule,collapse= " & ")
      metricTmp <- measureRule(restRule,X,target,pred)
      errNew <- as.numeric(metricTmp["err"])
      if(typeDecay == 1){
        decay <- (errNew-errOrig)/max(errOrig,0.000001)
      }else{
        decay <- (errNew-errOrig)
      }
      if( decay <= maxDecay){
        #if( errNew-errOrig <= maxDecay){
        ruleV <- ruleV[-i]
        # newRule saves the last changed rule and metrics
        newRuleMetric <- metricTmp
        if(length(ruleV)<=1)break
      }
    }
    return(newRuleMetric)
    #rule["condition"] <- paste(ruleV,collapse= " & ")
    #return(rule)
  }

rule2Table <-
  function(ruleExec,X,target){
    I <- rep(0,nrow(X))
    ruleExec <- paste("which(", ruleExec, ")")
    ixMatch <- eval(parse(text=ruleExec))
    if(length(ixMatch)>0) I[ixMatch] <- 1
    names(I) = NULL
    return(I)
  }
